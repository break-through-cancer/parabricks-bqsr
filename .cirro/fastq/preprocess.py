#!/usr/bin/env python3

import re

import pandas as pd

COLUMNS = ["patient", "sample", "status", "lane", "fastq_1", "fastq_2"]
FORM_ONLY = ("genome_source", "genome_index", "dbsnp", "known_indels", "custom_intervals", "use_intervals", "intervals_mode")


READ_TOKEN = re.compile(r"(?<=[._])(R?)([12])(?=[._])")


def pair_key(path: str) -> str:
    """The path with its last read-number token (R1/R2, _1/_2) masked, so mates share a key."""
    head, _, name = path.rpartition("/")
    matches = list(READ_TOKEN.finditer(name))
    if matches:
        m = matches[-1]
        name = f"{name[:m.start()]}{m.group(1)}#{name[m.end():]}"
    return f"{head}/{name}"


def build_fastq_samplesheet(files: pd.DataFrame, samplesheet: pd.DataFrame, log) -> pd.DataFrame:
    """One row per FASTQ (pair): read 1/2 become fastq_1/fastq_2; patient/status from sample metadata.

    Every file is kept: mates are paired by path with the read number masked, so top-up runs and
    lane-less chunks each become their own row. An unmatched read in a paired sample is an error.
    """
    if "readType" in files.columns:
        files = files.loc[files["readType"].fillna("R") == "R"]
    if files.empty:
        raise ValueError("No FASTQ files found in the input dataset(s)")
    if files["read"].isna().any():
        raise ValueError("Read number missing for: " + ", ".join(files.loc[files["read"].isna(), "file"]))

    rows, problems = [], []
    for sample, group in files.groupby("sample", sort=True):
        paired = (group["read"].astype(int) == 2).any()
        by_key = {}
        for rec in group.to_dict("records"):
            lane = "1" if pd.isna(rec.get("lane")) else str(rec["lane"])
            by_key.setdefault((lane, pair_key(rec["file"])), {})[int(rec["read"])] = rec["file"]
        for (lane, _key), mates in sorted(by_key.items()):
            if paired and set(mates) != {1, 2}:
                problems.append(f"{sample}: no mate found for {', '.join(mates.values())}")
                continue
            rows.append(dict(sample=sample, lane=lane, fastq_1=mates[1], fastq_2=mates.get(2, "")))
    if problems:
        raise ValueError("Cannot pair FASTQ files:\n  " + "\n  ".join(problems))

    wide = pd.DataFrame(rows)
    meta = samplesheet.reindex(columns=["sample", "patient", "status"]).set_index("sample")
    missing = meta["status"].isna().sum() + len(set(wide["sample"]) - set(meta.index))
    if missing:
        log.warning(f"status not provided for {missing} sample(s), defaulting to 0 (normal)")
    wide["patient"] = wide["sample"].map(meta["patient"]).fillna(wide["sample"])
    wide["status"] = wide["sample"].map(meta["status"]).fillna(0).astype(int)
    wide = wide.sort_values(["sample", "lane", "fastq_1"]).reset_index(drop=True)
    wide["lane"] = [str(i) for i in range(len(wide))]
    return wide[COLUMNS]


def apply_genome_params(params: dict) -> dict:
    """Translate Cirro form fields into pipeline genome parameters."""
    p = dict(params)
    source = p.get("genome_source", "igenomes")
    if source == "dataset":
        index = p.get("genome_index")
        if not index:
            raise ValueError("Custom genome selected but no BWA genome index dataset was provided")
        p["genome"] = "null"
        p["bwa_index"] = index
        p["ref_fasta"] = f"{index}/genome.fasta"
        p["ref_fasta_fai"] = f"{index}/genome.fasta.fai"
        sites = [p.get(k) for k in ("dbsnp", "known_indels") if p.get(k)]
        if sites:
            p["known_sites"] = ",".join(sites)
        if p.get("custom_intervals"):
            p["intervals"] = p["custom_intervals"]
    else:
        mode = p.get("intervals_mode") or ("none" if p.get("use_intervals") is False else "gatk_calling_regions")
        if mode == "none":
            p["no_intervals"] = True
        elif mode == "custom":
            if not p.get("custom_intervals"):
                raise ValueError("Intervals set to custom intervals but no intervals BED was selected")
            p["intervals"] = p["custom_intervals"]
    for k in FORM_ONLY:
        p.pop(k, None)
    return p


if __name__ == "__main__":
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()
    sheet = build_fastq_samplesheet(ds.files, ds.samplesheet, ds.logger)
    ds.logger.info(sheet.to_csv(index=False))
    sheet.to_csv(ds.params["input"], index=False)
    resolved = apply_genome_params(ds.params)
    for key in FORM_ONLY:
        ds.remove_param(key, force=True)
    for key in ("genome", "bwa_index", "ref_fasta", "ref_fasta_fai", "known_sites", "intervals", "no_intervals"):
        if key in resolved:
            ds.add_param(key, resolved[key], overwrite=True)
