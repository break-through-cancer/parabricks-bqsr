#!/usr/bin/env python3

import pandas as pd

COLUMNS = ["patient", "sample", "status", "lane", "fastq_1", "fastq_2"]
FORM_ONLY = ("genome_source", "genome_index", "dbsnp", "known_indels", "custom_intervals", "use_intervals", "intervals_mode")


def build_fastq_samplesheet(files: pd.DataFrame, samplesheet: pd.DataFrame, log) -> pd.DataFrame:
    """One row per FASTQ pair, as sarek_align's make_manifest: Cirro's ingest already pairs mates
    (sampleIndex), so rows are keyed by sampleIndex/sample/lane/dataset and read 1/2 become
    fastq_1/fastq_2. Nothing is dropped: a pair with two files for the same read, or a read 1 without
    its read 2 in a paired sample, is an error."""
    if "readType" in files.columns:
        files = files.loc[files["readType"].fillna("R") == "R"]
    if files.empty:
        raise ValueError("No FASTQ files found in the input dataset(s)")
    if files["read"].isna().any():
        raise ValueError("Read number missing for: " + ", ".join(files.loc[files["read"].isna(), "file"]))

    f = files.assign(read=files["read"].astype(int), lane=files.get("lane", pd.Series(index=files.index, dtype=object)).fillna("1").astype(str))
    key = [c for c in ("sampleIndex", "sample", "lane", "dataset") if c in f.columns]
    paired = set(f.loc[f["read"] == 2, "sample"])

    rows, problems = [], []
    for _, group in f.groupby(key, sort=True, dropna=False):
        sample, lane = group["sample"].iloc[0], group["lane"].iloc[0]
        dup = group["read"].value_counts()
        for read in dup[dup > 1].index:
            problems.append(f"{sample}: more than one read {read} for one pair: {', '.join(group.loc[group['read'] == read, 'file'])}")
        mates = dict(zip(group["read"], group["file"]))
        if sample in paired and set(mates) != {1, 2}:
            problems.append(f"{sample}: no mate found for {', '.join(mates.values())}")
            continue
        rows.append(dict(sample=sample, lane=lane, fastq_1=mates.get(1, ""), fastq_2=mates.get(2, "")))
    if problems:
        raise ValueError("Cannot build FASTQ pairs:\n  " + "\n  ".join(problems))

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
