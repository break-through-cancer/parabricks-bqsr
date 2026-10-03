#!/usr/bin/env python3

import pandas as pd

COLUMNS = ["patient", "sample", "status", "lane", "fastq_1", "fastq_2"]
FORM_ONLY = ("genome_source", "genome_index", "dbsnp", "known_indels", "custom_intervals", "use_intervals")


def build_fastq_samplesheet(files: pd.DataFrame, samplesheet: pd.DataFrame, log) -> pd.DataFrame:
    """One row per lane: read 1/2 become fastq_1/fastq_2; patient/status from sample metadata."""
    if "readType" in files.columns:
        files = files.loc[files["readType"].fillna("R") == "R"]
    if files.empty:
        raise ValueError("No FASTQ files found in the input dataset(s)")

    wide = (
        files.assign(read=files["read"].astype(int), lane=files["lane"].fillna("1").astype(str))
        .pivot_table(index=["sample", "lane"], columns="read", values="file", aggfunc="first")
        .rename(columns={1: "fastq_1", 2: "fastq_2"})
        .reset_index()
    )
    if "fastq_2" not in wide.columns:
        wide["fastq_2"] = ""
    wide["fastq_2"] = wide["fastq_2"].fillna("")

    meta = samplesheet.reindex(columns=["sample", "patient", "status"]).set_index("sample")
    missing = meta["status"].isna().sum() + len(set(wide["sample"]) - set(meta.index))
    if missing:
        log.warning(f"status not provided for {missing} sample(s), defaulting to 0 (normal)")
    wide["patient"] = wide["sample"].map(meta["patient"]).fillna(wide["sample"])
    wide["status"] = wide["sample"].map(meta["status"]).fillna(0).astype(int)
    wide = wide.sort_values(["sample", "lane"]).reset_index(drop=True)
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
    elif p.get("use_intervals") is False:
        p["no_intervals"] = True
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
