#!/usr/bin/env python3

import re
from collections import defaultdict
from typing import Iterable

import pandas as pd

PARABRICKS_FILE = re.compile(r"(?:^|/)preprocessing/parabricks/(?P<sample>[^/]+)/(?P<name>[^/]+)$")
MARKDUPLICATES_FILE = re.compile(r"(?:^|/)preprocessing/markduplicates/(?P<sample>[^/]+)/(?P<name>[^/]+)$")
RECAL_TABLE_FILE = re.compile(r"(?:^|/)preprocessing/recal_table/(?P<sample>[^/]+)/(?P<name>[^/]+)$")
INDEX_SUFFIX = {"bam": "bai", "cram": "crai"}
COLUMNS = ["sample", "alignment", "alignment_index", "recal_table"]


def build_samplesheet(paths: Iterable[str]) -> pd.DataFrame:
    """Build the pipeline samplesheet from an upstream dataset's file paths.

    Accepts either sarek_align's layout (preprocessing/parabricks/<sample>/: the
    pre-BQSR fq2bam alignment, its index and the recalibration table, co-located) or
    this pipeline's own apply_bqsr=false output (preprocessing/markduplicates/<sample>/:
    alignment and index; preprocessing/recal_table/<sample>/: the table, separately).
    Alignments from other stages (for example recalibrated/) are never used, since
    applying the table again would recalibrate those reads twice.
    """
    paths = list(paths)
    present = set(paths)
    found = defaultdict(lambda: {"alignments": [], "tables": []})

    for path in paths:
        match = PARABRICKS_FILE.search(path) or MARKDUPLICATES_FILE.search(path)
        if match:
            name = match["name"]
            ext = name.rsplit(".", 1)[-1]
            if ext in INDEX_SUFFIX:
                found[match["sample"]]["alignments"].append(path)
            elif ext == "table":
                found[match["sample"]]["tables"].append(path)
            continue
        match = RECAL_TABLE_FILE.search(path)
        if match and match["name"].endswith(".table"):
            found[match["sample"]]["tables"].append(path)

    if not found:
        raise ValueError(
            "No files under preprocessing/parabricks/<sample>/ (sarek_align, Parabricks "
            "aligner, known sites supplied, save_mapped on, baserecalibrator skipped) or "
            "preprocessing/markduplicates/<sample>/ + preprocessing/recal_table/<sample>/ "
            "(this pipeline's own apply_bqsr=false output) in the input dataset(s). This "
            "pipeline needs the pre-BQSR alignment, its index, and the recalibration table."
        )

    rows, problems = [], []
    for sample in sorted(found):
        alignments, tables = found[sample]["alignments"], found[sample]["tables"]
        if len(alignments) != 1:
            problems.append(f"{sample}: expected one BAM/CRAM, found {len(alignments)}")
        if len(tables) != 1:
            problems.append(
                f"{sample}: no .table recalibration table" if not tables
                else f"{sample}: expected one .table, found {len(tables)}"
            )
        if len(alignments) != 1 or len(tables) != 1:
            continue
        alignment = alignments[0]
        index = f"{alignment}.{INDEX_SUFFIX[alignment.rsplit('.', 1)[-1]]}"
        if index not in present:
            problems.append(f"{sample}: no index for {alignment.rsplit('/', 1)[-1]}")
            continue
        rows.append(dict(sample=sample, alignment=alignment, alignment_index=index, recal_table=tables[0]))

    if problems:
        raise ValueError("Cannot build the samplesheet:\n  " + "\n  ".join(problems))
    return pd.DataFrame(rows, columns=COLUMNS)


if __name__ == "__main__":
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()
    samplesheet = build_samplesheet(ds.files["file"])
    ds.logger.info(samplesheet.to_csv(index=False))
    samplesheet.to_csv(ds.params["input"], index=False)
    ds.logger.info(f"Wrote {samplesheet.shape[0]} row(s) to {ds.params['input']}")
