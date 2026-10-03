import logging

import pandas as pd
import pytest

from preprocess import apply_genome_params, build_fastq_samplesheet

LOG = logging.getLogger("test")


def files(rows):
    return pd.DataFrame(rows, columns=["sample", "lane", "read", "file"])


def test_paired_and_single_end_with_metadata_defaults():
    f = files([
        ("S1", "1", 1, "s3://b/S1_L1_R1.fastq.gz"), ("S1", "1", 2, "s3://b/S1_L1_R2.fastq.gz"),
        ("S1", "2", 1, "s3://b/S1_L2_R1.fastq.gz"), ("S1", "2", 2, "s3://b/S1_L2_R2.fastq.gz"),
        ("S2", "1", 1, "s3://b/S2_L1.fastq.gz"),
    ])
    meta = pd.DataFrame([{"sample": "S1", "patient": "P1", "status": 1}, {"sample": "S2"}])
    sheet = build_fastq_samplesheet(f, meta, LOG)
    assert list(sheet.columns) == ["patient", "sample", "status", "lane", "fastq_1", "fastq_2"]
    assert sheet["lane"].tolist() == ["0", "1", "2"]
    s2 = sheet[sheet["sample"] == "S2"].iloc[0]
    assert s2["patient"] == "S2" and s2["status"] == 0 and s2["fastq_2"] == ""
    assert sheet[sheet["sample"] == "S1"]["status"].tolist() == [1, 1]


def test_missing_lane_is_kept():
    f = files([("S1", None, 1, "s3://b/S1_R1.fastq.gz"), ("S1", None, 2, "s3://b/S1_R2.fastq.gz")])
    assert len(build_fastq_samplesheet(f, pd.DataFrame([{"sample": "S1"}]), LOG)) == 1


def test_index_files_are_ignored():
    f = files([("S1", "1", 1, "s3://b/S1_R1.fastq.gz"), ("S1", "1", 2, "s3://b/S1_R2.fastq.gz")])
    f["readType"] = ["R", "R"]
    f.loc[len(f)] = ["S1", "1", 1, "s3://b/S1_I1.fastq.gz", "I"]
    assert len(build_fastq_samplesheet(f, pd.DataFrame([{"sample": "S1"}]), LOG)) == 1


def test_no_files_is_an_error():
    with pytest.raises(ValueError, match="No FASTQ files"):
        build_fastq_samplesheet(files([]), pd.DataFrame(columns=["sample"]), LOG)


def test_igenomes_with_intervals_unticked():
    p = apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38", "use_intervals": False,
                             "igenomes_base": "s3://refs/igenomes/"})
    assert p["no_intervals"] is True and "genome_source" not in p and "use_intervals" not in p
    assert p["genome"] == "GATK.GRCh38"


def test_custom_genome_resolves_dataset_files_and_known_sites():
    p = apply_genome_params({"genome_source": "dataset", "genome_index": "s3://b/ds/data",
                             "dbsnp": "s3://r/a/germline_resource.vcf.gz", "known_indels": "s3://r/b/germline_resource.vcf.gz",
                             "custom_intervals": "s3://r/c/regions.bed", "use_intervals": True})
    assert p["genome"] == "null"
    assert p["bwa_index"] == "s3://b/ds/data"
    assert p["ref_fasta"] == "s3://b/ds/data/genome.fasta" and p["ref_fasta_fai"] == "s3://b/ds/data/genome.fasta.fai"
    assert p["known_sites"] == "s3://r/a/germline_resource.vcf.gz,s3://r/b/germline_resource.vcf.gz"
    assert p["intervals"] == "s3://r/c/regions.bed"
    assert not {"genome_source", "genome_index", "dbsnp", "known_indels", "custom_intervals", "use_intervals"} & p.keys()


def test_custom_genome_without_index_is_an_error():
    with pytest.raises(ValueError, match="BWA genome index"):
        apply_genome_params({"genome_source": "dataset"})
