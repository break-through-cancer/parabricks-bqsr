import logging

import pandas as pd
import pytest

import importlib.util
import pathlib

_spec = importlib.util.spec_from_file_location("fastq_preprocess", pathlib.Path(__file__).with_name("preprocess.py"))
_preprocess = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_preprocess)
apply_genome_params = _preprocess.apply_genome_params
build_fastq_samplesheet = _preprocess.build_fastq_samplesheet

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


def test_top_up_runs_in_the_same_lane_are_all_kept():
    f = files([
        ("S1", "1", 1, "s3://b/run1/S1_L001_R1_001.fastq.gz"), ("S1", "1", 2, "s3://b/run1/S1_L001_R2_001.fastq.gz"),
        ("S1", "1", 1, "s3://b/run2/S1_L001_R1_001.fastq.gz"), ("S1", "1", 2, "s3://b/run2/S1_L001_R2_001.fastq.gz"),
    ])
    sheet = build_fastq_samplesheet(f, pd.DataFrame([{"sample": "S1"}]), LOG)
    assert sorted(sheet["fastq_1"]) == ["s3://b/run1/S1_L001_R1_001.fastq.gz", "s3://b/run2/S1_L001_R1_001.fastq.gz"]
    assert sorted(sheet["fastq_2"]) == ["s3://b/run1/S1_L001_R2_001.fastq.gz", "s3://b/run2/S1_L001_R2_001.fastq.gz"]
    assert sheet["lane"].tolist() == ["0", "1"]


def test_laneless_chunks_are_all_kept_and_paired_by_name():
    f = files([
        ("S1", None, 1, "s3://b/S1_part1_R1.fastq.gz"), ("S1", None, 2, "s3://b/S1_part1_R2.fastq.gz"),
        ("S1", None, 1, "s3://b/S1_part2_R1.fastq.gz"), ("S1", None, 2, "s3://b/S1_part2_R2.fastq.gz"),
    ])
    sheet = build_fastq_samplesheet(f, pd.DataFrame([{"sample": "S1"}]), LOG)
    pairs = sorted(zip(sheet["fastq_1"], sheet["fastq_2"]))
    assert pairs == [("s3://b/S1_part1_R1.fastq.gz", "s3://b/S1_part1_R2.fastq.gz"),
                     ("s3://b/S1_part2_R1.fastq.gz", "s3://b/S1_part2_R2.fastq.gz")]


def test_read_1_without_a_matching_read_2_in_a_paired_sample_is_an_error():
    f = files([
        ("S1", "1", 1, "s3://b/S1_A_R1.fastq.gz"), ("S1", "1", 2, "s3://b/S1_A_R2.fastq.gz"),
        ("S1", "1", 1, "s3://b/S1_B_R1.fastq.gz"),
    ])
    with pytest.raises(ValueError, match="S1_B_R1.fastq.gz"):
        build_fastq_samplesheet(f, pd.DataFrame([{"sample": "S1"}]), LOG)


def test_every_references_library_field_requests_an_s3_path():
    import json
    form = json.loads(pathlib.Path(__file__).with_name("process-form.json").read_text())

    def fields(node):
        if isinstance(node, dict):
            if node.get("pathType") == "references":
                yield node
            for value in node.values():
                yield from fields(value)
        elif isinstance(node, list):
            for value in node:
                yield from fields(value)

    refs = list(fields(form))
    assert len(refs) == 3
    assert all(f.get("useS3Path") is True for f in refs), [f["title"] for f in refs if not f.get("useS3Path")]


def test_every_dataset_field_is_filtered_by_process():
    import json
    form = json.loads(pathlib.Path(__file__).with_name("process-form.json").read_text())

    def fields(node, key=None):
        if isinstance(node, dict):
            if node.get("pathType") == "dataset":
                yield key, node
            for k, value in node.items():
                yield from fields(value, k)
        elif isinstance(node, list):
            for value in node:
                yield from fields(value, key)

    datasets = dict(fields(form))
    assert datasets, "expected at least one dataset field"
    assert all(f.get("process") for f in datasets.values()), [k for k, f in datasets.items() if not f.get("process")]
    assert datasets["genome_index"]["process"] == ["process-cirro-genome-index-bwa-1-0", "genome_bwa_index"]


def test_form_offers_low_memory_off_by_default_and_hides_gpuwrite():
    import json
    here = pathlib.Path(__file__).parent
    form = (here / "process-form.json").read_text()
    advanced = json.loads(form)["form"]["properties"]["advanced"]["properties"]
    assert advanced["fq2bam_low_memory"]["default"] is False
    assert "fq2bam_gpuwrite" not in form
    assert "fq2bam_gpuwrite" not in (here / "process-input.json").read_text()
