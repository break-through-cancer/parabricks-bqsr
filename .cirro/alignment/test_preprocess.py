import pytest

from preprocess import build_samplesheet

ROOT = "s3://bucket/datasets/ds1/data/preprocessing"


def files(*names):
    return [f"{ROOT}/{n}" for n in names]


def test_one_sample_bam():
    sheet = build_samplesheet(files(
        "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai", "parabricks/S1/S1.table",
    ))
    assert list(sheet.columns) == ["sample", "alignment", "alignment_index", "recal_table"]
    assert sheet.to_dict("records") == [{
        "sample": "S1",
        "alignment": f"{ROOT}/parabricks/S1/S1.bam",
        "alignment_index": f"{ROOT}/parabricks/S1/S1.bam.bai",
        "recal_table": f"{ROOT}/parabricks/S1/S1.table",
    }]


def test_cram_and_multiple_samples_sorted():
    sheet = build_samplesheet(files(
        "parabricks/S2/S2.cram", "parabricks/S2/S2.cram.crai", "parabricks/S2/S2.table",
        "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai", "parabricks/S1/S1.table",
    ))
    assert list(sheet["sample"]) == ["S1", "S2"]
    assert sheet.loc[1, "alignment_index"].endswith("S2.cram.crai")


def test_other_stages_are_ignored():
    sheet = build_samplesheet(files(
        "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai", "parabricks/S1/S1.table",
        "recalibrated/S1/S1.recal.bam", "recalibrated/S1/S1.recal.bam.bai",
        "markduplicates/S1/S1.md.cram", "markduplicates/S1/S1.md.cram.crai",
    ))
    assert sheet["alignment"].tolist() == [f"{ROOT}/parabricks/S1/S1.bam"]


def test_no_parabricks_files():
    with pytest.raises(ValueError, match="preprocessing/parabricks"):
        build_samplesheet(files("recalibrated/S1/S1.recal.bam", "recalibrated/S1/S1.recal.bam.bai"))


def test_missing_table_and_index_are_reported_per_sample():
    with pytest.raises(ValueError) as err:
        build_samplesheet(files(
            "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai",
            "parabricks/S2/S2.cram", "parabricks/S2/S2.table",
        ))
    msg = str(err.value)
    assert "S1: no .table recalibration table" in msg
    assert "S2: no index for" in msg


def test_two_alignments_for_one_sample_is_an_error():
    with pytest.raises(ValueError, match="S1: expected one BAM/CRAM, found 2"):
        build_samplesheet(files(
            "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai",
            "parabricks/S1/S1.cram", "parabricks/S1/S1.cram.crai",
            "parabricks/S1/S1.table",
        ))


def test_same_sample_in_two_input_datasets_is_an_error():
    other = "s3://bucket/datasets/ds2/data/preprocessing"
    with pytest.raises(ValueError, match="S1: expected one BAM/CRAM, found 2"):
        build_samplesheet(files(
            "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai", "parabricks/S1/S1.table",
        ) + [f"{other}/parabricks/S1/S1.bam", f"{other}/parabricks/S1/S1.bam.bai", f"{other}/parabricks/S1/S1.table"])
