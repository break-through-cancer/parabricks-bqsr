import logging

import pandas as pd
import pytest

import importlib.util
import pathlib

_spec = importlib.util.spec_from_file_location("align_preprocess", pathlib.Path(__file__).with_name("preprocess.py"))
_preprocess = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_preprocess)
apply_genome_params = _preprocess.apply_genome_params
build_fastq_samplesheet = _preprocess.build_fastq_samplesheet

LOG = logging.getLogger("test")


def files(rows, with_index=True):
    cols = ["sampleIndex", "sample", "lane", "read", "file"]
    df = pd.DataFrame(rows, columns=cols)
    return df if with_index else df.drop(columns=["sampleIndex"])


META = pd.DataFrame([{"sample": "S1"}])


def test_pairs_by_sample_index_with_metadata_defaults():
    f = files([
        (0, "S1", "1", 1, "s3://b/a_R1.fq.gz"), (0, "S1", "1", 2, "s3://b/a_R2.fq.gz"),
        (1, "S1", "2", 1, "s3://b/b_R1.fq.gz"), (1, "S1", "2", 2, "s3://b/b_R2.fq.gz"),
        (2, "S2", "1", 1, "s3://b/c.fq.gz"),
    ])
    meta = pd.DataFrame([{"sample": "S1", "patient": "P1", "status": 1}, {"sample": "S2"}])
    sheet = build_fastq_samplesheet(f, meta, LOG)
    assert list(sheet.columns) == ["patient", "sample", "status", "lane", "fastq_1", "fastq_2"]
    assert sheet["lane"].tolist() == ["0", "1", "2"]
    s2 = sheet[sheet["sample"] == "S2"].iloc[0]
    assert s2["patient"] == "S2" and s2["status"] == 0 and s2["fastq_2"] == ""
    assert sheet[sheet["sample"] == "S1"]["status"].tolist() == [1, 1]


def test_pairing_does_not_depend_on_file_names():
    f = files([(0, "S1", "1", 1, "s3://b/x_1-merged.fastq.gz"), (0, "S1", "1", 2, "s3://b/x_2-merged.fastq.gz")])
    assert build_fastq_samplesheet(f, META, LOG)[["fastq_1", "fastq_2"]].values.tolist() == [["s3://b/x_1-merged.fastq.gz", "s3://b/x_2-merged.fastq.gz"]]


def test_top_up_runs_in_the_same_lane_are_all_kept():
    f = files([
        (0, "S1", "1", 1, "s3://b/run1/S1_R1.fq.gz"), (0, "S1", "1", 2, "s3://b/run1/S1_R2.fq.gz"),
        (1, "S1", "1", 1, "s3://b/run2/S1_R1.fq.gz"), (1, "S1", "1", 2, "s3://b/run2/S1_R2.fq.gz"),
    ])
    sheet = build_fastq_samplesheet(f, META, LOG)
    assert sorted(zip(sheet["fastq_1"], sheet["fastq_2"])) == [("s3://b/run1/S1_R1.fq.gz", "s3://b/run1/S1_R2.fq.gz"),
                                                               ("s3://b/run2/S1_R1.fq.gz", "s3://b/run2/S1_R2.fq.gz")]


def test_laneless_chunks_are_all_kept():
    f = files([(0, "S1", None, 1, "s3://b/p1_R1.fq.gz"), (0, "S1", None, 2, "s3://b/p1_R2.fq.gz"),
               (1, "S1", None, 1, "s3://b/p2_R1.fq.gz"), (1, "S1", None, 2, "s3://b/p2_R2.fq.gz")])
    assert len(build_fastq_samplesheet(f, META, LOG)) == 2


def test_read_1_without_read_2_in_a_paired_sample_is_an_error():
    f = files([(0, "S1", "1", 1, "s3://b/a_R1.fq.gz"), (0, "S1", "1", 2, "s3://b/a_R2.fq.gz"), (1, "S1", "1", 1, "s3://b/b_R1.fq.gz")])
    with pytest.raises(ValueError, match="b_R1.fq.gz"):
        build_fastq_samplesheet(f, META, LOG)


def test_two_files_for_the_same_read_of_one_pair_is_an_error_not_a_silent_drop():
    f = files([(0, "S1", "1", 1, "s3://b/a_R1.fq.gz"), (0, "S1", "1", 1, "s3://b/b_R1.fq.gz"), (0, "S1", "1", 2, "s3://b/a_R2.fq.gz")])
    with pytest.raises(ValueError, match="more than one read 1"):
        build_fastq_samplesheet(f, META, LOG)


def test_without_sample_index_one_pair_per_lane_still_works_and_extra_files_are_an_error():
    ok = files([("x", "S1", "1", 1, "s3://b/a_R1.fq.gz"), ("x", "S1", "1", 2, "s3://b/a_R2.fq.gz")], with_index=False)
    assert len(build_fastq_samplesheet(ok, META, LOG)) == 1
    bad = files([("x", "S1", "1", 1, "s3://b/a_R1.fq.gz"), ("x", "S1", "1", 1, "s3://b/b_R1.fq.gz"),
                 ("x", "S1", "1", 2, "s3://b/a_R2.fq.gz"), ("x", "S1", "1", 2, "s3://b/b_R2.fq.gz")], with_index=False)
    with pytest.raises(ValueError, match="more than one read 1"):
        build_fastq_samplesheet(bad, META, LOG)


def test_index_reads_are_ignored():
    f = files([(0, "S1", "1", 1, "s3://b/a_R1.fq.gz"), (0, "S1", "1", 2, "s3://b/a_R2.fq.gz"), (0, "S1", "1", 1, "s3://b/a_I1.fq.gz")])
    f["readType"] = ["R", "R", "I"]
    assert len(build_fastq_samplesheet(f, META, LOG)) == 1


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
    assert len(refs) == 4
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


def test_igenomes_with_custom_intervals_uses_them():
    p = apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38",
                             "intervals_mode": "custom", "custom_intervals": "s3://r/c/regions.bed"})
    assert p["intervals"] == "s3://r/c/regions.bed" and "no_intervals" not in p
    assert "intervals_mode" not in p and "custom_intervals" not in p


def test_igenomes_intervals_none_sets_no_intervals():
    p = apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38", "intervals_mode": "none"})
    assert p["no_intervals"] is True and "intervals" not in p


def test_igenomes_default_intervals_leaves_the_genome_default():
    p = apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38", "intervals_mode": "gatk_calling_regions"})
    assert "intervals" not in p and "no_intervals" not in p


def test_igenomes_custom_intervals_without_a_file_is_an_error():
    with pytest.raises(ValueError, match="custom intervals"):
        apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38", "intervals_mode": "custom"})


def test_gpu_field_defaults_to_2_and_caps_at_4():
    import json
    adv = json.loads((pathlib.Path(__file__).parent / "process-form.json").read_text())["form"]["properties"]["advanced"]["properties"]
    assert adv["fq2bam_gpus"]["default"] == 2 and adv["fq2bam_gpus"]["maximum"] == 4


def test_igenomes_branch_offers_three_interval_modes_and_a_custom_bed():
    import json
    form = json.loads((pathlib.Path(__file__).parent / "process-form.json").read_text())
    igenomes = form["form"]["properties"]["genome_selection"]["dependencies"]["genome_source"]["oneOf"][0]
    mode = igenomes["properties"]["intervals_mode"]
    assert mode["enum"] == ["gatk_calling_regions", "none", "custom"] and mode["default"] == "gatk_calling_regions"
    custom = [b for b in igenomes["dependencies"]["intervals_mode"]["oneOf"] if b["properties"]["intervals_mode"]["enum"] == ["custom"]][0]
    bed = custom["properties"]["custom_intervals"]
    assert bed["pathType"] == "references" and bed["useS3Path"] is True and custom["required"] == ["custom_intervals"]


def test_singleton_found_in_the_file_list_becomes_a_single_end_row_and_is_not_paired():
    f = files([(0, "S1", "1", 1, "s3://b/S1_R1.fq.gz"), (0, "S1", "1", 2, "s3://b/S1_R2.fq.gz"),
               (1, "S1", "1", 1, "s3://b/S1__singletons.fastq.gz")])
    meta = pd.DataFrame([{"sample": "S1", "fastq_singleton": "S1__singletons.fastq.gz"}])
    sheet = build_fastq_samplesheet(f, meta, LOG)
    assert sorted(zip(sheet["fastq_1"], sheet["fastq_2"])) == [("s3://b/S1_R1.fq.gz", "s3://b/S1_R2.fq.gz"), ("s3://b/S1__singletons.fastq.gz", "")]


def test_singleton_missing_from_the_file_list_is_resolved_next_to_fastq_1():
    f = files([(0, "S1", "1", 1, "s3://b/d/S1_R1.fq.gz"), (0, "S1", "1", 2, "s3://b/d/S1_R2.fq.gz")])
    meta = pd.DataFrame([{"sample": "S1", "fastq_singleton": "S1__singletons.fastq.gz"}])
    assert "s3://b/d/S1__singletons.fastq.gz" in build_fastq_samplesheet(f, meta, LOG)["fastq_1"].tolist()


def test_singleton_with_an_unset_read_number_is_still_handled():
    f = files([(0, "S1", "1", 1, "s3://b/S1_R1.fq.gz"), (0, "S1", "1", 2, "s3://b/S1_R2.fq.gz"),
               (1, "S1", "1", None, "s3://b/S1__singletons.fastq.gz")])
    meta = pd.DataFrame([{"sample": "S1", "fastq_singleton": "S1__singletons.fastq.gz"}])
    assert len(build_fastq_samplesheet(f, meta, LOG)) == 2


def test_samples_without_a_singleton_are_unchanged():
    f = files([(0, "S1", "1", 1, "s3://b/S1_R1.fq.gz"), (0, "S1", "1", 2, "s3://b/S1_R2.fq.gz"),
               (1, "S2", "1", 1, "s3://b/S2_R1.fq.gz"), (1, "S2", "1", 2, "s3://b/S2_R2.fq.gz")])
    meta = pd.DataFrame([{"sample": "S1", "fastq_singleton": "S1__singletons.fastq.gz"}, {"sample": "S2", "fastq_singleton": None}])
    sheet = build_fastq_samplesheet(f, meta, LOG)
    assert len(sheet[sheet["sample"] == "S2"]) == 1 and len(sheet[sheet["sample"] == "S1"]) == 2


def test_singleton_for_a_sample_without_fastqs_is_an_error():
    f = files([(0, "S1", "1", 1, "s3://b/S1_R1.fq.gz"), (0, "S1", "1", 2, "s3://b/S1_R2.fq.gz")])
    meta = pd.DataFrame([{"sample": "S9", "fastq_singleton": "S9__singletons.fastq.gz"}])
    with pytest.raises(ValueError, match="S9"):
        build_fastq_samplesheet(f, meta, LOG)


def test_singleton_named_in_cirros_file_list_column_is_resolved_next_to_its_mates():
    d = "s3://project-p/datasets/src/data"
    f = pd.DataFrame([
        {"sample": "S1", "file": f"{d}/S1_R1.fastq.gz", "dataset": "src", "sampleIndex": 1, "read": 1, "readType": "R", "singleton": "S1_singletons.fastq.gz"},
        {"sample": "S1", "file": f"{d}/S1_R2.fastq.gz", "dataset": "src", "sampleIndex": 1, "read": 2, "readType": "R", "singleton": "S1_singletons.fastq.gz"},
    ])
    sheet = build_fastq_samplesheet(f, pd.DataFrame([{"sample": "S1"}]), LOG)
    assert sorted(zip(sheet["fastq_1"], sheet["fastq_2"])) == [(f"{d}/S1_R1.fastq.gz", f"{d}/S1_R2.fastq.gz"), (f"{d}/S1_singletons.fastq.gz", "")]


def test_form_does_not_expose_the_fq2bam_intermediate_format():
    here = pathlib.Path(__file__).parent
    assert "fq2bam_intermediate_fmt" not in (here / "process-form.json").read_text()
    assert "fq2bam_intermediate_fmt" not in (here / "process-input.json").read_text()


def test_form_offers_an_optional_fq2bam_memory_override_and_maps_it():
    import json
    here = pathlib.Path(__file__).parent
    advanced = json.loads((here / "process-form.json").read_text())["form"]["properties"]["advanced"]
    field = advanced["properties"]["fq2bam_memory_gb"]
    assert field["type"] == "integer" and field["minimum"] == 16 and field["maximum"] == 768
    assert "default" not in field and "fq2bam_memory_gb" not in advanced.get("required", [])
    mapping = json.loads((here / "process-input.json").read_text())
    assert mapping["fq2bam_memory_gb"] == "$.dataset.params.advanced.fq2bam_memory_gb"


def test_form_offers_poly_g_trimming_auto_by_default_and_maps_it():
    import json
    here = pathlib.Path(__file__).parent
    trimming = json.loads((here / "process-form.json").read_text())["form"]["properties"]["read_trimming"]["properties"]
    assert trimming["poly_g_trimming"]["default"] == "auto"
    assert trimming["poly_g_trimming"]["enum"] == ["auto", "on", "off"]
    assert "trim_nextseq" not in trimming
    mapping = json.loads((here / "process-input.json").read_text())
    assert mapping["poly_g_trimming"] == "$.dataset.params.read_trimming.poly_g_trimming"
    assert "trim_nextseq" not in mapping


def test_form_offers_publish_markduplicates_off_by_default_and_maps_it():
    import json
    here = pathlib.Path(__file__).parent
    props = json.loads((here / "process-form.json").read_text())["form"]["properties"]
    field = props["publish_markduplicates"]
    assert field["type"] == "boolean" and field["default"] is False
    mapping = json.loads((here / "process-input.json").read_text())
    assert mapping["publish_markduplicates"] == "$.dataset.params.publish_markduplicates"


def test_form_offers_apply_bqsr_on_by_default_and_maps_it():
    import json
    here = pathlib.Path(__file__).parent
    props = json.loads((here / "process-form.json").read_text())["form"]["properties"]
    field = props["apply_bqsr"]
    assert field["type"] == "boolean" and field["default"] is True
    mapping = json.loads((here / "process-input.json").read_text())
    assert mapping["apply_bqsr"] == "$.dataset.params.apply_bqsr"


def test_quantization_section_only_offered_when_apply_bqsr_is_on():
    import json
    form = json.loads((pathlib.Path(__file__).parent / "process-form.json").read_text())["form"]
    assert "quantization" not in form["properties"]
    branches = form["dependencies"]["apply_bqsr"]["oneOf"]
    on_branch = [b for b in branches if b["properties"]["apply_bqsr"]["enum"] == [True]][0]
    off_branch = [b for b in branches if b["properties"]["apply_bqsr"]["enum"] == [False]][0]
    assert on_branch["properties"]["quantization"]["properties"]["quantize_quals_enabled"]["default"] is True
    assert "quantization" not in off_branch["properties"]
