import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
FASTQ_GPU = ("PARABRICKS_FQ2BAM", "PARABRICKS_FQ2BAM_PART", "PARABRICKS_MARKDUP", "PARABRICKS_APPLYBQSR")


def gpu_block(proc, name):
    text = (ROOT / ".cirro" / proc / "process-compute.config").read_text()
    for selector, body in re.findall(r"withName: '([^']*)' \{(.*?)\n    \}", text, re.S):
        if re.fullmatch(selector, name):
            return body
    raise AssertionError(f"{proc}: no selector matches {name}")


def test_parabricks_memory_grows_with_each_attempt():
    # A ~60x sample was killed for host memory at 88 GB; a retry with the same memory cannot succeed.
    for name in ("parabricks_fq2bam.nf", "parabricks_fq2bam_part.nf", "parabricks_markdup.nf", "parabricks_applybqsr.nf"):
        text = (ROOT / "modules" / "local" / name).read_text()
        memory = [l for l in text.splitlines() if l.strip().startswith("memory")][0]
        assert "task.attempt" in memory, f"{name}: {memory.strip()}"


def test_cirro_runs_every_parabricks_process_on_demand_and_retries_memory_kills():
    # AWS Batch reports a container memory kill without an exit code; Nextflow records Integer.MAX_VALUE.
    for proc, names in (("fastq", FASTQ_GPU), ("alignment", ("PARABRICKS_APPLYBQSR",))):
        for name in names:
            block = gpu_block(proc, name)
            assert "PW_ONDEMAND_JOB_QUEUE" in block, (proc, name)
            assert "2147483647" in block, (proc, name)


def test_fq2bam_memory_override_reaches_both_fq2bam_processes():
    for name in ("parabricks_fq2bam.nf", "parabricks_fq2bam_part.nf"):
        text = (ROOT / "modules" / "local" / name).read_text()
        memory = [l for l in text.splitlines() if l.strip().startswith("memory")][0]
        assert "params.fq2bam_memory_gb" in memory, f"{name}: {memory.strip()}"
