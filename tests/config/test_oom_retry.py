import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]


def test_parabricks_memory_grows_with_each_attempt():
    # A ~60x sample was killed for host memory at 88 GB; a retry with the same memory cannot succeed.
    for name in ("parabricks_fq2bam.nf", "parabricks_applybqsr.nf"):
        text = (ROOT / "modules" / "local" / name).read_text()
        memory = [l for l in text.splitlines() if l.strip().startswith("memory")][0]
        assert "task.attempt" in memory, f"{name}: {memory.strip()}"


def test_cirro_retries_parabricks_when_no_exit_code_is_reported():
    # AWS Batch reports a container memory kill without an exit code; Nextflow records Integer.MAX_VALUE.
    for proc in ("fastq", "alignment"):
        text = (ROOT / ".cirro" / proc / "process-compute.config").read_text()
        block = re.search(r"withName: '[^']*PARABRICKS_APPLYBQSR[^']*' \{(.*?)\n    \}", text, re.S).group(1)
        assert "2147483647" in block, proc
