import os
import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]


def flat_config(*extra):
    out = subprocess.run(["nextflow", *extra, "config", "-flat"], cwd=ROOT, capture_output=True, text=True,
                         env={**os.environ, "PW_ONDEMAND_JOB_QUEUE": "q"}, check=True)
    return out.stdout


def parabricks_processes():
    names = []
    for nf in (ROOT / "modules" / "local").glob("*.nf"):
        text = nf.read_text()
        if "clara-parabricks" in text:
            names += re.findall(r"^process\s+(\w+)", text, re.M)
    return sorted(names)


def test_every_parabricks_process_requests_a_gpu():
    procs = parabricks_processes()
    assert procs == ["PARABRICKS_APPLYBQSR", "PARABRICKS_FQ2BAM", "PARABRICKS_FQ2BAM_PART", "PARABRICKS_MARKDUP"]
    for cfg in (flat_config(), flat_config("-c", ".cirro/align/process-compute.config")):
        for p in procs:
            assert re.search(rf"withName:{p}'\.accelerator = ", cfg), f"{p} has no accelerator"


def test_applybqsr_requests_at_most_two_gpus():
    for cfg in (flat_config(), flat_config("-c", ".cirro/align/process-compute.config"),
                flat_config("-c", ".cirro/apply_bqsr/process-compute.config")):
        m = re.search(r"withName:PARABRICKS_APPLYBQSR'\.accelerator = (\d+)", cfg)
        assert m and 1 <= int(m.group(1)) <= 2
