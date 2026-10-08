import pathlib
import re
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
MODULES = ROOT / "modules" / "local"
KILLED = ("Process terminated with signal [SIGKILL: 9]. SIGKILL cannot be caught. A common reason for SIGKILL is running out of\n"
          "host memory.\nCould not run fq2bam\nExiting pbrun ...")


def wrapper():
    text = (MODULES / "pbrun_wrapper.nf").read_text()
    return re.search(r"'''\n(.*?)'''", text, re.S).group(1)


def run(tmp_path, output, code):
    fake = tmp_path / "bin" / "pbrun"
    fake.parent.mkdir()
    fake.write_text(f"#!/bin/bash\ncat <<'EOF'\n{output}\nEOF\nexit {code}\n")
    fake.chmod(0o755)
    script = f"set -euo pipefail\n{wrapper()}\npbrun fq2bam --ref x\necho after\n"
    return subprocess.run(["bash", "-c", script], cwd=tmp_path, capture_output=True, text=True,
                          env={"PATH": f"{fake.parent}:/usr/bin:/bin"})


def test_pbrun_killed_for_memory_exits_137_so_the_task_is_retried(tmp_path):
    # pbrun survives the kernel killing its BWA child and exits 255, which no retry rule matches.
    r = run(tmp_path, KILLED, 255)
    assert r.returncode == 137
    assert "SIGKILL" in r.stdout
    assert "after" not in r.stdout


def test_other_pbrun_errors_keep_their_exit_code(tmp_path):
    r = run(tmp_path, "[Parabricks Options Error]: Please specifiy one input from --in-fq, --in-fq-list, --in-se-fq, or --in-se-fq-list", 255)
    assert r.returncode == 255


def test_successful_pbrun_continues_the_script(tmp_path):
    r = run(tmp_path, "done", 0)
    assert r.returncode == 0
    assert "after" in r.stdout


@pytest.mark.parametrize("module", sorted(p.name for p in MODULES.glob("*.nf") if "clara-parabricks" in p.read_text()))
def test_every_parabricks_module_runs_pbrun_through_the_wrapper(module):
    text = (MODULES / module).read_text()
    script = text.split("script:", 1)[1].split("stub:", 1)[0]
    assert "${pbrunFunction()}" in script, module
    assert "pbrunFunction" in text.split("process ", 1)[0], f"{module} does not include pbrunFunction"
