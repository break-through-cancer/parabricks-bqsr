import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]
MODULE = (ROOT / "modules" / "local" / "mosdepth.nf").read_text()


def script(sample="S1", alignment="S1.cram"):
    body = MODULE.split("script:", 1)[1].split('"""', 2)[1]
    for var, value in {"task.cpus": "4", "ref_fasta": "ref.fa", "meta.sample": sample, "alignment": alignment}.items():
        body = body.replace("${" + var + "}", value)
    return body.replace("\\$", "$").replace("\\\\", "\\")


def run(tmp_path, fake):
    exe = tmp_path / "bin" / "mosdepth"
    exe.parent.mkdir()
    exe.write_text("#!/bin/bash\n" + fake)
    exe.chmod(0o755)
    return subprocess.run(["bash", "-c", script()], cwd=tmp_path, capture_output=True, text=True,
                          env={"PATH": f"{exe.parent}:/usr/bin:/bin"})


def test_mosdepth_killed_for_memory_keeps_its_exit_code_so_the_task_is_retried(tmp_path):
    # A 1M-read CRAM was killed at 4 GB; the old guard turned 137 into 1, which no retry rule matches.
    r = run(tmp_path, "echo 'Killed' >&2\nexit 137\n")
    assert r.returncode == 137


def test_decode_errors_still_fail_with_a_message(tmp_path):
    fake = ("echo '[E::cram_decode_slice] Slice decode failure' >&2\n"
            "printf 'chrom\\tlength\\tbases\\tmean\\tmin\\tmax\\ntotal\\t100\\t0\\t0\\t0\\t0\\n' > S1.mosdepth.summary.txt\n")
    r = run(tmp_path, fake)
    assert r.returncode == 1
    assert "could not read S1.cram" in r.stderr


def test_successful_mosdepth_passes(tmp_path):
    fake = "printf 'chrom\\tlength\\tbases\\tmean\\tmin\\tmax\\ntotal\\t100\\t500\\t5\\t0\\t9\\n' > S1.mosdepth.summary.txt\n"
    r = run(tmp_path, fake)
    assert r.returncode == 0, r.stderr


def test_mosdepth_memory_grows_with_each_attempt():
    memory = [l for l in MODULE.splitlines() if l.strip().startswith("memory")][0]
    assert "task.attempt" in memory, memory.strip()
