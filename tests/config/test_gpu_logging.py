import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]


def test_parabricks_tasks_log_gpu_model_and_memory_before_pbrun():
    # GPU type and memory decide whether --low-memory is needed and make A/B timings comparable.
    for nf in sorted((ROOT / "modules" / "local").glob("*.nf")):
        text = nf.read_text()
        if "clara-parabricks" not in text:
            continue
        script = text.split("script:", 1)[1].split("stub:", 1)[0]
        gpu = script.find("nvidia-smi --query-gpu=index,name,memory.total")
        pbrun = script.find("pbrun ")
        assert gpu != -1, f"{nf.name} does not log the GPU"
        assert gpu < pbrun, f"{nf.name} logs the GPU after pbrun"
