import json
import os
import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]


def config_params():
    out = subprocess.run(["nextflow", "config", "-flat"], cwd=ROOT, capture_output=True, text=True,
                         env={**os.environ}, check=True).stdout
    params = {}
    for line in out.splitlines():
        m = re.match(r"^params\.(\w+) = (.*)$", line)
        if m:
            raw = m.group(2)
            if raw == "null":
                value = None
            elif raw in ("true", "false"):
                value = raw == "true"
            elif re.fullmatch(r"-?\d+", raw):
                value = int(raw)
            else:
                value = raw.strip("'")
            params[m.group(1)] = value
    return params


def schema_params():
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    props = {}
    for group in schema.get("$defs", {}).values():
        props.update(group.get("properties", {}))
    props.update(schema.get("properties", {}))
    return props


def test_schema_declares_exactly_the_config_params():
    assert sorted(schema_params()) == sorted(config_params())


def test_schema_defaults_match_config():
    props = schema_params()
    for name, value in config_params().items():
        if value is None:
            assert "default" not in props[name], f"{name}: config default is null, schema sets {props[name].get('default')!r}"
        else:
            assert props[name].get("default") == value, f"{name}: config {value!r} != schema {props[name].get('default')!r}"


def test_pipeline_requires_nextflow_26_04():
    out = subprocess.run(["nextflow", "config", "-flat"], cwd=ROOT, capture_output=True, text=True,
                         env={**os.environ}, check=True).stdout
    assert "manifest.nextflowVersion = '!>=26.04.0'" in out
