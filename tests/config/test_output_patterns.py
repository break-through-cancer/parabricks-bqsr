import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]


def test_no_output_uses_an_either_or_index_pattern():
    # AWS Batch unstaging expands output globs with `ls`, so `{bai,crai}` logs
    # "cannot access" for whichever index does not exist.
    offenders = [f"{p.name}:{i}" for p in sorted((ROOT / "modules" / "local").glob("*.nf"))
                 for i, line in enumerate(p.read_text().splitlines(), 1) if "{bai,crai}" in line]
    assert offenders == []
