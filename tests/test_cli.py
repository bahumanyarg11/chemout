import json
from pathlib import Path

import pytest

from chemout.cli import main

DATA = Path(__file__).parent / "data" / "cclib-data"


def test_sniff(capsys):
    assert main(["sniff", str(DATA / "ORCA6.0" / "1618.out")]) == 0
    assert capsys.readouterr().out.strip() == "orca"


def test_info(capsys):
    assert main(["info", str(DATA / "Gaussian16" / "water_neutral_nbo_opt_freq.out")]) == 0
    out = capsys.readouterr().out
    assert "program" in out and "gaussian" in out
    assert "gibbs_hartree" in out
    assert "opt_converged" in out


def test_parse_json(tmp_path):
    out = tmp_path / "x.json"
    assert main(["parse", str(DATA / "xTB" / "1768_xyz.out"), "-o", str(out)]) == 0
    d = json.loads(out.read_text())
    assert d["program"]["name"] == "xtb"
    assert d["energies"]["total"] == pytest.approx(-39.629079572158)


def test_parse_jsonl_multi(tmp_path):
    out = tmp_path / "x.jsonl"
    files = [str(DATA / "xTB" / "1768_xyz.out"), str(DATA / "ORCA6.0" / "commercial_version.out")]
    assert main(["parse", *files, "-o", str(out), "-q"]) == 0
    lines = out.read_text().splitlines()
    assert len(lines) == 2
    assert json.loads(lines[1])["program"]["name"] == "orca"


def test_batch_csv_and_hdf5(tmp_path):
    pytest.importorskip("h5py")
    import h5py

    csv_path = tmp_path / "s.csv"
    h5 = tmp_path / "all.h5"
    rc = main(["batch", str(DATA / "xTB"), "--csv", str(csv_path), "-o", str(h5), "-q", "-j", "2"])
    assert rc == 0
    rows = csv_path.read_text().splitlines()
    assert rows[0].startswith("path,program,version,termination")
    assert len(rows) == 4  # header + 3 files
    with h5py.File(h5) as f:
        assert len(f.keys()) == 3
        g = f["1768_xyz"]
        assert g.attrs["program"] == "xtb"
        assert "geometries/coordinates" in g


def test_strict_exit_code(tmp_path):
    bad = tmp_path / "bad.out"
    bad.write_text("not a quantum chemistry file\n")
    assert main(["parse", str(bad), "--strict", "-q", "-o", str(tmp_path / "o.json")]) == 2
