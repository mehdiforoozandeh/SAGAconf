"""Byte-identity tests: this SAGAconf against the legacy code (tag v0-legacy) on synthetic data.

    pytest tests/

The legacy code comes from `git archive v0-legacy`, or from SAGACONF_LEGACY=/path/to/checkout.
Each case in tests/golden/cases.json must give byte-identical files, stdout and stderr; a
crashing case may differ only in the traceback's file/line frames (see tests/golden/gate.py).
The same gate runs on the paper's data with tests/golden/cases_real.json.
"""
import json
import os
import shutil
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN = os.path.join(REPO, "tests", "golden")
CASES = os.path.join(GOLDEN, "cases.json")


@pytest.fixture(scope="session")
def legacy(tmp_path_factory):
    if os.environ.get("SAGACONF_LEGACY"):
        return os.environ["SAGACONF_LEGACY"]
    dest = tmp_path_factory.mktemp("legacy")
    archive = subprocess.run(["git", "-C", REPO, "archive", "v0-legacy"],
                             capture_output=True, check=True).stdout
    subprocess.run(["tar", "-x", "-C", str(dest)], input=archive, check=True)
    return str(dest)


@pytest.fixture(scope="session")
def data(tmp_path_factory, legacy):
    dest = tmp_path_factory.mktemp("synthetic")
    subprocess.run([sys.executable, os.path.join(GOLDEN, "make_synthetic.py"), str(dest)], check=True)
    os.makedirs(dest / "parsed")
    # The parsed inputs come from the legacy parser, so every run case starts from its output.
    for rep in ("base", "verif"):
        for fmt in ("bed", "csv"):
            work = dest / ("tmp_" + rep)
            work.mkdir(exist_ok=True)
            subprocess.run([sys.executable, os.path.join(legacy, "SAGAconf_parser.py"), "--saga", "chmm",
                            "--out_format", fmt, str(dest / "chmm" / rep), "200", str(work)],
                           check=True, capture_output=True)
            shutil.move(str(work / ("parsed_posterior." + fmt)), str(dest / "parsed" / (rep + "." + fmt)))
    return str(dest)


@pytest.mark.parametrize("name", [c["name"] for c in json.load(open(CASES))])
def test_byte_identical(name, legacy, data, tmp_path):
    r = subprocess.run(
        [sys.executable, os.path.join(GOLDEN, "gate.py"), "--legacy", legacy, "--new", REPO,
         "--cases", CASES, "--data", data, "--workdir", str(tmp_path), "--only", name,
         "--allow-traceback-frames"],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
