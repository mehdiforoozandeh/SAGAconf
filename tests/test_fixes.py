"""Regression tests for bugs fixed after 1.1.0 whose fix changes behaviour on purpose,
so the byte-identity gate (tests/test_golden.py) cannot cover them.

    pytest tests/test_fixes.py
"""
import os
import subprocess
import sys

from sagaconf._utils import mp_inplace_binning

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HEADER = "chr\tstart\tend\tposterior0\tposterior1\n"


def _parsed(path, value):
    rows = ["chr21\t%d\t%d\t%s\t%s\n" % (s, s + 200, value, value) for s in range(0, 2000, 200)]
    with open(path, "w") as f:
        f.write(HEADER + "".join(rows))
    return str(path)


def test_full_mode_exits_nonzero_when_a_step_fails(tmp_path):
    # A verification file that cannot be read as numbers makes every full-mode step fail.
    base = _parsed(tmp_path / "base.bed", "0.5")
    verif = _parsed(tmp_path / "verif.bed", "x")
    for cmd in ([sys.executable, "-m", "sagaconf", "run"], [sys.executable, os.path.join(REPO, "SAGAconf.py")]):
        out = tmp_path / ("out_%d" % len(cmd))
        r = subprocess.run(cmd + [base, verif, str(out)], capture_output=True, text=True, cwd=REPO)
        assert r.returncode == 1, r.stderr
        assert "analysis step(s) failed" in r.stderr
        assert "r-values" in r.stderr


def _segway_dir(tmp_path, files):
    for i, text in enumerate(files):
        (tmp_path / ("posterior%d.bedGraph" % i)).write_text("track type=bedGraph\n" + text)
    return str(tmp_path)


def test_segway_last_line_without_newline(tmp_path):
    d = _segway_dir(tmp_path, ["chr1\t0\t100\t63\nchr1\t100\t200\t37\n",
                               "chr1\t0\t100\t37\nchr1\t100\t200\t63"])
    parsed = mp_inplace_binning(d, 100)
    assert list(parsed["posterior1"]) == [0.37, 0.63]
