"""Regression tests for bugs fixed after 1.1.0 whose fix changes behaviour on purpose,
so the byte-identity gate (tests/test_golden.py) cannot cover them.

    pytest tests/test_fixes.py
"""
import os
import subprocess
import sys

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
