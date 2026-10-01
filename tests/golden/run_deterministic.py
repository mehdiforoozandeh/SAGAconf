# Runs a SAGAconf entry script with deterministic plot output (fixed SVG ids, fixed PDF dates).
import os, sys, runpy
os.environ.setdefault("SOURCE_DATE_EPOCH", "0")
import matplotlib
matplotlib.use("Agg")
for rc in (matplotlib.rcParams, matplotlib.rcParamsDefault, matplotlib.rcParamsOrig):
    dict.__setitem__(rc, "svg.hashsalt", "sagaconf-golden")
script = os.path.abspath(sys.argv[1])
sys.argv = sys.argv[1:]
if os.path.basename(script) == "__main__.py":
    # A package's __main__.py: run it as `python -m package` so relative imports work.
    pkg_dir = os.path.dirname(script)
    sys.path.insert(0, os.path.dirname(pkg_dir))
    runpy.run_module(os.path.basename(pkg_dir), run_name="__main__", alter_sys=True)
else:
    sys.path.insert(0, os.path.dirname(script))
    runpy.run_path(script, run_name="__main__")
