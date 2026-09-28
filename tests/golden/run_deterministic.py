# Runs a SAGAconf entry script with deterministic plot output (fixed SVG ids, fixed PDF dates).
import os, sys, runpy
os.environ.setdefault("SOURCE_DATE_EPOCH", "0")
import matplotlib
matplotlib.use("Agg")
for rc in (matplotlib.rcParams, matplotlib.rcParamsDefault, matplotlib.rcParamsOrig):
    dict.__setitem__(rc, "svg.hashsalt", "sagaconf-golden")
script = sys.argv[1]
sys.argv = sys.argv[1:]
sys.path.insert(0, os.path.dirname(os.path.abspath(script)))
runpy.run_path(script, run_name="__main__")
