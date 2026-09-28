"""Byte-identity gate: run legacy SAGAconf and the new SAGAconf on the same inputs and
compare every file they write.

Usage:
    python tests/golden/gate.py --legacy /path/to/legacy/checkout --new /path/to/new/checkout \
        --cases tests/golden/cases.json --data /path/to/data --workdir /path/to/scratch

Each case runs in a fresh, empty working directory, so files written relative to the
working directory are compared too. A case passes only if both runs exit with the same
code, write the same set of files, every file is byte-identical, and every file listed in
the case's "expect" field exists. With --self-check, legacy is also run twice to prove the
comparison is deterministic.
"""
import argparse
import filecmp
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
RUNNER = os.path.join(HERE, "run_deterministic.py")

# Pin every source of run-to-run variation that is outside the code under test.
DETERMINISTIC_ENV = {
    "PYTHONHASHSEED": "0",
    "SOURCE_DATE_EPOCH": "0",
    "MPLBACKEND": "Agg",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
}


def run_case(checkout, script, args, rundir, logdir, tag):
    os.makedirs(rundir)
    for d in case_dirs_to_create(args):
        os.makedirs(os.path.join(rundir, d), exist_ok=True)
    env = dict(os.environ, **DETERMINISTIC_ENV)
    env.pop("PYTHONPATH", None)
    cmd = [sys.executable, RUNNER, os.path.join(checkout, script)] + args
    t0 = time.time()
    with open(os.path.join(logdir, tag + ".log"), "w") as log:
        rc = subprocess.call(cmd, cwd=rundir, env=env, stdout=log, stderr=subprocess.STDOUT)
    return rc, time.time() - t0


def case_dirs_to_create(args):
    # The parser writes into an existing directory; mark such args with a trailing "/".
    return [a.rstrip("/") for a in args if a.endswith("/") and not os.path.isabs(a)]


def list_files(root):
    out = []
    for dirpath, _, files in os.walk(root):
        for f in files:
            out.append(os.path.relpath(os.path.join(dirpath, f), root))
    return sorted(out)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compare_trees(a, b):
    fa, fb = list_files(a), list_files(b)
    problems = []
    only_a = sorted(set(fa) - set(fb))
    only_b = sorted(set(fb) - set(fa))
    problems += ["only in reference: " + f for f in only_a]
    problems += ["only in candidate: " + f for f in only_b]
    for f in sorted(set(fa) & set(fb)):
        if not filecmp.cmp(os.path.join(a, f), os.path.join(b, f), shallow=False):
            problems.append("bytes differ: " + f)
    return fa, problems


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--legacy", required=True, help="checkout of the legacy code (reference)")
    p.add_argument("--new", required=True, help="checkout of the new code (candidate)")
    p.add_argument("--cases", required=True, help="JSON list of cases")
    p.add_argument("--data", required=True, help="directory that {data} expands to")
    p.add_argument("--workdir", required=True, help="empty scratch directory for runs")
    p.add_argument("--only", nargs="*", help="run only these case names")
    p.add_argument("--self-check", action="store_true", help="also run legacy twice")
    a = p.parse_args()

    data = os.path.abspath(a.data)
    with open(a.cases) as fh:
        cases = json.load(fh)
    if a.only:
        cases = [c for c in cases if c["name"] in a.only]
    os.makedirs(a.workdir, exist_ok=True)
    logdir = os.path.join(a.workdir, "logs")
    os.makedirs(logdir, exist_ok=True)

    report, failed = [], 0
    for c in cases:
        args = [x.replace("{data}", data) for x in c["args"]]
        new_args = [x.replace("{data}", data) for x in c.get("new_args", c["args"])]
        new_script = c.get("new_script", c["script"])
        base = os.path.join(a.workdir, c["name"])
        if os.path.exists(base):
            shutil.rmtree(base)
        ref, cand = os.path.join(base, "legacy"), os.path.join(base, "new")
        rc_ref, t_ref = run_case(a.legacy, c["script"], args, ref, logdir, c["name"] + ".legacy")
        rc_new, t_new = run_case(a.new, new_script, new_args, cand, logdir, c["name"] + ".new")

        files, problems = compare_trees(ref, cand)
        if rc_ref != rc_new:
            problems.insert(0, "exit codes differ: legacy=%d new=%d" % (rc_ref, rc_new))
        missing = [e for e in c.get("expect", []) if not os.path.exists(os.path.join(ref, e))]
        problems += ["expected file missing from legacy output: " + e for e in missing]
        if not files:
            problems.append("legacy wrote no files")
        if a.self_check:
            ref2 = os.path.join(base, "legacy_again")
            run_case(a.legacy, c["script"], args, ref2, logdir, c["name"] + ".legacy_again")
            problems += ["legacy is not deterministic: " + x for x in compare_trees(ref, ref2)[1]]

        ok = not problems
        failed += not ok
        report.append({
            "case": c["name"], "ok": ok, "n_files": len(files),
            "exit_legacy": rc_ref, "exit_new": rc_new,
            "seconds_legacy": round(t_ref, 1), "seconds_new": round(t_new, 1),
            "problems": problems[:50], "n_problems": len(problems),
            "sha256": {f: sha256(os.path.join(ref, f)) for f in files},
        })
        print("%s %-28s files=%-4d legacy=%.0fs new=%.0fs%s" % (
            "PASS" if ok else "FAIL", c["name"], len(files), t_ref, t_new,
            "" if ok else "  (%d problems, first: %s)" % (len(problems), problems[0])), flush=True)

    with open(os.path.join(a.workdir, "gate_report.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    print("%d/%d cases passed; report: %s" % (
        len(cases) - failed, len(cases), os.path.join(a.workdir, "gate_report.json")))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
