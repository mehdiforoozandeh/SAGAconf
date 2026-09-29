"""Byte-identity gate: run legacy SAGAconf and the new SAGAconf on the same inputs and
compare every file they write and everything they print.

Usage:
    python tests/golden/gate.py --legacy /path/to/legacy/checkout --new /path/to/new/checkout \
        --cases tests/golden/cases.json --data /path/to/data --workdir /path/to/scratch

Each case runs in a fresh, empty working directory, so files written relative to the
working directory are compared too. A case passes only if both runs exit with the same
code, write the same set of files, every file is byte-identical, stdout and stderr are
byte-identical, and every file listed in the case's "expect" field exists. The one exception
is Python traceback frames ('File "...", line N, in f' and the source line under it): they
name file paths and line numbers, which move whenever code moves. They count as a difference
unless --allow-traceback-frames is given; the exception type and message must always match.
A case's optional "stage" maps paths inside the working directory to input files that are
symlinked there before the run (for code that reads cwd-relative paths); staged paths are
not compared. A case with a "fixed_by" note is one a bug fix changes on purpose: it must
differ from legacy, but only by added files, a crash that no longer happens, or extra printed
lines; no file that legacy wrote may change. With --self-check, legacy is also run twice to prove the
comparison is deterministic.
"""
import argparse
import filecmp
import hashlib
import json
import os
import platform
import re
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
# On AVX-512 CPUs, numpy 1.26 np.exp on a strided float64 view (e.g. bins[:, 5] in
# posterior_calibration) returns last-bit-different results depending on the array's memory
# alignment, so legacy is not repeatable run to run. Disabling the AVX-512 dispatch makes it
# repeatable; both legacy and new run under the same setting.
if platform.machine() in ("x86_64", "AMD64"):
    DETERMINISTIC_ENV["NPY_DISABLE_CPU_FEATURES"] = (
        "AVX512F AVX512CD AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL AVX512_SPR")


def run_case(checkout, script, args, rundir, logdir, tag, stage):
    os.makedirs(rundir)
    for d in case_dirs_to_create(args):
        os.makedirs(os.path.join(rundir, d), exist_ok=True)
    for rel, src in stage.items():
        os.makedirs(os.path.dirname(os.path.join(rundir, rel)), exist_ok=True)
        os.symlink(src, os.path.join(rundir, rel))
    env = dict(os.environ, **DETERMINISTIC_ENV)
    env.pop("PYTHONPATH", None)
    cmd = [sys.executable, RUNNER, os.path.join(checkout, script)] + args
    logs = (os.path.join(logdir, tag + ".stdout"), os.path.join(logdir, tag + ".stderr"))
    t0 = time.time()
    with open(logs[0], "w") as out, open(logs[1], "w") as err:
        rc = subprocess.call(cmd, cwd=rundir, env=env, stdout=out, stderr=err)
    return rc, time.time() - t0, logs


TRACEBACK_FRAME = re.compile(r'^  File "[^"]*", line \d+, in .*$')


def strip_traceback_frames(lines):
    out, skip_source = [], False
    for line in lines:
        if TRACEBACK_FRAME.match(line):
            skip_source = True
            continue
        if skip_source and line.startswith("    "):
            skip_source = False
            continue
        skip_source = False
        out.append(line)
    return out


def compare_streams(ref_logs, cand_logs):
    problems = []
    for name, a, b in zip(("stdout", "stderr"), ref_logs, cand_logs):
        if filecmp.cmp(a, b, shallow=False):
            continue
        la, lb = open(a).read().splitlines(), open(b).read().splitlines()
        if strip_traceback_frames(la) == strip_traceback_frames(lb):
            problems.append("TRACEBACK_FRAMES_ONLY %s differs only in traceback frames" % name)
            continue
        n = next((i for i, (x, y) in enumerate(zip(la, lb)) if x != y), min(len(la), len(lb)))
        problems.append("%s differs at line %d: legacy=%r new=%r" % (
            name, n + 1, la[n] if n < len(la) else "<end>", lb[n] if n < len(lb) else "<end>"))
    return problems


def case_dirs_to_create(args):
    # The parser writes into an existing directory; mark such args with a trailing "/".
    return [a.rstrip("/") for a in args if a.endswith("/") and not os.path.isabs(a)]


def list_files(root, skip=()):
    out = []
    for dirpath, _, files in os.walk(root):
        for f in files:
            rel = os.path.relpath(os.path.join(dirpath, f), root)
            if rel not in skip:
                out.append(rel)
    return sorted(out)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compare_trees(a, b, skip_a=(), skip_b=()):
    fa, fb = list_files(a, skip_a), list_files(b, skip_b)
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
    p.add_argument("--allow-traceback-frames", action="store_true",
                   help="accept stderr that differs only in traceback file/line frames")
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
        stage = {k: v.replace("{data}", data) for k, v in c.get("stage", {}).items()}
        new_stage = {k: v.replace("{data}", data) for k, v in c.get("new_stage", c.get("stage", {})).items()}
        base = os.path.join(a.workdir, c["name"])
        if os.path.exists(base):
            shutil.rmtree(base)
        ref, cand = os.path.join(base, "legacy"), os.path.join(base, "new")
        rc_ref, t_ref, log_ref = run_case(a.legacy, c["script"], args, ref, logdir,
                                          c["name"] + ".legacy", stage)
        rc_new, t_new, log_new = run_case(a.new, new_script, new_args, cand, logdir,
                                          c["name"] + ".new", new_stage)

        files, problems = compare_trees(ref, cand, stage, new_stage)
        problems += compare_streams(log_ref, log_new)
        if rc_ref != rc_new:
            problems.insert(0, "exit codes differ: legacy=%d new=%d" % (rc_ref, rc_new))
        missing = [e for e in c.get("expect", []) if not os.path.exists(os.path.join(ref, e))]
        problems += ["expected file missing from legacy output: " + e for e in missing]
        if not files:
            problems.append("legacy wrote no files")
        if a.self_check:
            ref2 = os.path.join(base, "legacy_again")
            _, _, log_ref2 = run_case(a.legacy, c["script"], args, ref2, logdir,
                                      c["name"] + ".legacy_again", stage)
            problems += ["legacy is not deterministic: " + x for x in
                         compare_trees(ref, ref2, stage, stage)[1] + compare_streams(log_ref, log_ref2)]
        if a.allow_traceback_frames:
            problems = [x for x in problems if "TRACEBACK_FRAMES_ONLY" not in x]
        if c.get("fixed_by"):
            # A bug fix changes this case on purpose. It must change something, and it may
            # only add files, stop a crash (exit code, stderr, leftover temporary copies) or
            # print more; every file legacy wrote must keep its exact bytes.
            allowed = ("only in candidate: ", "exit codes differ", "stdout differs", "stderr differs",
                       "TRACEBACK_FRAMES_ONLY", "only in reference: out/base_replicate/",
                       "only in reference: out/verification_replicate/")
            if not problems:
                problems = ["fixed_by case is still identical to legacy: " + c["fixed_by"]]
            else:
                problems = [x for x in problems if not x.startswith(allowed)]

        ok = not problems
        failed += not ok
        report.append({
            "case": c["name"], "ok": ok, "n_files": len(files),
            "exit_legacy": rc_ref, "exit_new": rc_new,
            "seconds_legacy": round(t_ref, 1), "seconds_new": round(t_new, 1),
            "problems": problems[:50], "n_problems": len(problems),
            "sha256": {f: sha256(os.path.join(ref, f)) for f in files},
            "sha256_stdout": sha256(log_ref[0]), "sha256_stderr": sha256(log_ref[1]),
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
