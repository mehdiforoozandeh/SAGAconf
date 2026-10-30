"""Check a Segway parse made with the fix for GitHub issue #4 against the reference parse.

    python tests/golden/verify_segway_binning.py POSTERIOR_DIR RESOLUTION REFERENCE CANDIDATE

The legacy parser split each interval wider than RESOLUTION into sub-bins and wrote int(p)
for each, so every probability strictly between 0 and 1 became 0. The fix writes p. This
script rebuilds, from the bedGraph files, which cell of the parsed file comes from such a
wide interval, then requires: the same header (column order), the same number of rows, the
same coordinate fields, the same text in every cell from an interval no wider than
RESOLUTION, and, in every cell from a wider interval, the candidate value equal to the
interval's probability (value/100) and the reference value equal to int(value/100).
The gate (gate.py) calls check() for cases with a "segway_binning" entry.
"""
import os
import sys


def interval_cells(posterior_dir, resolution):
    """For each posterior{i}.bedGraph, in the parser's order, a list of (value, wide) per row."""
    n = len([f for f in os.listdir(posterior_dir) if ".bedGraph" in f])
    cols = []
    for i in range(n):
        rows = []
        with open(os.path.join(posterior_dir, "posterior%d.bedGraph" % i)) as fh:
            next(fh)  # track line
            for line in fh:
                _, start, end, value = line.split("\t")
                start, end, value = int(start), int(end), int(value)
                if end - start > resolution:
                    rows += [(value, True)] * len(range(start, end, resolution))
                else:
                    rows.append((value, False))
        cols.append(rows)
    return cols


def check(posterior_dir, resolution, ref_file, cand_file):
    """Return a list of problems; empty means the candidate differs only as the fix intends."""
    sep, n_coord = ("\t", 3) if ref_file.endswith(".bed") else (",", 4)  # csv has an index column
    with open(ref_file) as fh:
        ref = fh.read().splitlines()
    with open(cand_file) as fh:
        cand = fh.read().splitlines()
    if ref[0] != cand[0]:
        return ["header differs: reference=%r candidate=%r" % (ref[0], cand[0])]
    if len(ref) != len(cand):
        return ["row count differs: reference=%d candidate=%d" % (len(ref) - 1, len(cand) - 1)]
    cols = interval_cells(posterior_dir, resolution)
    if len(ref[0].split(sep)) != n_coord + len(cols):
        return ["%d posterior files but header has %d value columns" % (len(cols), len(ref[0].split(sep)) - n_coord)]
    if any(len(c) != len(ref) - 1 for c in cols):
        return ["posterior files do not all bin to %d rows" % (len(ref) - 1)]

    problems = []
    for r in range(1, len(ref)):
        a, b = ref[r].split(sep), cand[r].split(sep)
        if a[:n_coord] != b[:n_coord]:
            problems.append("row %d: coordinates differ: %r vs %r" % (r, a[:n_coord], b[:n_coord]))
            continue
        for k, (value, wide) in enumerate(c[r - 1] for c in cols):
            x, y = a[n_coord + k], b[n_coord + k]
            if not wide:
                if x != y:
                    problems.append("row %d col %d: narrow interval changed: %r -> %r" % (r, k, x, y))
                continue
            if float(y) != value / 100:
                problems.append("row %d col %d: wide interval value %r, expected %r" % (r, k, y, value / 100))
            if float(x) != int(value / 100):
                problems.append("row %d col %d: reference wide value %r, expected int(%r)" % (r, k, x, value / 100))
        if len(problems) >= 50:
            break
    return problems


if __name__ == "__main__":
    out = check(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4])
    print("\n".join(out) or "OK")
    sys.exit(1 if out else 0)
