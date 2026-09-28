"""Write a small, seeded synthetic data set for a fast local run of the byte-identity gate.

    python tests/golden/make_synthetic.py /path/to/outdir

Layout written:
    chmm/{base,verif}/     ChromHMM-style *_chrN_posterior.txt files (chr21, chr22)
    segway/{base,verif}/   Segway-style posterior{i}.bedGraph files (integer percent)
    parsed/{base,verif}.bed and .csv   parsed posteriors, as SAGAconf_parser.py writes them
    mnemonics/{base,verif}.txt
The real-data gate uses the same cases on posteriors from the paper's runs.
"""
import os
import sys

import numpy as np

K = 6
RES = 200
CHROMS = {"chr21": 6000, "chr22": 4000}  # bins per chromosome


def states_and_posteriors(rng, n_bins, flip_from=None):
    if flip_from is None:
        seg = rng.integers(5, 40, size=n_bins)
        st = np.repeat(rng.integers(0, K, size=n_bins), seg)[:n_bins]
    else:
        st = flip_from.copy()
        flip = rng.random(n_bins) < 0.2
        st[flip] = rng.integers(0, K, flip.sum())
    p = np.full((n_bins, K), 0.2 / (K - 1))
    p[np.arange(n_bins), st] = 0.8
    p += rng.uniform(0, 0.05, size=p.shape)
    return st, p / p.sum(1, keepdims=True)


def main(out):
    rng = np.random.default_rng(20240321)
    post = {"base": {}, "verif": {}}
    for chrom, n in CHROMS.items():
        st, post["base"][chrom] = states_and_posteriors(rng, n)
        _, post["verif"][chrom] = states_and_posteriors(rng, n, flip_from=st)

    for rep, cell in (("base", "CellA"), ("verif", "CellB")):
        d = os.path.join(out, "chmm", rep)
        os.makedirs(d, exist_ok=True)
        for chrom, p in post[rep].items():
            with open(os.path.join(d, "%s_%d_%s_posterior.txt" % (cell, K, chrom)), "w") as f:
                f.write("%s\t%s\n" % (cell, chrom))
                f.write("\t".join("E%d" % (i + 1) for i in range(K)) + "\n")
                np.savetxt(f, p, fmt="%.4f", delimiter="\t")

        # Segway writes one bedGraph per label, run-length encoded, integer percent.
        d = os.path.join(out, "segway", rep)
        os.makedirs(d, exist_ok=True)
        for k in range(K):
            with open(os.path.join(d, "posterior%d.bedGraph" % k), "w") as f:
                f.write("track type=bedGraph name=posterior%d\n" % k)
                for chrom, p in post[rep].items():
                    pct = np.rint(p[:, k] * 100).astype(int)
                    start = 0
                    while start < len(pct):
                        end = start
                        while end + 1 < len(pct) and pct[end + 1] == pct[start] and end + 1 - start < 7:
                            end += 1
                        f.write("%s\t%d\t%d\t%d\n" % (chrom, start * RES, (end + 1) * RES, pct[start]))
                        start = end + 1

        d = os.path.join(out, "mnemonics")
        os.makedirs(d, exist_ok=True)
        names = ["Quie", "Prom", "Enha", "Tran", "Facu", "Enha_low"]
        with open(os.path.join(d, rep + ".txt"), "w") as f:
            f.write("old\tnew\n")
            order = rng.permutation(K)
            for i in range(K):
                f.write("%d\t%s\n" % (i, names[order[i]]))


if __name__ == "__main__":
    main(sys.argv[1])
