"""Command-line interface.

    sagaconf parse --saga {chmm,segway} POSTERIOR_DIR RESOLUTION OUTDIR
    sagaconf run BASE VERIF OUTDIR [--mode MODE] [options]

`legacy_run_main` and `legacy_parse_main` keep the original SAGAconf.py and
SAGAconf_parser.py interfaces (flags, help text and error messages) for existing pipelines.
"""
import argparse
import os
import shutil
from importlib.metadata import PackageNotFoundError, version

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

from ._chromhmm import ChrHMM_read_posteriordir
from ._utils import mp_inplace_binning
from .overall import is_repr_posterior
from .reports import (get_all_ct, get_overalls, get_rvals_activeregion, load_data,
                      overlap_vs_segment_length, post_clustering_keep_k_states, process_data,
                      quick_report, subset_data_to_activeregions)

MODES = ("full", "quick", "rvalues", "celltype", "seglength", "merge", "active-regions")
# Legacy boolean flags, in the order the original script checked them.
LEGACY_MODE_FLAGS = (("quick", "quick"), ("r_only", "rvalues"), ("active_regions", "active-regions"),
                     ("v_seglength", "seglength"), ("ct_only", "celltype"), ("merge_only", "merge"))
DEFAULT_CCRE_FILE = "sagaconf/biointerpret/GRCh38-cCREs.bed"
DEFAULT_MEULEMAN_FILE = "sagaconf/biointerpret/Meuleman.tsv"


def parse_posteriors(posteriordir, resolution, savedir, saga, out_format="bed"):
    """Bin a SAGA model's posterior files and write savedir/parsed_posterior.{bed,csv}."""
    if saga == "segway":
        binned_posterior = mp_inplace_binning(posteriordir, resolution)
    elif saga == "chmm":
        binned_posterior = ChrHMM_read_posteriordir(posteriordir, resolution)

    if out_format == "bed":
        binned_posterior.to_csv(savedir + "/parsed_posterior.bed", sep='\t', header=True, index=False)
    elif out_format == "csv":
        binned_posterior.to_csv(savedir + "/parsed_posterior.csv")


def run_sagaconf(base, verif, savedir, mode="full", base_mnemonics="NA", verif_mnemonics="NA",
                 chr21_only=False, window_size=1000, iou_threshold=0.75, repr_threshold=0.8,
                 merge_k=-1, verbose=False, ccre_file=DEFAULT_CCRE_FILE,
                 meuleman_file=DEFAULT_MEULEMAN_FILE):
    """Compare a base and a verification annotation; write the results to savedir.

    base, verif: parsed posterior files (.bed or .csv), as written by `sagaconf parse`.
    base_mnemonics, verif_mnemonics: state-name files, or "NA" for none.
    """
    os.makedirs(savedir, exist_ok=True)
    # Work on copies inside savedir; they are removed at the end.
    replicate_1_dir = savedir + "/base_replicate"
    replicate_2_dir = savedir + "/verification_replicate"
    os.makedirs(replicate_1_dir, exist_ok=True)
    os.makedirs(replicate_2_dir, exist_ok=True)

    posterior1_dir = _copy_posterior(base, replicate_1_dir)
    posterior2_dir = _copy_posterior(verif, replicate_2_dir)
    if base_mnemonics != "NA":
        shutil.copyfile(base_mnemonics, replicate_1_dir + "/mnemonics.txt")
    if verif_mnemonics != "NA":
        shutil.copyfile(verif_mnemonics, replicate_2_dir + "/mnemonics.txt")
    mnem = base_mnemonics != "NA" and verif_mnemonics != "NA"

    w = window_size

    def load(logit_transform):
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=chr21_only,
                                 logit_transform=logit_transform, force_WG=False)
        return process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem,
                            bm=base_mnemonics, vm=verif_mnemonics, match=False, custom_order=True)

    def r_values(loci1, loci2):
        return is_repr_posterior(loci1, loci2, ovr_threshold=iou_threshold, window_bp=w,
                                 matching="static", always_include_best_match=True, return_r=True)

    if mode == "quick":
        loci1, loci2 = load(True)
        quick_report(loci1, loci2, savedir, locis=True, w=w, to=iou_threshold, tr=repr_threshold)

    elif mode == "rvalues":
        loci1, loci2 = load(True)
        r_values(loci1, loci2).to_csv(savedir + "/r_values.bed", sep='\t', header=True, index=False)

    elif mode == "active-regions":
        loci1, loci2 = load(False)
        regions = dict(cCREs_file=ccre_file, Meuleman_file=meuleman_file, locis=True)
        loci1, loci2 = subset_data_to_activeregions(
            replicate_1_dir=loci1, replicate_2_dir=loci2, restrict_to="WG", **regions)
        get_rvals_activeregion(loci1, loci2, savedir, w=w, restrict_to="WG")
        for restrict_to in ("cCRE", "muel"):
            sub1, sub2 = subset_data_to_activeregions(
                replicate_1_dir=loci1, replicate_2_dir=loci2, restrict_to=restrict_to, **regions)
            get_rvals_activeregion(sub1, sub2, savedir, w=w, restrict_to=restrict_to)

    elif mode == "seglength":
        loci1, loci2 = load(True)
        overlap_vs_segment_length(replicate_1_dir=loci1, replicate_2_dir=loci2, savedir=savedir, locis=True)

    elif mode == "celltype":
        loci1, loci2 = load(True)
        get_all_ct(loci1, loci2, savedir, locis=True, w=w)

    elif mode == "merge":
        loci1, loci2 = load(False)
        post_clustering_keep_k_states(loci1, loci2, savedir, k=merge_k, locis=True, write_csv=False, w=w)

    else:  # full
        # Each step is independent: a failure is reported with -v and the next step still runs.
        try:
            loci1, loci2 = load(True)
            get_all_ct(loci1, loci2, savedir, locis=True, w=w)
        except Exception:
            if verbose:
                print("Failed to generated sample analysis.")

        if merge_k != -1:
            try:
                loci1, loci2 = load(False)
                post_clustering_keep_k_states(loci1, loci2, savedir, k=merge_k, locis=True, write_csv=False)
            except Exception:
                if verbose:
                    print("Failed to merge clusters up to k")

        try:
            loci1, loci2 = load(True)
            get_overalls(loci1, loci2, savedir, locis=True, w=w, to=iou_threshold, tr=repr_threshold)
        except Exception:
            if verbose:
                print("failed to get GW SAGAconf reproducibility results")

        try:
            rvalues = r_values(loci1, loci2)
            rvalues.to_csv(savedir + "/r_values.bed", sep='\t', header=True, index=False)
            _rval_hist(rvalues, savedir)
        except Exception:
            if verbose:
                print("failed to get GW SAGAconf reproducibility results")

    shutil.rmtree(replicate_1_dir, ignore_errors=True)
    shutil.rmtree(replicate_2_dir, ignore_errors=True)


def _copy_posterior(path, replicate_dir):
    for ext in (".bed", ".csv"):
        if ext in path.lower():
            dest = "%s/parsed_posterior%s" % (replicate_dir, ext)
            shutil.copyfile(path, dest)
            return dest
    raise ValueError("posterior file must be .bed or .csv: %s" % path)


def _rval_hist(rvalues, savedir):
    labels = rvalues.MAP.unique()
    fig, axs = plt.subplots(len(labels), 1, figsize=(15, 15), sharex=True, sharey=False)

    bin_edges = list(np.arange(0, 1, 0.01))
    for l in range(len(labels)):
        data = rvalues.loc[rvalues["MAP"] == labels[l], "r_value"]
        weights = np.ones_like(data) / len(data)

        axs[l].hist(data, bins=bin_edges, color="black", alpha=0.6,
                    label=labels[l], weights=weights)

        axs[l].text(0.02, 0.95, labels[l], transform=axs[l].transAxes,
                    horizontalalignment='left', verticalalignment='top',
                    fontsize=8)
        axs[l].tick_params(axis='y', labelsize=8)
        axs[l].set_xlim([0, 1])

    plt.tight_layout()
    plt.savefig('{}/rval_hist.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/rval_hist.svg'.format(savedir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()


# ---- `sagaconf` command ----------------------------------------------------------------------

def _hidden_alias(parser, flag, dest, **kw):
    """Accept an old flag name without listing it in --help."""
    parser.add_argument(flag, dest=dest, help=argparse.SUPPRESS, default=argparse.SUPPRESS, **kw)


def build_parser():
    try:
        pkg_version = version("sagaconf")
    except PackageNotFoundError:
        pkg_version = "unknown"
    parser = argparse.ArgumentParser(
        prog="sagaconf",
        description="Calibrated reproducibility scores (r-values) for chromatin state annotations.")
    parser.add_argument("--version", action="version", version="%(prog)s " + pkg_version)
    sub = parser.add_subparsers(dest="command", required=True, metavar="{parse,run}")

    p = sub.add_parser(
        "parse", help="convert ChromHMM or Segway posteriors into a parsed posterior file",
        description="Bin a SAGA model's posterior files into one parsed posterior file "
                    "(chr, start, end, one column per state).")
    p.add_argument("posteriordir", help="directory with the SAGA model's posterior files")
    p.add_argument("resolution", type=int, help="bin size of the SAGA model, in bp")
    p.add_argument("outdir", help="directory for parsed_posterior.bed (or .csv); created if missing")
    p.add_argument("--saga", required=True, choices=["chmm", "segway"], help="which SAGA model wrote the posteriors")
    p.add_argument("--out-format", dest="out_format", choices=["bed", "csv"], default="bed",
                   help="output format [default: bed]")
    _hidden_alias(p, "--out_format", "out_format", choices=["bed", "csv"])

    r = sub.add_parser(
        "run", help="score the reproducibility of a base annotation against a verification annotation",
        description="Compare a base and a verification annotation and write r-values, confident "
                    "segments and reports to OUTDIR.")
    r.add_argument("base", help="parsed posterior file (.bed or .csv) of the base annotation")
    r.add_argument("verif", help="parsed posterior file (.bed or .csv) of the verification annotation")
    r.add_argument("outdir", help="directory for the results; created if missing")
    r.add_argument("--mode", choices=MODES, default=None,
                   help="full: all analyses (default); quick: essential report only; rvalues: r_values.bed "
                        "only; celltype: per-annotation analyses; seglength: reproducibility vs position in "
                        "segment; merge: merge base states down to -k states; active-regions: r-values in "
                        "cCREs (needs --ccre-file)")
    r.add_argument("-bm", "--base-mnemonics", dest="base_mnemonics", default="NA", metavar="FILE",
                   help="state names for the base annotation (tab-separated: old, new)")
    r.add_argument("-vm", "--verif-mnemonics", dest="verif_mnemonics", default="NA", metavar="FILE",
                   help="state names for the verification annotation")
    r.add_argument("-s", "--chr21-only", dest="chr21_only", action="store_true",
                   help="analyse chr21 only (a quick check)")
    r.add_argument("-w", "--window-size", dest="window_size", type=int, default=1000, metavar="BP",
                   help="window around each bin, in bp [default: 1000]")
    r.add_argument("-to", "--iou-threshold", dest="iou_threshold", type=float, default=0.75,
                   help="IoU overlap needed to call two states corresponding [default: 0.75]")
    r.add_argument("-tr", "--repr-threshold", dest="repr_threshold", type=float, default=0.8,
                   help="r-value needed to call a bin reproduced (alpha in the paper) [default: 0.8]")
    r.add_argument("-k", "--merge-k", dest="merge_k", type=int, default=-1, metavar="K",
                   help="merge base states down to K states")
    r.add_argument("--ccre-file", default=DEFAULT_CCRE_FILE, help="cCRE BED file for --mode active-regions")
    r.add_argument("--meuleman-file", default=DEFAULT_MEULEMAN_FILE, help="Meuleman et al. DHS index for --mode active-regions")
    r.add_argument("-v", "--verbose", action="store_true", help="report which analysis steps failed")
    # Old SAGAconf.py flag names keep working.
    for flag, dest in (("--verbosity", "verbose"), ("--base_mnemonics", "base_mnemonics"),
                       ("--verif_mnemonics", "verif_mnemonics"), ("--windowsize", "window_size"),
                       ("--iou_threshold", "iou_threshold"), ("--repr_threshold", "repr_threshold"),
                       ("--merge_clusters", "merge_k")):
        kind = {"verbose": dict(action="store_true"), "window_size": dict(type=int), "merge_k": dict(type=int),
                "iou_threshold": dict(type=float), "repr_threshold": dict(type=float)}.get(dest, {})
        _hidden_alias(r, flag, dest, **kind)
    _hidden_alias(r, "--subset", "chr21_only", action="store_true")
    for flag, _ in LEGACY_MODE_FLAGS:
        _hidden_alias(r, "--" + flag, "legacy_" + flag, action="store_true")
    r.add_argument("-q", dest="legacy_quick", action="store_true", help=argparse.SUPPRESS, default=argparse.SUPPRESS)
    return parser


def _resolve_mode(args, error):
    legacy = [mode for flag, mode in LEGACY_MODE_FLAGS if getattr(args, "legacy_" + flag, False)]
    if args.mode and legacy and legacy[0] != args.mode:
        error("--mode %s conflicts with the legacy flag for mode %s" % (args.mode, legacy[0]))
    return args.mode or (legacy[0] if legacy else "full")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "parse":
        os.makedirs(args.outdir, exist_ok=True)
        parse_posteriors(args.posteriordir, args.resolution, args.outdir, args.saga, args.out_format)
    else:
        mode = _resolve_mode(args, parser.error)
        run_sagaconf(args.base, args.verif, args.outdir, mode=mode,
                     base_mnemonics=args.base_mnemonics, verif_mnemonics=args.verif_mnemonics,
                     chr21_only=args.chr21_only, window_size=args.window_size,
                     iou_threshold=args.iou_threshold, repr_threshold=args.repr_threshold,
                     merge_k=args.merge_k, verbose=args.verbose, ccre_file=args.ccre_file,
                     meuleman_file=args.meuleman_file)


# ---- legacy entry points (SAGAconf.py, SAGAconf_parser.py) ------------------------------------

def legacy_run_main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "base", help="path of parsed posterior CSV or BED file of base replicate annotation", type=str)
    parser.add_argument(
        "verif", help="path of parsed posterior CSV or BED file of verification replicate annotation", type=str)
    parser.add_argument(
        "savedir", help="directory to save SAGAconf results.", type=str)

    parser.add_argument(
        "-v", "--verbosity", help="increase output verbosity", action="store_true", default=False)

    parser.add_argument(
        "-bm", "--base_mnemonics", help="if specified, a txt file is used as mnemonics for base replicate.", type=str, default="NA")
    parser.add_argument(
        "-vm", "--verif_mnemonics", help="if specified, a txt file is used as mnemonics for verif replicate.", type=str, default="NA")

    parser.add_argument(
        "-s", "--subset", help="if specified, run SAGA on just one chromosome", action="store_true", default=False)

    parser.add_argument(
        "--active_regions", help="if specified, run SAGAconf on just active regions identified by cCREs or Meuleman et al.", action="store_true", default=False)

    parser.add_argument(
        "--v_seglength", help="if specified, get reproducibility as a function of position relative to segment length", action="store_true", default=False)

    parser.add_argument(
        "-w", "--windowsize",  help="window size (bp) to account for around each genomic bin [default=1000bp]", type=int, default=1000)
    parser.add_argument(
        "-to", "--iou_threshold",  help="Threshold on the IoU of overlap for considering a pair of labels as corresponding.", type=float, default=0.75)
    parser.add_argument(
        "-tr", "--repr_threshold",  help="Threshold on the reproducibility score for considering a segment as reproduced.", type=float, default=0.8)

    parser.add_argument(
        "-k", "--merge_clusters",  help="specify k value to merge base annotation states until k states. ", type=int, default=-1)

    parser.add_argument(
        "-q", "--quick",  help="if True, only a subset of essential analysis are performed for quick report.", action="store_true", default=False)

    parser.add_argument(
        "--r_only",  help="if True, only r_values are computed.", action="store_true", default=False)

    parser.add_argument(
        "--ct_only",  help="if True, only general cellype analysis are performed.", action="store_true", default=False)

    parser.add_argument(
        "--merge_only",  help="if specified, only the specified k value is used to merge base annotation states until k states. ", action="store_true", default=False)

    args = parser.parse_args()
    mode = next((mode for flag, mode in LEGACY_MODE_FLAGS if getattr(args, flag)), "full")
    run_sagaconf(args.base, args.verif, args.savedir, mode=mode,
                 base_mnemonics=args.base_mnemonics, verif_mnemonics=args.verif_mnemonics,
                 chr21_only=args.subset, window_size=args.windowsize,
                 iou_threshold=args.iou_threshold, repr_threshold=args.repr_threshold,
                 merge_k=args.merge_clusters, verbose=args.verbosity)


def legacy_parse_main():
    parser = argparse.ArgumentParser()
    parser.add_argument("posteriordir", help="directory with all posterior files inside.", type=str)
    parser.add_argument("resolution", help="resolution of the SAGA model (bp).", type=int)
    parser.add_argument("savedir", help="the directory to save the parsed posterior file.", type=str)
    parser.add_argument("--out_format", help="the format for saving parsed posteriors.", choices=["bed", "csv"], type=str, default="bed")

    parser.add_argument('--saga', required=True, choices=['segway', 'chmm'])
    args = parser.parse_args()
    parse_posteriors(args.posteriordir, args.resolution, args.savedir, args.saga, args.out_format)


if __name__ == "__main__":
    main()
