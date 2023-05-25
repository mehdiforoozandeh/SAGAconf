# import argparse
import os, argparse
import multiprocessing as mp
from _utils import *
from _chromhmm import *

parser = argparse.ArgumentParser()
parser.add_argument("posteriordir", help="directory with all posterior files inside.", type=str)
parser.add_argument("resolution", help="resolution of the SAGA model (bp).", type=int)
parser.add_argument("savedir", help="the directory to save the parsed posterior file.", type=str)
parser.add_argument("--out_format", help="the format for saving parsed posteriors.", choices=["bed", "csv"], type=str, default="bed")

parser.add_argument('--saga', required=True, choices=['segway', 'chmm'])
args = parser.parse_args()

if args.saga == "segway":
    binned_posterior = mp_inplace_binning(args.posteriordir, args.resolution)
elif args.saga == "chmm":
    binned_posterior = ChrHMM_read_posteriordir(args.posteriordir, args.resolution)

if args.out_format == "bed":
    # Write the DataFrame to a BED file
    binned_posterior.to_csv(args.savedir + "/parsed_posterior.bed", sep='\t', header=True, index=False)
elif args.out_format == "csv":
    binned_posterior.to_csv(args.savedir + "/parsed_posterior.csv")
