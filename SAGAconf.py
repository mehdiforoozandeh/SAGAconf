from reports import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("base", help="path of parsed posterior CSV file of base replicate annotation.", type=str)
parser.add_argument("verif", help="path of parsed posterior CSV file of verification replicate annotation.", type=str)
parser.add_argument("savedir", help="directory to save SAGAconf results.", type=str)

parser.add_argument("-v", "--verbosity", help="increase output verbosity", action="store_true")

parser.add_argument("-bm", "--base_mnemonics", help="increase output verbosity")
parser.add_argument("-vm", "--verif_mnemonics", help="increase output verbosity")

parser.add_argument("--genecode", help="path to parsed genecode data", type=str, default="biovalidation/parsed_genecode_data_hg38_release42.csv")
parser.add_argument("--rnaseq", help="path to rnaseq TPM per gene data", type=str)

parser.add_argument("--contour", help="generate contours plots.", action="store_true")
parser.add_argument("--subset", help="run SAGA on just one chromosome", action="store_true")

parser.add_argument("-w", "--windowsize",  help="window size (bp) to account for around each genomic bin", type=int)

args = parser.parse_args()

# TODO: figure out the how to use load_data and process_data without while keeping the mnemonics and subset options.


if os.path.exists(args.savedir)==False:
        os.mkdir(args.savedir)

try:
    get_all_ct(args.base, args.verif, args.savedir)

except:
    pass

try:
    get_all_labels(args.base, args.verif, args.savedir)
except:
    pass

try:
    post_clustering(args.base, args.verif, args.savedir)

except:
    pass

try:
    get_all_bioval(
        args.base, args.verif, 
        args.savedir,
        genecode_dir=args.genecode, 
        rnaseq=args.rnaseq)
except:
    pass

try:
    get_overalls(args.base, args.verif, args.savedir)

except:
    pass

if args.contour:
    try:
        get_contour(args.base, args.verif, args.savedir)
    except:
        pass

try:
    gather_labels(args.base, args.savedir, contour=args.contour)

except:
    pass

try:
    after_SAGAconf_metrics(args.base, args.verif, args.genecode, args.savedir, rnaseq=None)
    before_after_saga(args.savedir)

except:
    pass