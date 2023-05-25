from reports import *
import argparse

"""
# TODO:
1. the input files should be in BED3+k format
2. the main output file should be in BED format of confident segments
3. the savedir should look sth like this:
    a. the BED of confident segments
    b. analysis/ directory where everything else reside
4. remove the replicate dirs that are temporarily created
5. add the option to post_cluster the annotations up to a certain point to get better reproducibility.
"""

parser = argparse.ArgumentParser()
parser.add_argument(
    "base", help="path of parsed posterior CSV or BED file of base replicate annotation", type=str)
parser.add_argument(
    "verif", help="path of parsed posterior CSV or BED file of verification replicate annotation", type=str)
parser.add_argument(
    "savedir", help="directory to save SAGAconf results.", type=str)

parser.add_argument(
    "-v", "--verbosity", help="increase output verbosity", action="store_true")

parser.add_argument(
    "-bm", "--base_mnemonics", help="if specified, a txt file is used as mnemonics for base replicate.", type=str, default="NA")
parser.add_argument(
    "-vm", "--verif_mnemonics", help="if specified, a txt file is used as mnemonics for verif replicate.", type=str, default="NA")

# parser.add_argument(
#     "-g", "--genecode", help="path to parsed genecode data", 
#     type=str, default="biovalidation/parsed_genecode_data_hg38_release42.csv")

# parser.add_argument(
#     "--rnaseq", help="path to rnaseq TPM per gene data", type=str, default=None)

parser.add_argument(
    "-c", "--contour", help="if specified, generate contours plots.", action="store_true", default=False)
parser.add_argument(
    "-s", "--subset", help="if specified, run SAGA on just one chromosome", action="store_true", default=True)

parser.add_argument(
    "-w", "--windowsize",  help="window size (bp) to account for around each genomic bin [default=1000bp]", type=int, default=1000)
parser.add_argument(
    "-to", "--iou_threshold",  help="Threshold on the IoU of overlap for considering a pair of labels as corresponding.", type=float, default=0.75)
parser.add_argument(
    "-tr", "--repr_threshold",  help="Threshold on the reproducibility score for considering a segment as reproduced.", type=float, default=0.8)

parser.add_argument(
    "-k", "--merge_clusters",  help="specify k value to merge base annotation states until k states. ", type=int, default=-1)

args = parser.parse_args()

##############################################################################################################################

if os.path.exists(args.savedir)==False:
        os.mkdir(args.savedir)

# we need replicate_1_dir and replicate_2_dir and savedir

if os.path.exists(args.savedir + "/base_replicate")==False:
    os.mkdir(args.savedir + "/base_replicate")

replicate_1_dir = args.savedir + "/base_replicate"

if os.path.exists(args.savedir + "/verification_replicate")==False:
    os.mkdir(args.savedir + "/verification_replicate")

replicate_2_dir = args.savedir + "/verification_replicate"

# copy parsed_posterior files to replicate_*_dir
if ".bed" in args.base.lower():
    os.system(f"cp {args.base} {replicate_1_dir}/parsed_posterior.bed")
    posterior1_dir = f"{replicate_1_dir}/parsed_posterior.bed"
elif ".csv" in args.base.lower():
    os.system(f"cp {args.base} {replicate_1_dir}/parsed_posterior.csv")
    posterior1_dir = f"{replicate_1_dir}/parsed_posterior.csv"

if ".bed" in args.verif.lower():
    os.system(f"cp {args.verif} {replicate_2_dir}/parsed_posterior.bed")
    posterior2_dir = f"{replicate_2_dir}/parsed_posterior.bed"
elif ".csv" in args.verif.lower():
    os.system(f"cp {args.verif} {replicate_2_dir}/parsed_posterior.csv")
    posterior2_dir = f"{replicate_2_dir}/parsed_posterior.csv"

# optionally, copy mnemonics into replicate folders.
if args.base_mnemonics != "NA":
    os.system(f"cp {args.base_mnemonics} {replicate_1_dir}/mnemonics.txt")

if args.verif_mnemonics != "NA":
    os.system(f"cp {args.verif_mnemonics} {replicate_2_dir}/mnemonics.txt")

if args.base_mnemonics != "NA" and args.verif_mnemonics != "NA":
    mnem = True
else:
    mnem = False

w = args.windowsize
issubset = args.subset

##############################################################################################################################

try:
    loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=True, force_WG=False)
    loci1, loci2 = process_data(
        loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
        vm=args.verif_mnemonics, match=False, custom_order=True)

    get_all_ct(loci1, loci2, args.savedir, locis=True, w=w)

except:
    if args.verbosity:
        print("Failed to generated some of the general celltype analysis.")
try:
    loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
    loci1, loci2 = process_data(
        loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
        vm=args.verif_mnemonics, match=False, custom_order=True)

    get_all_labels(loci1, loci2, args.savedir, locis=True)
except:
    if args.verbosity:
        print("Failed to generated some of the label-specific analysis.")

try:
    loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
    loci1, loci2 = process_data(
        loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
        vm=args.verif_mnemonics, match=False, custom_order=True)

    post_clustering(loci1, loci2, args.savedir, locis=True, to=args.iou_threshold, tr=args.repr_threshold)

except:
    if args.verbosity:
        print("Failed to perform post-clustering")

if args.merge_clusters !=-1:
    try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        post_clustering_keep_k_states(loci1, loci2, args.savedir, k=args.merge_clusters, locis=True, write_csv=True)

    except:
        if args.verbosity:
            print("Failed to merge clusters")
    
        
# try:
#     loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
#     loci1, loci2 = process_data(
#         loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
#         vm=args.verif_mnemonics, match=False, custom_order=True)

#     get_all_bioval(
#         loci1, loci2, 
#         args.savedir,
#         genecode_dir=args.genecode, 
#         rnaseq=args.rnaseq, locis=True)
# except:
#     if args.verbosity:
#         print("failed to perform biological validation results")

try:
    loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=True, force_WG=False)
    loci1, loci2 = process_data(
        loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
        vm=args.verif_mnemonics, match=False, custom_order=True)

    get_overalls(loci1, loci2, args.savedir, locis=True, w=w, to=args.iou_threshold, tr=args.repr_threshold)

except:
    if args.verbosity:
        print("failed to get overall SAGAconf reproducibility results")

if args.contour:
    try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        get_contour(loci1, loci2, args.savedir, locis=True)
    except:
        if args.verbosity:
            print("failed to generate overall reproducibility contours")

try:
    gather_labels(args.base, args.savedir, contour=args.contour)

except:
    if args.verbosity:
        print("failed to gather per-label results")

try:
    loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=True, force_WG=False)
    loci1, loci2 = process_data(
        loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
        vm=args.verif_mnemonics, match=False, custom_order=True)

    after_SAGAconf_metrics(loci1, loci2, args.genecode, args.savedir, rnaseq=None, locis=True, w=w, to=args.iou_threshold, tr=args.repr_threshold)
    before_after_saga(args.savedir)

except:
    if args.verbosity:
        print("failed to generate before vs. after SAGAconf results")


os.system(f"rm -rf {replicate_1_dir}")
os.system(f"rm -rf {replicate_2_dir}")

listofres = os.listdir(args.savedir)
main_res = ["confident_segments", "states_post_clustered_posterior"]
os.mkdir(args.savedir + "/analysis")
for r in listofres:
    if main_res[0] not in r and main_res[1] not in r:
        os.system(f"mv {args.savedir}/{r} {args.savedir}/analysis/")