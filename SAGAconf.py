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
    "-v", "--verbosity", help="increase output verbosity", action="store_true", default=False)

parser.add_argument(
    "-bm", "--base_mnemonics", help="if specified, a txt file is used as mnemonics for base replicate.", type=str, default="NA")
parser.add_argument(
    "-vm", "--verif_mnemonics", help="if specified, a txt file is used as mnemonics for verif replicate.", type=str, default="NA")

# parser.add_argument(
#     "-g", "--genecode", help="path to parsed genecode data", 
#     type=str, default="biovalidation/parsed_genecode_data_hg38_release42.csv")

# parser.add_argument(
#     "--rnaseq", help="path to rnaseq TPM per gene data", type=str, default=None)

# parser.add_argument(
#     "-c", "--contour", help="if specified, generate contours plots.", action="store_true", default=False)

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

##############################################################################################################################
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
##############################################################################################################################

if args.quick:
    # try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=True, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        quick_report(loci1, loci2, args.savedir, locis=True,  w=w, to=args.iou_threshold, tr=args.repr_threshold)

    # except:
    #     if args.verbosity:
    #         print("Failed to generated quick report.")

elif args.r_only:
    # try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=True, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        # get_overalls(loci1, loci2, args.savedir, locis=True, w=w, to=args.iou_threshold, tr=args.repr_threshold)
        rvalues = is_repr_posterior(
            loci1, loci2, ovr_threshold=args.iou_threshold, window_bp=w, matching="static",
            always_include_best_match=True, return_r=True)

        rvalues.to_csv(args.savedir+f"/r_values.bed", sep='\t', header=True, index=False)

    # except:
    #     if args.verbosity:
    #         print("failed to get GW SAGAconf reproducibility results")

elif args.active_regions:
    # try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        loci1, loci2 = subset_data_to_activeregions(
            replicate_1_dir=loci1, replicate_2_dir=loci2,
            cCREs_file="src/biointerpret/GRCh38-cCREs.bed",
            Meuleman_file="src/biointerpret/Meuleman.tsv", restrict_to="WG", locis=True)

        get_rvals_activeregion(loci1, loci2, args.savedir, w=w, restrict_to="WG")

        # subset the loci_1 and loci_2 to only active regions
        ccre_loci1, ccre_loci2 = subset_data_to_activeregions(
            replicate_1_dir=loci1, replicate_2_dir=loci2,
            cCREs_file="src/biointerpret/GRCh38-cCREs.bed",
            Meuleman_file="src/biointerpret/Meuleman.tsv", restrict_to="cCRE", locis=True)

        get_rvals_activeregion(ccre_loci1, ccre_loci2, args.savedir, w=w, restrict_to="cCRE")

        muel_loci1, muel_loci2 = subset_data_to_activeregions(
            replicate_1_dir=loci1, replicate_2_dir=loci2,
            cCREs_file="src/biointerpret/GRCh38-cCREs.bed",
            Meuleman_file="src/biointerpret/Meuleman.tsv", restrict_to="muel", locis=True)

        get_rvals_activeregion(muel_loci1, muel_loci2, args.savedir, w=w, restrict_to="muel")

        # heatmaps_on_active_regions(
        #     replicate_1_dir=loci1, 
        #     replicate_2_dir=loci2, 
        #     savedir=args.savedir,
        #     cCREs_file="src/biointerpret/GRCh38-cCREs.bed",
        #     Meuleman_file="src/biointerpret/Meuleman.tsv", 
        #     locis=True, w=w)

    # except:
    #     if args.verbosity:
    #         print("failed to get GW SAGAconf reproducibility results on active regions!")

elif args.v_seglength:
    # try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=True, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        overlap_vs_segment_length(
            replicate_1_dir=loci1, 
            replicate_2_dir=loci2, 
            savedir=args.savedir, 
            locis=True)

    # except:
    #     if args.verbosity:
    #         print("failed to get reproducibility results relative to segment length!")

elif args.ct_only: 
    # try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=True, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        get_all_ct(loci1, loci2, args.savedir, locis=True, w=w)

    # except:
    #     if args.verbosity:
    #         print("Failed to generated sample analysis.")

elif args.merge_only:
    # try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        post_clustering_keep_k_states(loci1, loci2, args.savedir, k=args.merge_clusters, locis=True, write_csv=False, w=w)

    # except:
    #     if args.verbosity:
    #         print("Failed to merge clusters up to k")

else:
    # try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=True, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        get_all_ct(loci1, loci2, args.savedir, locis=True, w=w)

    # except:
    #     if args.verbosity:
    #         print("Failed to generated sample analysis.")


    # try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        post_clustering(loci1, loci2, args.savedir, locis=True, to=args.iou_threshold, tr=args.repr_threshold)

    # except:
    #     if args.verbosity:
    #         print("Failed to perform post-clustering")

    if args.merge_clusters !=-1:
        # try:
            loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
            loci1, loci2 = process_data(
                loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
                vm=args.verif_mnemonics, match=False, custom_order=True)

            post_clustering_keep_k_states(loci1, loci2, args.savedir, k=args.merge_clusters, locis=True, write_csv=False)

        # except:
        #     if args.verbosity:
        #         print("Failed to merge clusters up to k")
        
    # try:
        loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=True, force_WG=False)
        loci1, loci2 = process_data(
            loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
            vm=args.verif_mnemonics, match=False, custom_order=True)

        get_overalls(loci1, loci2, args.savedir, locis=True, w=w, to=args.iou_threshold, tr=args.repr_threshold)

    # except:
    #     if args.verbosity:
    #         print("failed to get GW SAGAconf reproducibility results")

##############################################################################################################################

os.system(f"rm -rf {replicate_1_dir}")
os.system(f"rm -rf {replicate_2_dir}")

# listofres = os.listdir(args.savedir)
# main_res = ["confident_segments", "states_post_clustered_posterior"]
# os.mkdir(args.savedir + "/analysis")
# for r in listofres:
#     if main_res[0] not in r and main_res[1] not in r:
#         os.system(f"mv {args.savedir}/{r} {args.savedir}/analysis/")

##############################################################################################################################

# try:
#     loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
#     loci1, loci2 = process_data(
#         loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
#         vm=args.verif_mnemonics, match=False, custom_order=True)

#     get_all_labels(loci1, loci2, args.savedir, locis=True)
# except:
#     if args.verbosity:
#         print("Failed to generated some of the label-specific analysis.")
# 
#    
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


# if args.contour:
#     try:
#         loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=False, force_WG=False)
#         loci1, loci2 = process_data(
#             loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
#             vm=args.verif_mnemonics, match=False, custom_order=True)

#         get_contour(loci1, loci2, args.savedir, locis=True)
#     except:
#         if args.verbosity:
#             print("failed to generate overall reproducibility contours")
# try:
#     gather_labels(args.base, args.savedir, contour=args.contour)

# except:
#     if args.verbosity:
#         print("failed to gather per-label results")

# try:
#     loci1, loci2 = load_data(posterior1_dir, posterior2_dir, subset=issubset, logit_transform=True, force_WG=False)
#     loci1, loci2 = process_data(
#         loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=mnem, bm=args.base_mnemonics, 
#         vm=args.verif_mnemonics, match=False, custom_order=True)

#     after_SAGAconf_metrics(loci1, loci2, args.genecode, args.savedir, rnaseq=None, locis=True, w=w, to=args.iou_threshold, tr=args.repr_threshold)
#     before_after_saga(args.savedir)

# except:
#     if args.verbosity:
#         print("failed to generate before vs. after SAGAconf results")