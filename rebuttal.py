import os, pybedtools, random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.interpolate import UnivariateSpline
from sklearn.linear_model import LinearRegression
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.stats import gaussian_kde
from matplotlib.lines import Line2D
from src.bio_valid import *
from scipy.stats import pearsonr
from reports import *
from scipy.stats import ttest_ind, mannwhitneyu, ks_2samp, wilcoxon, zscore
from numpy import mean, var
from math import sqrt
from statsmodels.stats.weightstats import ztest
from matplotlib.ticker import FuncFormatter

def get_listofruns(maindir="rebuttal"):
    listofruns = [
        {"replicate_1_dir":"chromhmm_runs/GM12878_rep1/", 
        "replicate_2_dir":"chromhmm_runs/GM12878_rep2/", 
        "savedir":"{}/r1vsr2/chmm/GM12878/".format(maindir)},
        
        {"replicate_1_dir":"segway_runs/GM12878_rep1/", 
        "replicate_2_dir":"segway_runs/GM12878_rep2/", 
        "savedir":"{}/r1vsr2/segway/GM12878/".format(maindir)},
        ######## MCF-7 ########
    
        {"replicate_1_dir":"chromhmm_runs/MCF-7_rep1/", 
        "replicate_2_dir":"chromhmm_runs/MCF-7_rep2/", 
        "savedir":"{}/r1vsr2/chmm/MCF-7/".format(maindir)},

        {"replicate_1_dir":"segway_runs/MCF-7_rep1/", 
        "replicate_2_dir":"segway_runs/MCF-7_rep2/", 
        "savedir":"{}/r1vsr2/segway/MCF-7/".format(maindir)},

        ######## CD14 ########
    
        {"replicate_1_dir":"chromhmm_runs/CD14-positive_monocyte_rep1/", 
        "replicate_2_dir":"chromhmm_runs/CD14-positive_monocyte_rep2/", 
        "savedir":"{}/r1vsr2/chmm/CD14/".format(maindir)},
    
        {"replicate_1_dir":"segway_runs/CD14-positive_monocyte_rep1/", 
        "replicate_2_dir":"segway_runs/CD14-positive_monocyte_rep2/", 
        "savedir":"{}/r1vsr2/segway/CD14/".format(maindir)},
        
        ######## K562 ########
    
        {"replicate_1_dir":"chromhmm_runs/K562_rep1/", 
        "replicate_2_dir":"chromhmm_runs/K562_rep2/", 
        "savedir":"{}/r1vsr2/chmm/K562/".format(maindir)},
    
        {"replicate_1_dir":"segway_runs/K562_rep1/", 
        "replicate_2_dir":"segway_runs/K562_rep2/", 
        "savedir":"{}/r1vsr2/segway/K562/".format(maindir)},

        ######## HeLa-S3 ########
    
        {"replicate_1_dir":"chromhmm_runs/HeLa-S3_rep1/", 
        "replicate_2_dir":"chromhmm_runs/HeLa-S3_rep2/", 
        "savedir":"{}/r1vsr2/chmm/HeLa-S3/".format(maindir)},

        {"replicate_1_dir":"segway_runs/HeLa-S3_rep1/", 
        "replicate_2_dir":"segway_runs/HeLa-S3_rep2/", 
        "savedir":"{}/r1vsr2/segway/HeLa-S3/".format(maindir)},

        ######## GM12878 ########
        {"replicate_1_dir":"chromhmm_runs/GM12878_concat/", 
        "replicate_2_dir":"chromhmm_runs/GM12878_concat/", 
        "savedir":"{}/concat/chmm/GM12878/".format(maindir)},

        {"replicate_1_dir":"segway_runs/GM12878_concat_rep1/", 
        "replicate_2_dir":"segway_runs/GM12878_concat_rep2/", 
        "savedir":"{}/concat/segway/GM12878/".format(maindir)},
            
        ######## MCF-7 ########

        {"replicate_1_dir":"chromhmm_runs/MCF-7_concat/", 
        "replicate_2_dir":"chromhmm_runs/MCF-7_concat/", 
        "savedir":"{}/concat/chmm/MCF-7/".format(maindir)},

        {"replicate_1_dir":"segway_runs/MCF-7_concat_rep1/", 
        "replicate_2_dir":"segway_runs/MCF-7_concat_rep2/", 
        "savedir":"{}/concat/segway/MCF-7/".format(maindir)},

        ######## CD14 ########

        {"replicate_1_dir":"chromhmm_runs/CD14-positive_monocyte_concat/", 
        "replicate_2_dir":"chromhmm_runs/CD14-positive_monocyte_concat/", 
        "savedir":"{}/concat/chmm/CD14/".format(maindir)},

        {"replicate_1_dir":"segway_runs/CD14-positive_monocyte_concat_rep1/", 
        "replicate_2_dir":"segway_runs/CD14-positive_monocyte_concat_rep2/", 
        "savedir":"{}/concat/segway/CD14/".format(maindir)},
        
        ######## K562 ########

        {"replicate_1_dir":"chromhmm_runs/K562_concat/", 
        "replicate_2_dir":"chromhmm_runs/K562_concat/", 
        "savedir":"{}/concat/chmm/K562/".format(maindir)},

        {"replicate_1_dir":"segway_runs/K562_concat_rep1/", 
        "replicate_2_dir":"segway_runs/K562_concat_rep2/", 
        "savedir":"{}/concat/segway/K562/".format(maindir)},

        ######## HeLa-S3 ########

        {"replicate_1_dir":"chromhmm_runs/HeLa-S3_concat/", 
        "replicate_2_dir":"chromhmm_runs/HeLa-S3_concat/", 
        "savedir":"{}/concat/chmm/HeLa-S3/".format(maindir)},

        {"replicate_1_dir":"segway_runs/HeLa-S3_concat_rep1/", 
        "replicate_2_dir":"segway_runs/HeLa-S3_concat_rep2/", 
        "savedir":"{}/concat/segway/HeLa-S3/".format(maindir)},
        ######## GM12878 ########

        {"replicate_1_dir":"chromhmm_runs/GM12878_rep2_rs5/", 
        "replicate_2_dir":"chromhmm_runs/GM12878_rep2_rs27/", 
        "savedir":"{}/paraminit/chmm/GM12878/".format(maindir)},


        {"replicate_1_dir":"segway_runs/GM12878_rep1/", 
        "replicate_2_dir":"segway_runs/GM12878_rep1_rs7/", 
        "savedir":"{}/paraminit/segway/GM12878/".format(maindir)},
            
        ######## MCF-7 ########

        {"replicate_1_dir":"chromhmm_runs/MCF-7_rep1_rs5/", 
        "replicate_2_dir":"chromhmm_runs/MCF-7_rep1_rs27/", 
        "savedir":"{}/paraminit/chmm/MCF-7/".format(maindir)},


        {"replicate_1_dir":"segway_runs/MCF-7_rep1_rs5/", 
        "replicate_2_dir":"segway_runs/MCF-7_rep1_rs7/", 
        "savedir":"{}/paraminit/segway/MCF-7/".format(maindir)},

        ######## CD14 ########

        {"replicate_1_dir":"chromhmm_runs/CD14-positive_monocyte_rep1_rs5/", 
        "replicate_2_dir":"chromhmm_runs/CD14-positive_monocyte_rep1_rs27/", 
        "savedir":"{}/paraminit/chmm/CD14/".format(maindir)},


        {"replicate_1_dir":"segway_runs/CD14-positive_monocyte_rep1_rs5/", 
        "replicate_2_dir":"segway_runs/CD14-positive_monocyte_rep1_rs7/", 
        "savedir":"{}/paraminit/segway/CD14/".format(maindir)},
        
        ######## K562 ########

        {"replicate_1_dir":"chromhmm_runs/K562_rep2_rs5/", 
        "replicate_2_dir":"chromhmm_runs/K562_rep2_rs27/", 
        "savedir":"{}/paraminit/chmm/K562/".format(maindir)},

        {"replicate_1_dir":"segway_runs/K562_rep1/", 
        "replicate_2_dir":"segway_runs/K562_rep1_rs7/", 
        "savedir":"{}/paraminit/segway/K562/".format(maindir)},

        ######## HeLa-S3 ########

        {"replicate_1_dir":"chromhmm_runs/HeLa-S3_rep2_rs5/", 
        "replicate_2_dir":"chromhmm_runs/HeLa-S3_rep2_rs27/", 
        "savedir":"{}/paraminit/chmm/HeLa-S3/".format(maindir)},

        {"replicate_1_dir":"segway_runs/HeLa-S3_rep1/", 
        "replicate_2_dir":"segway_runs/HeLa-S3_rep1_rs7/", 
        "savedir":"{}/paraminit/segway/HeLa-S3/".format(maindir)}
        ]
    return listofruns

def get_single_run(r): # r is run_dict
    savedir = r["savedir"]
    if "chromhmm_runs" in r["replicate_1_dir"] and "concat" in r["replicate_1_dir"]:
        base_mnemonics = r["replicate_1_dir"] + "/mnemonics_rep1.txt"
        verif_mnemonics = r["replicate_2_dir"] + "/mnemonics_rep2.txt"

        replicate_1_dir = r["replicate_1_dir"] + "/parsed_posterior_rep1.csv"
        replicate_2_dir = r["replicate_2_dir"] + "/parsed_posterior_rep2.csv"
    else:
        base_mnemonics = r["replicate_1_dir"] + "/mnemonics.txt"
        verif_mnemonics = r["replicate_2_dir"] + "/mnemonics.txt"

        replicate_1_dir = r["replicate_1_dir"] + "/parsed_posterior.csv"
        replicate_2_dir = r["replicate_2_dir"] + "/parsed_posterior.csv"

    try:
        print(f"trying to get original r-values for {savedir}")
        os.system(f"python SAGAconf.py --r_only -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
    except:
        print(f"failed to get original r-values for {savedir}")

    if "GM12878" in replicate_1_dir:
        expression_data = "src/biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv"
    elif "MCF-7" in replicate_1_dir:
        expression_data = "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl"
    elif "K562" in replicate_1_dir:
        expression_data = "src/biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv"
    else:
        expression_data = "EMPTY"
    
    if expression_data != "EMPTY":
        print(f"trying to get expression analysis for {savedir}")
        os.system(f"python SAGAconf.py --merge_only -k 14 -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
        os.system(f"python SAGAconf.py --merge_only -k 12 -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
        os.system(f"python SAGAconf.py --merge_only -k 10 -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")

        try:
            r_dist_vs_expression3(f"{savedir}/r_values.bed", expression_data, savedir+"/16_states/", interpret=True)
            r_dist_vs_expression3_genebody(f"{savedir}/r_values.bed", expression_data, savedir+"/16_states/", interpret=True)
            conf_v_nonconf_vs_expression(f"{savedir}/r_values.bed", expression_data, savedir+"/16_states/", interpret=True)
            r_dist_vs_expression_boxplot(f"{savedir}/r_values.bed", expression_data, savedir+"/16_states/", interpret=True)
        except:
            pass

        try:
            r_dist_vs_expression3(f"{savedir}/r_values_14_states.bed", expression_data, savedir+"/14_states/", interpret=True)
            r_dist_vs_expression3_genebody(f"{savedir}/r_values_14_states.bed", expression_data, savedir+"/14_states/", interpret=True)
            conf_v_nonconf_vs_expression(f"{savedir}/r_values_14_states.bed", expression_data, savedir+"/14_states/", interpret=True)
            r_dist_vs_expression_boxplot(f"{savedir}/r_values_14_states.bed", expression_data, savedir+"/14_states/", interpret=True)
        except:
            pass

        
        try:
            r_dist_vs_expression3(f"{savedir}/r_values_12_states.bed", expression_data, savedir+"/12_states/", interpret=True)
            r_dist_vs_expression3_genebody(f"{savedir}/r_values_12_states.bed", expression_data, savedir+"/12_states/", interpret=True)
            conf_v_nonconf_vs_expression(f"{savedir}/r_values_12_states.bed", expression_data, savedir+"/12_states/", interpret=True)
            r_dist_vs_expression_boxplot(f"{savedir}/r_values_12_states.bed", expression_data, savedir+"/12_states/", interpret=True)
        except:
            pass
        
        try:
            r_dist_vs_expression3(f"{savedir}/r_values_10_states.bed", expression_data, savedir+"/10_states/", interpret=True)
            r_dist_vs_expression3_genebody(f"{savedir}/r_values_10_states.bed", expression_data, savedir+"/10_states/", interpret=True)
            conf_v_nonconf_vs_expression(f"{savedir}/r_values_10_states.bed", expression_data, savedir+"/10_states/", interpret=True)
            r_dist_vs_expression_boxplot(f"{savedir}/r_values_10_states.bed", expression_data, savedir+"/10_states/", interpret=True)
        except:
            pass

    try:
        print(f"trying to get per-segment analysis for {savedir}")
        r_distribution_over_segment(f"{savedir}/r_values.bed", savedir)
        os.system(f"python SAGAconf.py --v_seglength -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
    except:
        print(f"failed to get per-segment analysis for {savedir}")

    try:
        print(f"trying to get ccre analysis for {savedir}")
        os.system(f"python SAGAconf.py --active_regions -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
        r_distribution_activeregions2(
        f"{savedir}/r_values_WG.bed", 
        f"{savedir}/r_values_cCRE.bed", 
        f"{savedir}/r_values_muel.bed", 
        savedir)
    except:
        print(f"failed to get ccre analysis for {savedir}")
    
    """
    get original r_values 
    get original r_distribution_over_segment 
    get original naive overlap v_seglength
    
    get SAGAconf on WG and active_regions 
    get original r_distribution_activeregions2

    get original r_vs_expression + r_vs_expression_genebody
    get post_clustered r_value 14, 12, 10
    get r_vs_expression + r_vs_expression_genebody for 14, 12, 10 states
    """

def get_subset_transc(r):
    savedir = r["savedir"]
    if "chromhmm_runs" in r["replicate_1_dir"] and "concat" in r["replicate_1_dir"]:
        base_mnemonics = r["replicate_1_dir"] + "/mnemonics_rep1.txt"
        verif_mnemonics = r["replicate_2_dir"] + "/mnemonics_rep2.txt"

        replicate_1_dir = r["replicate_1_dir"] + "/parsed_posterior_rep1.csv"
        replicate_2_dir = r["replicate_2_dir"] + "/parsed_posterior_rep2.csv"
    else:
        base_mnemonics = r["replicate_1_dir"] + "/mnemonics.txt"
        verif_mnemonics = r["replicate_2_dir"] + "/mnemonics.txt"

        replicate_1_dir = r["replicate_1_dir"] + "/parsed_posterior.csv"
        replicate_2_dir = r["replicate_2_dir"] + "/parsed_posterior.csv"

    if "GM12878" in replicate_1_dir:
        expression_data = "src/biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv"
    elif "MCF-7" in replicate_1_dir:
        expression_data = "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl"
    elif "K562" in replicate_1_dir:
        expression_data = "src/biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv"
    else:
        expression_data = "EMPTY"
    
    if expression_data != "EMPTY":
        print(f"trying to get expression analysis for {savedir}")
        # try:
        if os.path.exists(f"{savedir}/r_values.bed") == False:
            os.system(f"python SAGAconf.py -s --r_only -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
        # r_dist_vs_expression3(f"{savedir}/r_values.bed", expression_data, savedir+"/16_states/", interpret=True)
        # r_dist_vs_expression3_genebody(f"{savedir}/r_values.bed", expression_data, savedir+"/16_states/", interpret=True)
        r_dist_vs_expression_boxplot(f"{savedir}/r_values.bed", expression_data, savedir+"/16_states/", interpret=True)
        conf_v_nonconf_vs_expression(f"{savedir}/r_values.bed", expression_data, savedir+"/16_states/", interpret=True)
        return
        # except Exception as e:
        #     print("ERROR:   ", e)

        try:
            os.system(f"python SAGAconf.py --merge_only -s -k 14 -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
            # r_dist_vs_expression3(f"{savedir}/r_values_14_states.bed", expression_data, savedir+"/14_states/", interpret=True)
            # r_dist_vs_expression3_genebody(f"{savedir}/r_values_14_states.bed", expression_data, savedir+"/14_states/", interpret=True)
            conf_v_nonconf_vs_expression(f"{savedir}/r_values_14_states.bed", expression_data, savedir+"/14_states/", interpret=True)
            r_dist_vs_expression_boxplot(f"{savedir}/r_values_14_states.bed", expression_data, savedir+"/14_states/", interpret=True)
        except Exception as e:
            print("ERROR:   ", e)
     
        try:
            os.system(f"python SAGAconf.py --merge_only -s -k 12 -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
            # r_dist_vs_expression3(f"{savedir}/r_values_12_states.bed", expression_data, savedir+"/12_states/", interpret=True)
            # r_dist_vs_expression3_genebody(f"{savedir}/r_values_12_states.bed", expression_data, savedir+"/12_states/", interpret=True)
            conf_v_nonconf_vs_expression(f"{savedir}/r_values_12_states.bed", expression_data, savedir+"/12_states/", interpret=True)
            r_dist_vs_expression_boxplot(f"{savedir}/r_values_12_states.bed", expression_data, savedir+"/12_states/", interpret=True)
        except Exception as e:
            print("ERROR:   ", e)
        
        try:
            os.system(f"python SAGAconf.py --merge_only -s -k 10 -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
            # r_dist_vs_expression3(f"{savedir}/r_values_10_states.bed", expression_data, savedir+"/10_states/", interpret=True)
            # r_dist_vs_expression3_genebody(f"{savedir}/r_values_10_states.bed", expression_data, savedir+"/10_states/", interpret=True)
            conf_v_nonconf_vs_expression(f"{savedir}/r_values_10_states.bed", expression_data, savedir+"/10_states/", interpret=True)
            r_dist_vs_expression_boxplot(f"{savedir}/r_values_10_states.bed", expression_data, savedir+"/10_states/", interpret=True)
        except Exception as e:
            print("ERROR:   ", e)

def get_runs(maindir = "rebuttal", mp=True, n_processes=10, subset=False):
    list_of_runs = get_listofruns(maindir)
    random.shuffle(list_of_runs)

    if os.path.exists(maindir)==False:
        os.mkdir(maindir)

    if os.path.exists(maindir+"/r1vsr2")==False:
        os.mkdir(maindir+"/r1vsr2")

    if os.path.exists(maindir+"/r1vsr2/chmm")==False:
        os.mkdir(maindir+"/r1vsr2/chmm")

    if os.path.exists(maindir+"/r1vsr2/segway")==False:
        os.mkdir(maindir+"/r1vsr2/segway")
    
    if os.path.exists(maindir+"/concat")==False:
        os.mkdir(maindir+"/concat")
    
    if os.path.exists(maindir+"/concat/chmm")==False:
        os.mkdir(maindir+"/concat/chmm")

    if os.path.exists(maindir+"/concat/segway")==False:
        os.mkdir(maindir+"/concat/segway")

    if os.path.exists(maindir+"/paraminit")==False:
        os.mkdir(maindir+"/paraminit")
    
    if os.path.exists(maindir+"/paraminit/chmm")==False:
        os.mkdir(maindir+"/paraminit/chmm")

    if os.path.exists(maindir+"/paraminit/segway")==False:
        os.mkdir(maindir+"/paraminit/segway")

    if mp:
        with Pool(n_processes) as p:
            if subset:
                p.map(get_subset_transc, list_of_runs)
            else:
                p.map(get_single_run, list_of_runs)

            p.map(corresp_emiss_v_iou, list_of_runs)
    else:
        for r in list_of_runs:
            if subset:
                get_subset_transc(r)
            else:
                get_single_run(r)
            corresp_emiss_v_iou(r)

def corresp_emiss_v_iou(r):
    savedir = r["savedir"]
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)
    replicate_1_dir = r["replicate_1_dir"]
    replicate_2_dir = r["replicate_2_dir"]
    try:
        if "chromhmm_runs" in replicate_1_dir:
            compare_corresp_methods(replicate_1_dir, replicate_2_dir, savedir, saga="chmm")
        else:
            compare_corresp_methods(replicate_1_dir, replicate_2_dir, savedir, saga="segway")
    except:
        pass

def r_distribution_over_segment(r_value_file, savedir, custom_bin=True):
    interpretation_terms = ["Prom", "Prom_fla", "Enha", "Enha_low", "Biva", "Tran", "Cons", "Facu", "K9K3", "Quie", "Unkn"]
    r_vals = pd.read_csv(r_value_file, sep="\t").to_numpy()
    all_segs = {term: [] for term in interpretation_terms}

    resolution = r_vals[0, 2] - r_vals[0, 1]
    current_map = r_vals[0, 3]
    current_start = r_vals[0, 1]
    current_seg = []

    for i in range(1, len(r_vals)):
        if r_vals[i, 0] == r_vals[i-1, 0] and r_vals[i, 3] == current_map:
            #middle of the segment
            current_seg.append(
                [r_vals[i, 1] - current_start, 
                r_vals[i, 4]])
            
        else:
            #last_segment_ends
            seg_length = r_vals[i-1, 2] - current_start
            current_seg = np.reshape(np.array(current_seg), (-1, 2))
            if len(current_seg) > 1:
                current_seg[:, 0] = [j/(len(current_seg)-1) for j in range(len(current_seg))]
            else:
                current_seg[:, 0] = np.array([0])
            translate_to_term = max([x for x in interpretation_terms if x in current_map], key=len)
            all_segs[translate_to_term].append(current_seg)

            #new_segment
            current_map = r_vals[i, 3]
            current_start = r_vals[i, 1]
            current_seg = []

            current_seg.append(
                [r_vals[i, 1] - current_start, 
                r_vals[i, 4]])

    def custom_binning(x, y, n_bins):
        """
        This function bins the input array x into n_bins and averages the corresponding values in y for each bin.
        
        Parameters:
        x (array-like): Input array to be binned.
        y (array-like): Values to be averaged over each bin of x.
        n_bins (int): Number of bins to divide x into.

        Returns:
        binned_x (array-like): The center value of each bin.
        avg_y (array-like): The average y value for each bin.
        """
        
        # Define the range of x
        x_range = np.linspace(np.min(x), np.max(x), n_bins+1)
        
        # Bin x using numpy's digitize function
        indices = np.digitize(x, x_range)
        
        # Initialize an empty list to hold the average y values for each bin
        avg_y = []
        
        # Initialize an empty list to hold the center value of each bin
        binned_x = []
        
        # For each bin, calculate the average y value and the center value of the bin
        for i in range(1, n_bins+1):
            mean_y = np.mean(y[indices == i])

            if np.isnan(mean_y):
                avg_y.append(avg_y[-1])
            else:
                avg_y.append(mean_y)

            binned_x.append((x_range[i] + x_range[i-1]) / 2)
        
        return binned_x, avg_y

    if not custom_bin:
        splines_100_1k = {}
        splines_1k_10k = {}
        splines_10k_plus = {}

    else:
        n_bins = 5
        min_samples = 100
        binned_100_1k = {}
        binned_1k_10k = {}
        binned_10k_plus = {}

    for k in all_segs.keys():
        try:
            if len([seg for seg in all_segs[k] if 100 < (len(seg)*resolution) <= 1000]) > 0:
                subset1 = np.concatenate([seg for seg in all_segs[k] if 100 < (len(seg)*resolution) <= 1000])
                sorted_indices = np.argsort(subset1[:, 0])
                subset1 = subset1[sorted_indices]
                x1 = subset1[:, 0]
                y1 = subset1[:, 1]
                if len(x1) >= min_samples:
                    if custom_bin:
                        x1, y1 = custom_binning(x1, y1, n_bins)
                        # print(len(x1), len(y1))
                        binned_100_1k[k] = [x1, y1]
                    else:
                        splines_100_1k[k] = UnivariateSpline(x1, y1, k=5)

        except:
            pass

        try:
            if len([seg for seg in all_segs[k] if 1000 < (len(seg)*resolution) <= 10000]) > 0:
                subset2 = np.concatenate([seg for seg in all_segs[k] if 1000 < (len(seg)*resolution) <= 10000])
                sorted_indices = np.argsort(subset2[:, 0])
                subset2 = subset2[sorted_indices]
                x2 = subset2[:, 0]
                y2 = subset2[:, 1]
                if len(x2) >= min_samples:
                    if custom_bin:
                        x2, y2 = custom_binning(x2, y2, n_bins)
                        # print(len(x2), len(y2))
                        binned_1k_10k[k] = [x2, y2]
                    else:
                        splines_1k_10k[k] = UnivariateSpline(x2, y2, k=5)

        except:
            pass

        try:
            if len([seg for seg in all_segs[k] if 10000 < (len(seg)*resolution)]):
                subset3 = np.concatenate([seg for seg in all_segs[k] if 10000 < (len(seg)*resolution)])
                sorted_indices = np.argsort(subset3[:, 0])
                subset3 = subset3[sorted_indices]
                x3 = subset3[:, 0]
                y3 = subset3[:, 1]
                if len(x3) >= min_samples:
                    if custom_bin:
                        x3, y3 = custom_binning(x3, y3, n_bins)
                        # print(len(x3), len(y3))
                        binned_10k_plus[k] = [x3, y3]
                    else:
                        splines_10k_plus[k] = UnivariateSpline(x3, y3, k=5)

        except:
            pass
    
    # Create a new figure with 3 subplots
    fig = plt.figure(figsize=(15, 5))
    gs = gridspec.GridSpec(2, 3, height_ratios=[0.5, 5])

    ax0 = plt.subplot(gs[1, 0])
    ax1 = plt.subplot(gs[1, 1], sharex=ax0, sharey=ax0)
    ax2 = plt.subplot(gs[1, 2], sharex=ax0, sharey=ax0)

    # Create a color map
    colors = plt.cm.get_cmap('rainbow', len(all_segs.keys()))
    lines = []  # list to store the lines for legend
    labels = []  # list to store the labels for legend

    for i, k in enumerate(all_segs.keys()):
        legend_added = False

        # Generate x values
        if not custom_bin:
            x_values = np.linspace(0, 1, 100)

        try:
            # Subplot 1
            if custom_bin:
                x_values, y_values = binned_100_1k[k][0], binned_100_1k[k][1]
            else:
                y_values = splines_100_1k[k](x_values)

            line, = ax0.plot(x_values, y_values, label=k, color=colors(i))

            if not legend_added:
                lines.append(line)
                labels.append(k)
                legend_added = True

            ax0.set_title('Segments with length < 1kb')
            ax0.set_xlabel('Position relative to segment')
            ax0.set_ylabel('r_value')
        except:
            pass

        try:
            # Subplot 2
            if custom_bin:
                x_values, y_values = binned_1k_10k[k][0], binned_1k_10k[k][1]
            else:
                y_values = splines_1k_10k[k](x_values)

            line, = ax1.plot(x_values, y_values, label=k, color=colors(i))
            
            if not legend_added:
                lines.append(line)
                labels.append(k)
                legend_added = True
            
            ax1.set_title('Segments with length 1kb - 10kb')
            ax1.set_xlabel('Position relative to segment')
            ax1.set_ylabel('r_value')
        except:
            pass

        try:
            # Subplot 3
            if custom_bin:
                x_values, y_values = binned_10k_plus[k][0], binned_10k_plus[k][1]
            else:
                y_values = splines_10k_plus[k](x_values)

            line, = ax2.plot(x_values, y_values, label=k, color=colors(i))

            if not legend_added:
                lines.append(line)
                labels.append(k)
                legend_added = True
            
            ax2.set_title('Segments with length > 10kb')
            ax2.set_xlabel('Position relative to segment')
            ax2.set_ylabel('r_value')
        except:
            pass

    # Show the plot
    # Create a separate subplot for the legend at the top
    ax_legend = plt.subplot(gs[0, :])
    ax_legend.axis('off')  # Hide the axes

    # Show the legend in this subplot
    fig.legend(lines, labels, loc='center', ncol=len(labels), bbox_to_anchor=(0.5, 0.5), bbox_transform=ax_legend.transAxes)
    plt.tight_layout()
    plt.savefig(f"{savedir}/v_segment_length.pdf", format='pdf')
    plt.savefig(f"{savedir}/v_segment_length.svg", format='svg')

def get_all_r_distribution_over_segment(maindir="rebuttal"):
    list_of_runs = get_listofruns(maindir)

    for r in list_of_runs:
        r_val_file = r["savedir"] + "/r_value.bed"
        r_distribution_over_segment(r_val_file)

def get_overlap_cCRE(file1, file2):
    # Load the bed files
    bed1 = pybedtools.BedTool(file1)
    bed2 = pd.read_csv(file2, sep="\t")
    bed2.columns = ["chrom",'start', 'end', 'MAP', 'r_value']

    # Convert DataFrame to BedTool
    bed2 = pybedtools.BedTool.from_dataframe(bed2)

    # Get the intersection
    intersected = bed2.intersect(bed1, wa=True, wb=True)

    # Filter the columns
    df = intersected.to_dataframe()
    df = df[['chrom', 'start', 'end', 'name', 'score']]

    # Rename the columns
    df.columns = ['chr', 'start', 'end', 'MAP', 'r_value']

    return df

def get_overlap_Meuleman(file1, file2):
    # Load the bed files
    bed1 = pd.read_csv(file1, sep="\t")
    bed2 = pd.read_csv(file2, sep="\t")
    bed1.columns = ["chr", "start", "end", "identifier", "mean_signal", "numsamples", "summit", "core_start", "core_end", "component"]
    bed2.columns = ["chr",'start', 'end', 'MAP', 'r_value']

    # Convert DataFrame to BedTool
    bed1 = pybedtools.BedTool.from_dataframe(bed1)
    bed2 = pybedtools.BedTool.from_dataframe(bed2)

    # Get the intersection
    intersected = bed2.intersect(bed1, wa=True, wb=True)

    # Filter the columns
    df = pd.DataFrame(np.array(intersected.to_dataframe())[:, :5], columns= ['chr', 'start', 'end', 'MAP', 'r_value'])
    return df

def calculate_coverage(DF, map_value):
    return sum([1 for m in DF["MAP"] if m==map_value])/len(DF)

def custom_histplot(*args, **kwargs):
    # Extract data and variable from args
    data = kwargs.pop("data")
    variable = args[0]
    color = kwargs.pop("color")  # Extract color
    label = kwargs.pop("label")  # Extract label

    # Ensure the data is numeric and drop any NaN values
    numeric_data = pd.to_numeric(data[variable], errors='coerce').dropna()

    density = gaussian_kde(numeric_data)
    xs = np.linspace(np.min(numeric_data), np.max(numeric_data), 200)
    density_values = density(xs)

    # Normalize the density values so they sum up to 1
    probabilities = density_values / np.sum(density_values)

    plt.plot(xs, probabilities, color=color, label=label) 

    plt.xlabel("r_value")
    plt.ylabel("Probability")

def r_distribution_activeregions(r_value_file, cCREs_file, Meuleman_file, savedir):
    WG = pd.read_csv(r_value_file, sep="\t")

    Meuleman = get_overlap_Meuleman(Meuleman_file, r_value_file)

    cCREs = get_overlap_cCRE(cCREs_file, r_value_file)

    # Combine the dataframes into one for easier plotting
    WG['source'] = 'WG'
    Meuleman['source'] = 'Meuleman'
    cCREs['source'] = 'cCREs'
    combined = pd.concat([WG, Meuleman, cCREs])

    # Create a FacetGrid object
    g = sns.FacetGrid(combined, col="MAP", hue="source", col_wrap=4, sharey=False, legend_out=True)

    # Map a custom plot to each subplot
    g.map_dataframe(custom_histplot, "r_value")

    # Calculate coverage for each category and add as text to the plots
    for ax in g.axes.flat:
        map_val = ax.get_title().split('=')[-1].strip()  # Extract MAP value from title
        coverage_WG = calculate_coverage(WG, map_val)  # You need to define this function
        coverage_cCRE = calculate_coverage(cCREs, map_val)  # You need to define this function
        coverage_Meuleman = calculate_coverage(Meuleman, map_val)  # You need to define this function

        # Add text to the plot
        ax.text(0.02, 0.98, f'cvg_WG = {coverage_WG:.3f}', color='blue', transform=ax.transAxes, fontsize=8, verticalalignment='top')
        ax.text(0.02, 0.92, f'cvg_cCRE = {coverage_cCRE:.3f}', color='green', transform=ax.transAxes, fontsize=8, verticalalignment='top')
        ax.text(0.02, 0.86, f'cvg_Meuleman = {coverage_Meuleman:.3f}', color='orange', transform=ax.transAxes, fontsize=8, verticalalignment='top')

    # Add the legend after the plots are drawn
    g.add_legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3)

    # Show the plot
    plt.tight_layout()
    plt.savefig(f"{savedir}/rvalue_activeregions.pdf", format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f"{savedir}/rvalue_activeregions.svg", format='svg', dpi=300, bbox_inches='tight')

def custom_histplot2(*args, **kwargs):
    # Extract data and variable from args
    data = kwargs.pop("data")
    variable = args[0]
    color = kwargs.pop("color")  # Extract color
    label = kwargs.pop("label")  # Extract label
    coverage = kwargs.pop("coverage")  # Extract coverage
    ax = kwargs.pop("ax")

    # Ensure the data is numeric and drop any NaN values
    numeric_data = pd.to_numeric(data[variable], errors='coerce').dropna()

    density = gaussian_kde(numeric_data)
    xs = np.linspace(np.min(numeric_data), np.max(numeric_data), 200)
    density_values = density(xs)

    # Normalize the density values so they sum up to 1
    probabilities = density_values / np.sum(density_values)

    # Multiply the probabilities by the coverage
    probabilities *= coverage

    ax.plot(xs, probabilities, color=color, label=label) 

    ax.set_xlabel("r_value")
    ax.set_ylabel("Probability")

def r_distribution_activeregions2(r_value_file, r_value_cCREs, r_value_Meuleman, savedir):
    WG = pd.read_csv(r_value_file, sep="\t")
    Meuleman = pd.read_csv(r_value_Meuleman, sep="\t")
    cCREs = pd.read_csv(r_value_cCREs, sep="\t")

    # Combine the dataframes into one for easier plotting
    WG['source'] = 'WG'
    Meuleman['source'] = 'Meuleman'
    cCREs['source'] = 'cCREs'
    combined = pd.concat([WG, Meuleman, cCREs])

    # Create a FacetGrid object
    g = sns.FacetGrid(combined, col="MAP", hue="source", col_wrap=4, sharey=False)

    # Calculate coverage for each category and add as text to the plots
    for ax in g.axes.flat:
        map_val = ax.get_title().split('=')[-1].strip()  # Extract MAP value from title
        coverage_WG = calculate_coverage(WG, map_val)  # You need to define this function
        coverage_cCRE = calculate_coverage(cCREs, map_val)  # You need to define this function
        coverage_Meuleman = calculate_coverage(Meuleman, map_val)  # You need to define this function

        # Create a DataFrame for each source with only the rows for the current MAP value
        WG_current = WG[WG["MAP"] == map_val]
        cCREs_current = cCREs[cCREs["MAP"] == map_val]
        Meuleman_current = Meuleman[Meuleman["MAP"] == map_val]

        # Map a custom plot to each subplot for each source separately
        custom_histplot2("r_value", data=WG_current, color=g._colors[0], label='WG', coverage=coverage_WG, ax=ax)
        custom_histplot2("r_value", data=cCREs_current, color=g._colors[1], label='cCREs', coverage=coverage_cCRE, ax=ax)
        custom_histplot2("r_value", data=Meuleman_current, color=g._colors[2], label='Meuleman', coverage=coverage_Meuleman, ax=ax)

        # Add text to the plot
        ax.text(0.02, 0.98, f'cvg_WG = {coverage_WG:.3f}', color=g._colors[0], transform=ax.transAxes, fontsize=8, verticalalignment='top')
        ax.text(0.02, 0.92, f'cvg_cCRE = {coverage_cCRE:.3f}', color=g._colors[1], transform=ax.transAxes, fontsize=8, verticalalignment='top')
        ax.text(0.02, 0.86, f'cvg_Meuleman = {coverage_Meuleman:.3f}', color=g._colors[2], transform=ax.transAxes, fontsize=8, verticalalignment='top')

    # Create dummy lines for the legend
    line_wg = Line2D([0], [0], color=g._colors[0], lw=2)
    line_ccres = Line2D([0], [0], color=g._colors[1], lw=2)
    line_meuleman = Line2D([0], [0], color=g._colors[2], lw=2)

    # Add the legend after all plots are drawn
    g.add_legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3)

    # Show the plot
    plt.tight_layout()
    plt.savefig(f"{savedir}/rvalue_activeregions.pdf", format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f"{savedir}/rvalue_activeregions.svg", format='svg', dpi=300, bbox_inches='tight')

def r_dist_vs_expression(r_value_file, expression_file, savedir):
    gene_coords = load_gene_coords("src/biovalidation/parsed_genecode_data_hg38_release42.csv")
    r_vals = pd.read_csv(r_value_file, sep="\t")
    if "tsv" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=True)

    elif "pkl" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=False)
    
    # Convert DataFrame to BedTool
    bed1 = pybedtools.BedTool.from_dataframe(data)
    bed2 = pybedtools.BedTool.from_dataframe(r_vals)

    # Get the intersection
    df = bed2.intersect(bed1, wa=True, wb=True).to_dataframe()
    df = df[["chrom", "start", "end", "name", "score", "blockSizes"]]
    df.columns = ["chr", "start", "end", "MAP", "r_value", "TPM"]
    # df["TPM"] = np.log10(df["TPM"])

    unique_labels = df['MAP'].unique()
    n = len(unique_labels)

    # Find the factors of n that are closest to each other to get a 'square' layout
    sqrt_n = np.sqrt(n)
    cols = int(np.ceil(sqrt_n))
    rows = int(np.ceil(n / cols))

    fig, axs = plt.subplots(rows, cols, figsize=(4*cols, 4*rows), sharex=True, sharey=True) # Changed figsize here

    # Flatten the axes array if there's more than one row and column
    if rows > 1 and cols > 1:
        axs = axs.flatten()

    for i, label in enumerate(unique_labels):
        subset = df[df['MAP'] == label]
        sns.scatterplot(data=subset, x='r_value', y='TPM', ax=axs[i])
        sns.regplot(data=subset, x='r_value', y='TPM', ax=axs[i], scatter=False, color='red')
        correlation_coefficient = pearsonr(subset['r_value'], subset['TPM'])[0]
        axs[i].set_title(f'{label} | Pearson_r = {correlation_coefficient:.2f}')#, fontsize=9)

    # Remove any unused subplots
    if n < rows * cols:
        for i in range(n, rows * cols):
            fig.delaxes(axs[i])

    plt.tight_layout()
    plt.savefig(f"{savedir}/rvalue_v_expression.pdf", format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f"{savedir}/rvalue_v_expression.svg", format='svg', dpi=300, bbox_inches='tight')

def r_dist_vs_expression2(r_value_file, expression_file, savedir, n_bins=20):
    interpretation_terms = ["Prom", "Prom_fla", "Enha", "Enha_low", "Biva", "Tran", "Cons", "Facu", "K9K3", "Quie", "Unkn"]
    gene_coords = load_gene_coords("src/biovalidation/parsed_genecode_data_hg38_release42.csv")
    r_vals = pd.read_csv(r_value_file, sep="\t")
    if "tsv" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=True)

    elif "pkl" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=False)
    
    # Convert DataFrame to BedTool
    bed1 = pybedtools.BedTool.from_dataframe(data)
    bed2 = pybedtools.BedTool.from_dataframe(r_vals)

    # Get the intersection
    df = bed2.intersect(bed1, wa=True, wb=True).to_dataframe()
    df = df[["chrom", "start", "end", "name", "score", "blockSizes", "itemRgb"]]
    df.columns = ["chr", "start", "end", "MAP", "r_value", "TPM", "GeneID"]
    # df["TPM"] = np.log10(df["TPM"])

    unique_labels = df['MAP'].unique()
    binned_exp = {}

    expression_range = df["TPM"].max() - df["TPM"].min()
    for b in range(int(df["TPM"].min()), int(df["TPM"].max())+1, int(expression_range/n_bins)):
        bin_range = (b, b + int(expression_range/n_bins))

        subset = df.loc[
            (bin_range[0] < df["TPM"]) & 
            (df["TPM"] < bin_range[1]), :]

        if len(subset) < 1:
            continue

        listofgenes = []
        grouped = subset.groupby('GeneID')
        for _, group in grouped:
            listofgenes.append(group)

        binned_exp[bin_range] = listofgenes
    
    parsed = [] # each entry is a [GeneID, TPM_bin_center, exact_TPM, MAP, Fraction r>0.9, mean r_val]
    for b, l in binned_exp.items():
        for g in l:
            grouped = g.groupby('MAP')
            for _, group in grouped:
                expression_bin = (b[0]+b[1])/2
                expression_tpm = group["TPM"].unique()[0]
                robust = len(group.loc[(group["r_value"] >= 0.9), :]) / len(group)
                mean_r = group["r_value"].mean()
                seg_name = group["MAP"].unique()[0]
                seg_name = max([x for x in interpretation_terms if x in group["MAP"].unique()[0]], key=len)
                entry = [g["GeneID"].unique()[0], expression_bin, expression_tpm, seg_name, robust, mean_r]
                parsed.append(entry)

    # Convert parsed list to DataFrame for easier manipulation
    parsed_df = pd.DataFrame(parsed, columns=['GeneID', 'TPM_bin_center', "exact_TPM",  'MAP', 'Fraction_r>0.9', 'mean_r_val'])
    parsed_df = parsed_df.sort_values(by='MAP').reset_index(drop=True)

    ################################################################################################################
    unique_labels = parsed_df['MAP'].unique()
    n = len(unique_labels)

    # Find the factors of n that are closest to each other to get a 'square' layout
    sqrt_n = np.sqrt(n)
    cols = int(np.ceil(sqrt_n))
    rows = int(np.ceil(n / cols))

    fig, axs = plt.subplots(rows, cols, figsize=(3*cols, 3*rows), sharex=True, sharey=True) # Changed figsize here

    # Flatten the axes array if there's more than one row and column
    if rows > 1 and cols > 1:
        axs = axs.flatten()

    for i, label in enumerate(unique_labels):
        subset = parsed_df[parsed_df['MAP'] == label]
        sns.scatterplot(data=subset, x='Fraction_r>0.9', y='exact_TPM', ax=axs[i])
        sns.regplot(data=subset, x='Fraction_r>0.9', y='exact_TPM', ax=axs[i], scatter=False, color='red')
        correlation_coefficient = pearsonr(subset['Fraction_r>0.9'], subset['TPM_bin_center'])[0]
        axs[i].set_title(f'{label} | Pearson_r = {correlation_coefficient:.2f}')#, fontsize=9)

    # Remove any unused subplots
    if n < rows * cols:
        for i in range(n, rows * cols):
            fig.delaxes(axs[i])

    plt.tight_layout()
    plt.savefig(f"{savedir}/_robust_v_expression.pdf", format='pdf', dpi=150, bbox_inches='tight')
    plt.savefig(f"{savedir}/_robust_v_expression.svg", format='svg', dpi=150, bbox_inches='tight')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    fig, axs = plt.subplots(rows, cols, figsize=(3*cols, 3*rows), sharex=True, sharey=True) # Changed figsize here

    # Flatten the axes array if there's more than one row and column
    if rows > 1 and cols > 1:
        axs = axs.flatten()

    for i, label in enumerate(unique_labels):
        subset = parsed_df[parsed_df['MAP'] == label]
        sns.scatterplot(data=subset, x='mean_r_val', y='exact_TPM', ax=axs[i])
        sns.regplot(data=subset, x='mean_r_val', y='exact_TPM', ax=axs[i], scatter=False, color='red')
        correlation_coefficient = pearsonr(subset['mean_r_val'], subset['TPM_bin_center'])[0]
        axs[i].set_title(f'{label} | Pearson_r = {correlation_coefficient:.2f}')#, fontsize=9)

    # Remove any unused subplots
    if n < rows * cols:
        for i in range(n, rows * cols):
            fig.delaxes(axs[i])

    plt.tight_layout()
    plt.savefig(f"{savedir}/mean_r_v_expression.pdf", format='pdf', dpi=150, bbox_inches='tight')
    plt.savefig(f"{savedir}/mean_r_v_expression.svg", format='svg', dpi=150, bbox_inches='tight')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()
    ################################################################################################################

    # Create a color dictionary for each unique MAP value
    color_dict = {label: color for label, color in zip(parsed_df['MAP'].unique(), plt.cm.rainbow(np.linspace(0, 1, len(parsed_df['MAP'].unique()))))}

    # Plot 1: Fraction labeled as (e.g.) Trans with r > 0.9
    plt.figure(figsize=(10, 6))
    for label in parsed_df['MAP'].unique():
        subset = parsed_df[parsed_df['MAP'] == label]
        subset = subset.sort_values(by='TPM_bin_center')  # Sort by 'TPM_bin_center'
        
        # Calculate average and standard deviation for each bin
        avg = subset.groupby('TPM_bin_center')['Fraction_r>0.9'].mean()
        std = subset.groupby('TPM_bin_center')['Fraction_r>0.9'].std()
        
        plt.errorbar(avg.index.values, avg.values, fmt='o', color=color_dict[label], label=label)
        # plt.errorbar(avg.index.values, avg.values, yerr=std.values, fmt='o', color=color_dict[label], label=label)
        
        # Fit and plot a spline
        try:
            spline = UnivariateSpline(subset['TPM_bin_center'], subset['Fraction_r>0.9'], k=1)
            xs = np.linspace(subset['TPM_bin_center'].min(), subset['TPM_bin_center'].max(), n_bins)
            plt.plot(xs, spline(xs), color=color_dict[label])
        except:
            pass

    plt.xlabel('Gene Expression')
    plt.ylabel('Fraction labeled as Trans with r > 0.9')
    plt.legend()
    plt.savefig(r_value_file.replace("r_values.bed", "binned_exp_vs_ratio_robust.pdf"), format='pdf')
    plt.savefig(r_value_file.replace("r_values.bed", "binned_exp_vs_ratio_robust.svg"), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    # Plot 2: Mean Trans r-value
    plt.figure(figsize=(10, 6))
    for label in parsed_df['MAP'].unique():
        subset = parsed_df[parsed_df['MAP'] == label]
        subset = subset.sort_values(by='TPM_bin_center')  # Sort by 'TPM_bin_center'
        
        # Calculate average and standard deviation for each bin
        avg = subset.groupby('TPM_bin_center')['mean_r_val'].mean()
        std = subset.groupby('TPM_bin_center')['mean_r_val'].std()
        
        plt.errorbar(avg.index.values, avg.values,  fmt='o', color=color_dict[label], label=label)
        # plt.errorbar(avg.index.values, avg.values, yerr=std.values, fmt='o', color=color_dict[label], label=label)
        
        # Fit and plot a spline
        try:
            spline = UnivariateSpline(subset['TPM_bin_center'], subset['mean_r_val'], k=1)
            xs = np.linspace(subset['TPM_bin_center'].min(), subset['TPM_bin_center'].max(), n_bins)
            plt.plot(xs, spline(xs), color=color_dict[label])
        except:
            pass

    plt.xlabel('Gene Expression')
    plt.ylabel('Mean r-value')
    plt.legend()
    plt.savefig(r_value_file.replace("r_values.bed", "binned_exp_vs_mean_r.pdf"), format='pdf')
    plt.savefig(r_value_file.replace("r_values.bed", "binned_exp_vs_mean_r.svg"), format='svg')
    plt.clf()

    # Plot 1: Fraction labeled as (e.g.) Trans with r > 0.9
    plt.figure(figsize=(10, 6))
    for label in parsed_df['MAP'].unique():
        subset = parsed_df[parsed_df['MAP'] == label]
        subset = subset.sort_values(by='TPM_bin_center')  # Sort by 'TPM_bin_center'
        plt.scatter(subset['TPM_bin_center'], subset['Fraction_r>0.9'], color=color_dict[label], label=label, s=10)
        
        # Fit and plot a spline
        spline = UnivariateSpline(subset['TPM_bin_center'], subset['Fraction_r>0.9'], k=1)
        xs = np.linspace(subset['TPM_bin_center'].min(), subset['TPM_bin_center'].max(), n_bins)
        plt.plot(xs, spline(xs), color=color_dict[label])

    plt.xlabel('Gene Expression')
    plt.ylabel('Fraction labeled as Trans with r > 0.9')
    plt.legend()
    plt.savefig(r_value_file.replace("r_values.bed", "exp_vs_ratio_robust.pdf"), format='pdf')
    plt.savefig(r_value_file.replace("r_values.bed", "exp_vs_ratio_robust.svg"), format='svg')
    plt.clf()

    # Plot 2: Mean Trans r-value
    plt.figure(figsize=(10, 6))
    for label in parsed_df['MAP'].unique():
        subset = parsed_df[parsed_df['MAP'] == label]
        subset = subset.sort_values(by='TPM_bin_center')  # Sort by 'TPM_bin_center'
        plt.scatter(subset['TPM_bin_center'], subset['mean_r_val'], color=color_dict[label], label=label, s=10)
        
        # Fit and plot a spline
        spline = UnivariateSpline(subset['TPM_bin_center'], subset['mean_r_val'], k=1)
        xs = np.linspace(subset['TPM_bin_center'].min(), subset['TPM_bin_center'].max(), n_bins)
        plt.plot(xs, spline(xs), color=color_dict[label])

    plt.xlabel('Gene Expression')
    plt.ylabel('Mean r-value')
    plt.legend()
    plt.savefig(r_value_file.replace("r_values.bed", "exp_vs_mean_r.pdf"), format='pdf')
    plt.savefig(r_value_file.replace("r_values.bed", "exp_vs_mean_r.svg"), format='svg')
    plt.clf()

def r_dist_vs_expression3(r_value_file, expression_file, savedir, n_bins=20 , interpret=True):
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)

    interpretation_terms = ["Prom", "Prom_fla", "Enha", "Enha_low", "Biva", "Tran", "Cons", "Facu", "K9K3", "Quie", "Unkn"]
    gene_coords = load_gene_coords("src/biovalidation/parsed_genecode_data_hg38_release42.csv")
    r_vals = pd.read_csv(r_value_file, sep="\t")
    if "tsv" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=True)

    elif "pkl" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=False)
    
    # Convert DataFrame to BedTool
    bed1 = pybedtools.BedTool.from_dataframe(data)
    bed2 = pybedtools.BedTool.from_dataframe(r_vals)

    # Get the intersection
    df = bed2.intersect(bed1, wa=True, wb=True).to_dataframe()
    df2 = bed2.intersect(bed1, wa=True, wb=True, v=True).to_dataframe()
    df2["TPM"] = 0

    df = df[["chrom", "start", "end", "name", "score", "blockSizes"]]
    df.columns = ["chr", "start", "end", "MAP", "r_value", "TPM"]
    df2.columns = ["chr", "start", "end", "MAP", "r_value", "TPM"]

    # Merge df and df2
    df = pd.concat([df, df2])

    # Sort the resulting DataFrame by 'chr' and 'start'
    df = df.sort_values(by=['chr', 'start'])

    # df["TPM"] = np.log10(df["TPM"] + 1e-6)
    del df2

    unique_labels = df['MAP'].unique()
    binned_positions = {}

    r_range = 1
    for b in np.arange(0, 1, float(r_range/n_bins)):
        bin_range = (b, b + float(r_range/n_bins))

        subset = df.loc[
            (bin_range[0] < df["r_value"]) & 
            (df["r_value"] < bin_range[1]), :]

        if len(subset) < 1:
            continue

        perlabel = []
        grouped = subset.groupby('MAP')
        for _, group in grouped:
            perlabel.append(group)

        binned_positions[bin_range] = perlabel
    
    parsed = [] # each entry is a [r_bin_center, mean_r, MAP, mean_exp, non_zero_exp_fraction]
    for b, p in binned_positions.items():
        for l in p:
            r_bin_center = (b[0]+b[1])/2
            mean_r = l["r_value"].mean()
            # seg_name = l["MAP"].unique()[0]
            if interpret:
                seg_name_parts = l["MAP"].unique()[0].split("+")
                seg_name = "+".join([max([x for x in interpretation_terms if x in part], key=len) for part in seg_name_parts])
            else:
                seg_name = l["MAP"].unique()[0]

            mean_TPM = l["TPM"].mean()
            non_zero_exp_fraction = len(l.loc[(l["TPM"] > 0), :]) / len(l)

            entry = [r_bin_center, mean_r, seg_name, mean_TPM, non_zero_exp_fraction]
            parsed.append(entry)

    # Convert parsed list to DataFrame for easier manipulation
    parsed_df = pd.DataFrame(parsed, columns=['r_value', 'mean_r', "MAP",  'mean_exp', 'non_zero_exp_fraction'])
    parsed_df = parsed_df.sort_values(by='MAP').reset_index(drop=True)

    map_to_color = {label: color for label, color in zip(parsed_df['MAP'].unique(), plt.cm.rainbow(np.linspace(0, 1, len(parsed_df['MAP'].unique()))))}

    # Create a new figure for each plot
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    # For each unique MAP value
    for map_value in parsed_df['MAP'].unique():
        # Filter data for the current MAP value
        data = parsed_df[parsed_df['MAP'] == map_value]

        # Sort data by 'mean_r'
        data = data.sort_values('r_value')

        # Get the color for the current MAP value
        color = map_to_color[map_value]

        # Plot 1: mean_r vs mean_exp
        sns.regplot(x=data['r_value'], y=data['mean_exp'], scatter=True, color=color, label=map_value)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # plt.legend(loc='center', ncol=len(parsed_df['MAP'].unique()), bbox_to_anchor=(0.5, 0.5))
    
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.savefig(f"{savedir}/mean_expression_v_r.pdf", format='pdf')
    plt.savefig(f"{savedir}/mean_expression_v_r.svg", format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()
   
    ################################################################################################################
    ################################################################################################################    

    # Get the unique MAP values
    map_values = parsed_df['MAP'].unique()

    # Calculate the number of columns and rows for the subplots
    num_labels = len(map_values)
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    # Create a new figure with subplots. Adjust the figsize and layout as needed.
    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[16, 9])
    fig.tight_layout(pad=5.0)

    # Open a text file to store the metrics
    with open(f"{savedir}/mean_expression_v_r_metrics.txt", "w") as f:
        # For each unique MAP value
        for i, map_value in enumerate(map_values):
            # Filter data for the current MAP value
            data = parsed_df[parsed_df['MAP'] == map_value]

            # Sort data by 'mean_r'
            data = data.sort_values('r_value')

            # Fit a linear regression model and calculate R-squared value
            model = LinearRegression().fit(data[['r_value']], data['mean_exp'])
            r2 = model.score(data[['r_value']], data['mean_exp'])
            f.write(f"{map_value}|R2 = {r2}\n")

            try:
                # Calculate Pearson correlation coefficient
                pearson_r, _ = pearsonr(data['r_value'], data['mean_exp'])
                f.write(f"{map_value}|Pearson_r = {pearson_r}\n")
            except:
                pearson_r = 0
            
            try:
                # Calculate Spearman's rank correlation coefficient
                spearman_rho, _ = spearmanr(data['r_value'], data['mean_exp'])
                f.write(f"{map_value}|Spearman's rho = {spearman_rho}\n")
            except:
                pass

            try:
                # Calculate Kendall's tau
                kendall_tau, _ = kendalltau(data['r_value'], data['mean_exp'])
                f.write(f"{map_value}|Kendall's tau = {kendall_tau}\n")
            except:
                pass

            # Plot mean_r vs mean_exp in a subplot using sns.regplot
            ax = axs[i // n_cols, i % n_cols]
            sns.regplot(x=data['r_value'], y=data['mean_exp'], color=map_to_color[map_value], ax=ax)
            ax.set_title(f"{map_value} | R2: {r2:.2f} | Pearson_r: {pearson_r:.2f}")

    # Save the figure
    plt.ylim(bottom=0)
    plt.savefig(f"{savedir}/mean_expression_v_r_subplots.pdf", format='pdf')
    plt.savefig(f"{savedir}/mean_expression_v_r_subplots.svg", format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

def r_dist_vs_expression3_genebody(r_value_file, expression_file, savedir, n_bins=20, interpret=True):
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)

    interpretation_terms = ["Prom", "Prom_fla", "Enha", "Enha_low", "Biva", "Tran", "Cons", "Facu", "K9K3", "Quie", "Unkn"]
    gene_coords = load_gene_coords("src/biovalidation/parsed_genecode_data_hg38_release42.csv")
    r_vals = pd.read_csv(r_value_file, sep="\t")
    if "tsv" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=True)

    elif "pkl" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=False)
    
    # Convert DataFrame to BedTool
    bed1 = pybedtools.BedTool.from_dataframe(data)
    bed2 = pybedtools.BedTool.from_dataframe(r_vals)

    # Get the intersection
    df = bed2.intersect(bed1, wa=True, wb=True).to_dataframe()

    df = df[["chrom", "start", "end", "name", "score", "blockSizes"]]
    df.columns = ["chr", "start", "end", "MAP", "r_value", "TPM"]

    # Merge df and df2

    # Sort the resulting DataFrame by 'chr' and 'start'
    df = df.sort_values(by=['chr', 'start'])

    unique_labels = df['MAP'].unique()
    binned_positions = {}

    r_range = 1
    for b in np.arange(0, 1, float(r_range/n_bins)):
        bin_range = (b, b + float(r_range/n_bins))

        subset = df.loc[
            (bin_range[0] < df["r_value"]) & 
            (df["r_value"] < bin_range[1]), :]

        if len(subset) < 1:
            continue

        perlabel = []
        grouped = subset.groupby('MAP')
        for _, group in grouped:
            perlabel.append(group)

        binned_positions[bin_range] = perlabel
    
    parsed = [] # each entry is a [r_bin_center, mean_r, MAP, mean_exp, non_zero_exp_fraction]
    for b, p in binned_positions.items():
        for l in p:
            r_bin_center = (b[0]+b[1])/2
            mean_r = l["r_value"].mean()
            # seg_name = l["MAP"].unique()[0]
            if interpret:
                seg_name_parts = l["MAP"].unique()[0].split("+")
                seg_name = "+".join([max([x for x in interpretation_terms if x in part], key=len) for part in seg_name_parts])
            else:
                seg_name = l["MAP"].unique()[0]

            mean_TPM = l["TPM"].mean()
            non_zero_exp_fraction = len(l.loc[(l["TPM"] > 0), :]) / len(l)

            entry = [r_bin_center, mean_r, seg_name, mean_TPM, non_zero_exp_fraction]
            parsed.append(entry)

    # Convert parsed list to DataFrame for easier manipulation
    parsed_df = pd.DataFrame(parsed, columns=['r_value', 'mean_r', "MAP",  'mean_exp', 'non_zero_exp_fraction'])
    parsed_df = parsed_df.sort_values(by='MAP').reset_index(drop=True)

    map_to_color = {label: color for label, color in zip(parsed_df['MAP'].unique(), plt.cm.rainbow(np.linspace(0, 1, len(parsed_df['MAP'].unique()))))}

    # Create a new figure for each plot
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    # For each unique MAP value
    for map_value in parsed_df['MAP'].unique():
        # Filter data for the current MAP value
        data = parsed_df[parsed_df['MAP'] == map_value]

        # Sort data by 'mean_r'
        data = data.sort_values('r_value')

        # Get the color for the current MAP value
        color = map_to_color[map_value]

        # Plot 1: mean_r vs mean_exp
        sns.regplot(x=data['r_value'], y=data['mean_exp'], scatter=True, color=color, label=map_value)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # plt.legend(loc='center', ncol=len(parsed_df['MAP'].unique()), bbox_to_anchor=(0.5, 0.5))

    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.savefig(f"{savedir}/mean_expression_v_r_genebody.pdf", format='pdf')
    plt.savefig(f"{savedir}/mean_expression_v_r_genebody.svg", format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    ################################################################################################################
    ################################################################################################################    

    # Get the unique MAP values
    map_values = parsed_df['MAP'].unique()

    # Calculate the number of columns and rows for the subplots
    num_labels = len(map_values)
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    # Create a new figure with subplots. Adjust the figsize and layout as needed.
    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[16, 9])
    fig.tight_layout(pad=5.0)

    # Open a text file to store the metrics
    with open(f"{savedir}/mean_expression_v_r_genebody_metrics.txt", "w") as f:
        # For each unique MAP value
        for i, map_value in enumerate(map_values):
            # Filter data for the current MAP value
            data = parsed_df[parsed_df['MAP'] == map_value]

            # Sort data by 'mean_r'
            data = data.sort_values('r_value')

            # Fit a linear regression model and calculate R-squared value
            model = LinearRegression().fit(data[['r_value']], data['mean_exp'])
            r2 = model.score(data[['r_value']], data['mean_exp'])
            f.write(f"{map_value}|R2 = {r2}\n")

            try:
                # Calculate Pearson correlation coefficient
                pearson_r, _ = pearsonr(data['r_value'], data['mean_exp'])
                f.write(f"{map_value}|Pearson_r = {pearson_r}\n")
            except:
                pearson_r = 0
            
            try:
                # Calculate Spearman's rank correlation coefficient
                spearman_rho, _ = spearmanr(data['r_value'], data['mean_exp'])
                f.write(f"{map_value}|Spearman's rho = {spearman_rho}\n")
            except:
                pass

            try:
                # Calculate Kendall's tau
                kendall_tau, _ = kendalltau(data['r_value'], data['mean_exp'])
                f.write(f"{map_value}|Kendall's tau = {kendall_tau}\n")
            except:
                pass

            # Plot mean_r vs mean_exp in a subplot using sns.regplot
            ax = axs[i // n_cols, i % n_cols]
            sns.regplot(x=data['r_value'], y=data['mean_exp'], color=map_to_color[map_value], ax=ax)
            ax.set_title(f"{map_value} | R2: {r2:.2f} | Pearson_r: {pearson_r:.2f}")

    # Save the figure
    plt.ylim(bottom=0)
    plt.savefig(f"{savedir}/mean_expression_v_r_genebody_subplots.pdf", format='pdf')
    plt.savefig(f"{savedir}/mean_expression_v_r_genebody_subplots.svg", format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

def conf_v_nonconf_vs_expression(r_value_file, expression_file, savedir, n_bins=20 , interpret=True, alpha=0.9):
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)

    print("performing conf_v_nonconf_vs_expression")
    interpretation_terms = ["Prom", "Prom_fla", "Enha", "Enha_low", "Biva", "Tran", "Cons", "Facu", "K9K3", "Quie", "Unkn"]
    gene_coords = load_gene_coords("src/biovalidation/parsed_genecode_data_hg38_release42.csv")
    r_vals = pd.read_csv(r_value_file, sep="\t")
    if "tsv" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=True)

    elif "pkl" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=False)
    
    # Convert DataFrame to BedTool
    bed1 = pybedtools.BedTool.from_dataframe(data)
    bed2 = pybedtools.BedTool.from_dataframe(r_vals)

    # Get the intersection
    df = bed2.intersect(bed1, wa=True, wb=True).to_dataframe()
    df2 = bed2.intersect(bed1, wa=True, wb=True, v=True).to_dataframe()
    df2["TPM"] = 0

    df = df[["chrom", "start", "end", "name", "score", "blockSizes"]]
    df.columns = ["chr", "start", "end", "MAP", "r_value", "TPM"]
    df2.columns = ["chr", "start", "end", "MAP", "r_value", "TPM"]

    # Merge df and df2
    df = pd.concat([df, df2])

    # Sort the resulting DataFrame by 'chr' and 'start'
    df = df.sort_values(by=['chr', 'start'])
    del df2
    
    df_confident = df.loc[df["r_value"] >= alpha, :].reset_index(drop=True)
    df_non_confident = df.loc[df["r_value"] < alpha, :].reset_index(drop=True)

    def comparison_metrics(a, b):
        mean_diff = a.mean() - b.mean()
        try:
            t_stat, t_pval = ttest_ind(a, b)
        except:
            t_pval = np.nan 

        try:
            u_stat, u_pval = mannwhitneyu(a, b, alternative='two-sided')
        except:
            u_pval = np.nan 

        try:
            w, p = wilcoxon(a, b)
        except:
            p = np.nan 

        try:
            z_stat, zp_val = ztest(a, b, alternative='two-sided')
        except:
            zp_val = np.nan 

        return {"mean_diff":mean_diff, "t_pval":t_pval, "u_pval":u_pval, "wilcoxon":p, "zp_val":zp_val}

    map_values = [t for t in df['MAP'].unique() if "Tran" in t]
    # map_values = df['MAP'].unique()

    if len(map_values) > 1:
        print(map_values)
        # Calculate the number of columns and rows for the subplots
        num_labels = len(map_values)
        n_cols = math.floor(math.sqrt(num_labels))
        n_rows = math.ceil(num_labels / n_cols)

        # Create a new figure with subplots. Adjust the figsize and layout as needed.
        fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[n_rows*4, n_cols*3])
        fig.tight_layout(pad=5.0)

        agg_metrics = {}
        for i, label in enumerate(map_values):
            if axs.ndim == 1:
                ax = axs[i]
            else:
                ax = axs[i // n_cols, i % n_cols]

            df_label = df[df['MAP'] == label]
            df_confident_label = df_confident[df_confident['MAP'] == label]
            df_non_confident_label = df_non_confident[df_non_confident['MAP'] == label]

            # data_to_plot = pd.concat(
            #     [pd.Series(np.array(df_label['TPM'])), pd.Series(np.array(df_confident_label['TPM'])), pd.Series(np.array(df_non_confident_label['TPM']))], axis=1)
            # data_to_plot.columns = ['All', 'Confident', 'Non-confident']

            df_all = df_label[['TPM']].copy()
            df_all['Category'] = 'All'

            df_confident = df_confident_label[['TPM']].copy()
            df_confident['Category'] = 'Confident'

            df_non_confident = df_non_confident_label[['TPM']].copy()
            df_non_confident['Category'] = 'Non-confident'

            data_to_plot = pd.concat([df_all, df_confident, df_non_confident])

            sns.boxplot(x='Category', y='TPM', data=data_to_plot, palette=['grey', 'mediumaquamarine', 'lightcoral'], showfliers=False, ax=ax)
            
            try:
                sns.boxplot(data_to_plot, palette=['grey', 'mediumaquamarine', 'lightcoral'], showfliers=False, ax=ax)
            except:
                print(data_to_plot)

            metrics = comparison_metrics(df_confident_label['TPM'], df_non_confident_label['TPM'])

            for ii,kk in metrics.items():
                if ii not in agg_metrics:
                    agg_metrics[ii] = kk
                else:
                    agg_metrics[ii] = agg_metrics[ii] + kk

            title = label + f" | Z-test -log(p) = {(-1 * np.log( metrics['zp_val'] )):.2f}" 
            ax.set_title(title, fontsize=8)
            ax.set_ylabel("mean_expression")
        
        for ii in agg_metrics.keys():
            agg_metrics[ii] = agg_metrics[ii] / len(map_values)

    else:
        label = map_values[0]
        df_label = df[df['MAP'] == label]
        df_confident_label = df_confident[df_confident['MAP'] == label]
        df_non_confident_label = df_non_confident[df_non_confident['MAP'] == label]

        # data_to_plot = pd.concat(
        #     [pd.Series(np.array(df_label['TPM'])), pd.Series(np.array(df_confident_label['TPM'])), pd.Series(np.array(df_non_confident_label['TPM']))], axis=1)
        # data_to_plot.columns = ['All', 'Confident', 'Non-confident']
        
        df_all = df_label[['TPM']].copy()
        df_all['Category'] = 'All'

        df_confident = df_confident_label[['TPM']].copy()
        df_confident['Category'] = 'Confident'

        df_non_confident = df_non_confident_label[['TPM']].copy()
        df_non_confident['Category'] = 'Non-confident'

        data_to_plot = pd.concat([df_all, df_confident, df_non_confident])
        sns.boxplot(x='Category', y='TPM', data=data_to_plot, palette=['grey', 'mediumaquamarine', 'lightcoral'], showfliers=False)

        # sns.boxplot(data_to_plot, palette=['grey', 'mediumaquamarine', 'lightcoral'], showfliers=False)

        agg_metrics = comparison_metrics(df_confident_label['TPM'], df_non_confident_label['TPM'])
        title = label + f" | Z-test -log(p) = {(-1 * np.log( agg_metrics['zp_val'] )):.2f}" 
        plt.title(title)
        plt.ylabel("mean_expression")
    
    print("saving metrics in text format for conf_v_nonconf_vs_expression")
    with open(f"{savedir}/conf_vs_nonconf_meanEXP_metrics.txt", "w") as f:
        for ii, kk in agg_metrics.items():
            f.write(f"{ii}={kk}\n")

    print("creating plots for conf_v_nonconf_vs_expression")
    plt.tight_layout()
    plt.savefig(f"{savedir}/conf_vs_nonconf_meanEXP.pdf", format='pdf')
    plt.savefig(f"{savedir}/conf_vs_nonconf_meanEXP.svg", format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

def r_dist_vs_expression_boxplot(r_value_file, expression_file, savedir, n_bins=20 , interpret=True):
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)

    print("performing r_dist_vs_expression_boxplot")
    interpretation_terms = ["Prom", "Prom_fla", "Enha", "Enha_low", "Biva", "Tran", "Cons", "Facu", "K9K3", "Quie", "Unkn"]
    gene_coords = load_gene_coords("src/biovalidation/parsed_genecode_data_hg38_release42.csv")
    r_vals = pd.read_csv(r_value_file, sep="\t")
    if "tsv" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=True)

    elif "pkl" in expression_file:
        data = load_transcription_data(
            expression_file, 
            gene_coords, csv=False)
    
    # Convert DataFrame to BedTool
    bed1 = pybedtools.BedTool.from_dataframe(data)
    bed2 = pybedtools.BedTool.from_dataframe(r_vals)

    # Get the intersection
    df = bed2.intersect(bed1, wa=True, wb=True).to_dataframe()
    df2 = bed2.intersect(bed1, wa=True, wb=True, v=True).to_dataframe()
    df2["TPM"] = 0

    df = df[["chrom", "start", "end", "name", "score", "blockSizes"]]
    df.columns = ["chr", "start", "end", "MAP", "r_value", "TPM"]
    df2.columns = ["chr", "start", "end", "MAP", "r_value", "TPM"]

    # Merge df and df2
    df = pd.concat([df, df2])

    # Sort the resulting DataFrame by 'chr' and 'start'
    df = df.sort_values(by=['chr', 'start'])

    # df["TPM"] = np.log10(df["TPM"] + 1e-6)
    del df2

    unique_labels = df['MAP'].unique()
    binned_positions = {}

    r_range = 1
    for b in np.arange(0, 1, float(r_range/n_bins)):
        bin_range = (b, b + float(r_range/n_bins))

        subset = df.loc[
            (bin_range[0] < df["r_value"]) & 
            (df["r_value"] < bin_range[1]), :]

        if len(subset) < 1:
            continue

        perlabel = []
        grouped = subset.groupby('MAP')
        for _, group in grouped:
            perlabel.append(group)

        binned_positions[bin_range] = perlabel
    
    parsed = [] # each entry is a [r_bin_center, mean_r, MAP, mean_exp, non_zero_exp_fraction]
    for b, p in binned_positions.items():
        for l in p:
            r_bin_center = float(str(f"{((b[0]+b[1])/2):.3f}"))
            mean_r = l["r_value"].mean()
            # seg_name = l["MAP"].unique()[0]
            if interpret:
                seg_name_parts = l["MAP"].unique()[0].split("+")
                seg_name = "+".join([max([x for x in interpretation_terms if x in part], key=len) for part in seg_name_parts])
            else:
                seg_name = l["MAP"].unique()[0]

            mean_TPM = l["TPM"].mean()
            non_zero_exp_fraction = len(l.loc[(l["TPM"] > 0), :]) / len(l)

            entry = [r_bin_center, mean_r, seg_name, mean_TPM, non_zero_exp_fraction]
            parsed.append(entry)

    # Convert parsed list to DataFrame for easier manipulation
    parsed_df = pd.DataFrame(parsed, columns=['r_value', 'mean_r', "MAP",  'mean_exp', 'non_zero_exp_fraction'])
    parsed_df = parsed_df.sort_values(by='MAP').reset_index(drop=True)

    map_to_color = {label: color for label, color in zip(parsed_df['MAP'].unique(), plt.cm.rainbow(np.linspace(0, 1, len(parsed_df['MAP'].unique()))))} 

    # Get the unique MAP values
    map_values = parsed_df['MAP'].unique()

    # Calculate the number of columns and rows for the subplots
    num_labels = len(map_values)
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    # Create a new figure with subplots. Adjust the figsize and layout as needed.
    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[16, 9])
    fig.tight_layout(pad=5.0)

    # Open a text file to store the metrics
    for i, map_value in enumerate(map_values):
        # Filter data for the current MAP value
        data = parsed_df[parsed_df['MAP'] == map_value]

        # Sort data by 'mean_r'
        data = data.sort_values('r_value')

        # Plot mean_r vs mean_exp in a subplot using sns.regplot
        ax = axs[i // n_cols, i % n_cols]

        #### prompt: [[[ plot a boxplot here. particularly, for each bin, I want a boxplot with x="r_bin_center" and y="mean_exp" and the color should be according to map_value]]] ####
        sns.boxplot(x="r_value", y="mean_exp", data=data, ax=ax, color=map_to_color[map_value], showfliers=False)
        ax.set_title(f"{map_value}")

    # Save the figure
    print("making plots for r_dist_vs_expression_boxplot")
    plt.ylim(bottom=0)
    plt.savefig(f"{savedir}/mean_exp_boxplot_v_r_subplots.pdf", format='pdf')
    plt.savefig(f"{savedir}/mean_exp_boxplot_v_r_subplots.svg", format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

if __name__ == "__main__":
    # savedir = "tests/rebuttal_example/rebuttal_test_run/"
    # r_dist_vs_expression4("tests/rebuttal_example/rebuttal_test_run/r_values.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)
    # r_dist_vs_expression4("tests/rebuttal_example/rebuttal_test_run/r_values_14_states.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)
    # r_dist_vs_expression3("tests/rebuttal_example/rebuttal_test_run/r_values.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)
    # conf_v_nonconf_vs_expression("tests/rebuttal_example/rebuttal_test_run/r_values.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)
    # r_dist_vs_expression_boxplot("tests/rebuttal_example/rebuttal_test_run/r_values.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)
    # r_dist_vs_expression3_genebody("tests/rebuttal_example/rebuttal_test_run/r_values.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)
    # get_runs(maindir = "rebuttal", mp=True, n_processes=10)

    # get_runs(maindir = "rebuttal", mp=True, n_processes=15)
    get_runs(maindir = "rebuttal_subset", mp=False, n_processes=10, subset=True)
    exit()

    # savedir = "tests/rebuttal_example/rebuttal_test_run/"
    # r_dist_vs_expression3("tests/rebuttal_example/rebuttal_test_run/r_values_14_states.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)
    # r_dist_vs_expression3_genebody("tests/rebuttal_example/rebuttal_test_run/r_values_14_states.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)

    # replicate_1_dir = "tests/rebuttal_example/GM12878_R1/parsed_posterior.csv"
    # replicate_2_dir = "tests/rebuttal_example/GM12878_R2/parsed_posterior.csv"
    # savedir = "tests/rebuttal_example/rebuttal_test_run/"
    # base_mnemonics = "tests/rebuttal_example/GM12878_R1/mnemonics.txt"
    # verif_mnemonics = "tests/rebuttal_example/GM12878_R2/mnemonics.txt"

    # os.system(f"python SAGAconf.py --r_only -v -s -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
    # r_distribution_over_segment(f"{savedir}/r_values.bed", savedir)
    # os.system(f"python SAGAconf.py --v_seglength -v -s -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
    
    # os.system(f"python SAGAconf.py --active_regions -v -s -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
    # r_distribution_activeregions2(
    # f"{savedir}/r_values_WG.bed", 
    # f"{savedir}/r_values_cCRE.bed", 
    # f"{savedir}/r_values_muel.bed", 
    # savedir)

    # if "GM12878" in replicate_1_dir:
    #     expression_data = "src/biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv"
    # elif "MCF-7" in replicate_1_dir:
    #     expression_data = "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl"
    # elif "K562" in replicate_1_dir:
    #     expression_data = "src/biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv"
    # else:
    #     expression_data = "EMPTY"
    
    # if expression_data != "EMPTY":
    #     os.system(f"python SAGAconf.py --merge_only -k 14 -v -s -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
    #     os.system(f"python SAGAconf.py --merge_only -k 12 -v -s -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
    #     os.system(f"python SAGAconf.py --merge_only -k 10 -v -s -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")

    #     r_dist_vs_expression3(f"{savedir}/r_values.bed", expression_data, savedir+"/16_states/", interpret=True)
    #     r_dist_vs_expression3_genebody(f"{savedir}/r_values.bed", expression_data, savedir+"/16_states/", interpret=True)

    #     r_dist_vs_expression3(f"{savedir}/r_values_14_states.bed", expression_data, savedir+"/14_states/", interpret=True)
    #     r_dist_vs_expression3_genebody(f"{savedir}/r_values_14_states.bed", expression_data, savedir+"/14_states/", interpret=True)

    #     r_dist_vs_expression3(f"{savedir}/r_values_12_states.bed", expression_data, savedir+"/12_states/", interpret=True)
    #     r_dist_vs_expression3_genebody(f"{savedir}/r_values_12_states.bed", expression_data, savedir+"/12_states/", interpret=True)

    #     r_dist_vs_expression3(f"{savedir}/r_values_10_states.bed", expression_data, savedir+"/10_states/", interpret=True)
    #     r_dist_vs_expression3_genebody(f"{savedir}/r_values_10_states.bed", expression_data, savedir+"/10_states/", interpret=True)

    # os.system(f"python SAGAconf.py --merge_only -k 12 -v -s -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
    # os.system(f"python SAGAconf.py --v_seglength -v -s -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")
    # os.system(f"python SAGAconf.py --active_regions -v -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")

    # r_distribution_activeregions2(
    #     "tests/rebuttal_example/sagaconf_res_new/r_values_WG.bed", 
    #     "tests/rebuttal_example/sagaconf_res_new/r_values_cCRE.bed", 
    #     "tests/rebuttal_example/sagaconf_res_new/r_values_muel.bed", 
    #     savedir)
    
    # r_distribution_over_segment("tests/r_values.bed", savedir)
    # r_distribution_activeregions("tests/r_values.bed", "src/biointerpret/GRCh38-cCREs.bed", "src/biointerpret/Meuleman.tsv", savedir)
    
    # r_dist_vs_expression3("tests/rebuttal_example/sagaconf_res_new/r_values_12_states.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)
    # r_dist_vs_expression3_genebody("tests/rebuttal_example/sagaconf_res_new/r_values_12_states.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)
    # r_dist_vs_expression3_genebody("tests/r_values.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir, interpret=True)

    # r_dist_vs_expression2("tests/r_values.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir)
    # r_dist_vs_expression("tests/r_values.bed", "src/biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl", savedir)




