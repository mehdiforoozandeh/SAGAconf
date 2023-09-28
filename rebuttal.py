import os, pybedtools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import matplotlib.gridspec as gridspec


def get_listofruns(maindir="rebuttal"):
    listofruns = [
        {"replicate_1_dir":"chromhmm_runs/GM12878_rep1", 
        "replicate_2_dir":"chromhmm_runs/GM12878_rep2/", 
        "savedir":"{}/r1vsr2/chmm/MCF-7/".format(maindir)},
        
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

def get_runs(maindir = "rebuttal"):
    list_of_runs = get_listofruns(maindir)

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


    for r in list_of_runs:
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

        os.system(f"python SAGAconf.py --r_only -s -tr 0.9 -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")

def r_distribution_over_segment(r_value_file):
    interpretation_terms = ["Prom", "Prom_fla", "Enha", "Enha_low", "Biva", "Tran", "Cons", "Facu", "K9K3", "Quie", "Unkn"]
    r_vals = pd.read_csv(r_value_file, sep="\t").to_numpy()
    all_segs = {term: [] for term in interpretation_terms}

    resolution = r_vals[0, 2] - r_vals[0, 1]
    current_map = r_vals[0, 3]
    current_start = r_vals[0, 1]
    current_seg = []

    for i in range(len(r_vals)):
        if r_vals[i, 3] == current_map:
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

    splines_100_1k = {}
    splines_1k_10k = {}
    splines_10k_plus = {}

    for k in all_segs.keys():
        try:
            subset1 = np.concatenate([seg for seg in all_segs[k] if 100 < (len(seg)*resolution) <= 1000])
            sorted_indices = np.argsort(subset1[:, 0])
            subset1 = subset1[sorted_indices]
            x1 = subset1[:,0]
            y1 = subset1[:,1]
            splines_100_1k[k] = UnivariateSpline(x1, y1, k=5)

        except:
            pass

        try:
            subset2 = np.concatenate([seg for seg in all_segs[k] if 1000 < (len(seg)*resolution) <= 10000])
            sorted_indices = np.argsort(subset2[:, 0])
            subset2 = subset2[sorted_indices]
            x2 = subset2[:,0]
            y2 = subset2[:,1]
            splines_1k_10k[k] = UnivariateSpline(x2, y2, k=5)

        except:
            pass

        try:
            subset3 = np.concatenate([seg for seg in all_segs[k] if 10000 < (len(seg)*resolution)])
            sorted_indices = np.argsort(subset3[:, 0])
            subset3 = subset3[sorted_indices]
            x3 = subset3[:,0]
            y3 = subset3[:,1]
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
        # Generate x values
        x_values = np.linspace(0, 1, 100)

        try:
            # Subplot 1
            y_values = splines_100_1k[k](x_values)
            line, = ax0.plot(x_values, y_values, label=k, color=colors(i))
            lines.append(line)
            labels.append(k)
            ax0.set_title('Segments with length < 1kb')
            ax0.set_xlabel('Position relative to segment')
            ax0.set_ylabel('r_value')
        except:
            pass

        try:
            # Subplot 2
            y_values = splines_1k_10k[k](x_values)
            ax1.plot(x_values, y_values, label=k, color=colors(i))
            ax1.set_title('Segments with length 1kb - 10kb')
            ax1.set_xlabel('Position relative to segment')
            ax1.set_ylabel('r_value')
        except:
            pass

        try:
            # Subplot 3
            y_values = splines_10k_plus[k](x_values)
            ax2.plot(x_values, y_values, label=k, color=colors(i))
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
    plt.savefig(r_value_file.replace(".bed", "v_segment_length.pdf"), format='pdf')
    plt.savefig(r_value_file.replace(".bed", "v_segment_length.svg"), format='svg')

def get_all_r_distribution_over_segment(maindir="rebuttal"):
    list_of_runs = get_listofruns(maindir)

    for r in list_of_runs:
        r_val_file = r["savedir"] + "/r_value.bed"
        r_distribution_over_segment(r_val_file)

def get_overlap_with_activeregions(file1, file2):
    # Load the bed files
    a = pybedtools.BedTool(file1)
    b = pybedtools.BedTool(file2)

    # Intersect the bed files
    intersected = a.intersect(b, wa=True, wb=True)

    # Select the required columns
    df = intersected.to_dataframe()
    df = df[['chrom', 'start', 'end', 'name', 'score', 'itemRgb']]

    return df

if __name__ == "__main__":
    # get_runs(maindir = "rebuttal")
    get_all_r_distribution_over_segment(maindir = "rebuttal")
