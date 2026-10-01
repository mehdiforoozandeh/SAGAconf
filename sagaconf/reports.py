import math
import matplotlib.gridspec as gridspec
import numpy as np
import os
import pandas as pd
import pybedtools
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pybedtools import BedTool
from scipy.interpolate import UnivariateSpline
from sklearn import metrics
from ._cluster_matching import Hungarian_algorithm, IoU_overlap, connect_bipartite, joint_overlap_prob, overall_overlap_ratio, overlap_matrix
from ._reproducibility import NMI_from_matrix, joint_prob_MAP_with_posterior, joint_prob_with_binned_posterior, posterior_calibration
from ._utils import logit_array
from .granul import granularity_vs_agreement_nonsymmetric, merge_clusters
from .overall import condense_segments, is_repr_posterior, keep_reproducible_annotations, perlabel_is_reproduced, single_point_repr, write_MAPloci_in_BED
from .run import intersect_parsed_posteriors, read_mnemonics


def load_data(posterior1_dir, posterior2_dir, subset=False, logit_transform=False, force_WG=False):
    print("loading and intersecting")
    loci_1, loci_2 = intersect_parsed_posteriors(
        posterior1_dir, 
        posterior2_dir)

    if subset and force_WG==False:
        loci_1 = loci_1.loc[loci_1["chr"]=="chr21"].reset_index(drop=True)
        loci_2 = loci_2.loc[loci_2["chr"]=="chr21"].reset_index(drop=True)

    print("the shapes of the input matrices are: {}, {}".format(str(loci_1.shape), str(loci_2.shape)))

    if logit_transform:
        loci_1.iloc[:,3:] = logit_array(np.array(loci_1.iloc[:,3:]))
        loci_2.iloc[:,3:] = logit_array(np.array(loci_2.iloc[:,3:]))

    return loci_1, loci_2


def process_data(loci_1, loci_2, replicate_1_dir, replicate_2_dir, mnemons=True, vm="NA", bm="NA", match=False, custom_order=True):
    # print('generating confmat 1 ...')
    num_labels = loci_1.shape[1]-3

    loci_1.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    loci_2.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    
    if mnemons:
        print("loading mnemonics...")
        if os.path.exists(
            "/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics_rep1.txt") and os.path.exists(
            "/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics_rep2.txt"):

            print("reading concat mnemonics")
            loci_1_mnemon = read_mnemonics("/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics_rep1.txt")
            loci_2_mnemon = read_mnemonics("/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics_rep2.txt")
        else:
            print("reading mnemonics")
            if bm == "NA":
                bm = "/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics.txt"
            if vm == "NA":
                vm = "/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics.txt"

            loci_1_mnemon = read_mnemonics(bm)
            loci_2_mnemon = read_mnemonics(vm)

        mnemon1_dict = {}
        for i in loci_1_mnemon:
            if len(i.split("_")) == 2:
                mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]
            if len(i.split("_")) == 3:
                mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4] + '_' + i.split("_")[2][:3]

        mnemon2_dict = {}
        for i in loci_2_mnemon:
            if len(i.split("_")) == 2:
                mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]
            if len(i.split("_")) == 3:
                mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4] + '_' + i.split("_")[2][:3]

        #handle missing mnemonics
        for i in range(num_labels):
            if str(i) not in mnemon1_dict.keys():
                mnemon1_dict[str(i)] = str(f"{i}_Unkn")
            if str(i) not in mnemon2_dict.keys():
                mnemon2_dict[str(i)] = str(f"{i}_Unkn")
        
        if match:
            conf_mat = overlap_matrix(loci_1, loci_2, type="IoU")

            assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')
            for i in range(len(assignment_pairs)):
                assignment_pairs[i] = (mnemon1_dict[str(assignment_pairs[i][0])], mnemon2_dict[str(assignment_pairs[i][1])])
            print(assignment_pairs)



            loci_1, loci_2 = \
                connect_bipartite(loci_1, loci_2, assignment_pairs, mnemon=True)
            
            print('connected barpartite')

        else:
            loci_1.columns = list(loci_1.columns[:3]) + [mnemon1_dict[c.replace("posterior","")] for c in loci_1.columns[3:]]
            loci_2.columns = list(loci_2.columns[:3]) + [mnemon2_dict[c.replace("posterior","")] for c in loci_2.columns[3:]]

    else:
        if match:
            conf_mat = overlap_matrix(loci_1, loci_2, type="IoU")

            assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')
            loci_1, loci_2 = \
                connect_bipartite(loci_1, loci_2, assignment_pairs, mnemon=False)

            print('connected barpartite')
    
    if mnemons and custom_order:
        SORT_ORDER = {"Prom": 0, "Prom_fla":1, "Enha":2, "Enha_low":3, "Biva":4, "Tran":5, "Cons":6, "Facu":7, "K9K3":8, "Quie":9, "Unkn":10}
        try:
            new_columns = []
            for c in loci_1.columns[3:]:
                l = "_".join(c.split("_")[1:])
                new_columns.append(str(SORT_ORDER[l])+"_"+c)
                
            new_columns.sort()
            for i in range(len(new_columns)):
                new_columns[i] = new_columns[i][2:]

            loci_1 = loci_1[["chr", "start", "end"] + new_columns]
        except:
            pass

        ##########################################################################################
        ##########################################################################################
        try:
            new_columns = []
            for c in loci_2.columns[3:]:
                l = "_".join(c.split("_")[1:])
                new_columns.append(str(SORT_ORDER[l])+"_"+c)
                
            new_columns.sort()
            for i in range(len(new_columns)):
                new_columns[i] = new_columns[i][2:]

            loci_2 = loci_2[["chr", "start", "end"] + new_columns]
        except:
            pass

    return loci_1, loci_2


def convert_to_GenomeBrowser_viewable_BED(initial_rvalue_bed):
    r_vals = pd.read_csv(initial_rvalue_bed, sep="\t")
    r_vals.columns = ["chrom", "chromStart", "chromEnd", "name", "score"]
    r_vals["name"] = r_vals["name"].str.lower()
    LABEL_COLOR_MAP = {
        'prom': (1.0, 0.0, 0.0),
        'prom_fla': (1.0, 0.26666666666666666, 0.0),
        'enha': (1.0, 0.7647058823529411, 0.30196078431372547),
        'enha_low': (1.0, 1.0, 0.0),
        'biva': (0.7411764705882353, 0.7176470588235294, 0.4196078431372549),
        'ctcf': (0.7686274509803922, 0.8823529411764706, 0.0196078431372549),
        'tran': (0.0, 0.5019607843137255, 0.0),
        'k9k3': (0.4, 0.803921568627451, 0.6666666666666666),
        'facu': (0.5019607843137255, 0.0, 0.5019607843137255),
        'cons': (0.5411764705882353, 0.5686274509803921, 0.8156862745098039),
        'quie': (1.0, 1.0, 1.0),
        'Unkn': (0.0, 0.0, 0.0)
    }
    r_vals['itemRgb'] = r_vals['name'].apply(lambda x: LABEL_COLOR_MAP["_".join(x.split('_')[1:])])
    r_vals['itemRgb'] = r_vals['itemRgb'].apply(lambda x: ','.join([str(int(i*255)) for i in x]))
    r_vals.insert(5, 'strand', '.')
    r_vals['thickStart'] = r_vals['chromStart']
    r_vals['thickEnd'] = r_vals['chromEnd']
    r_vals = r_vals[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']]
    bed = BedTool.from_dataframe(r_vals)
    bed.saveas(initial_rvalue_bed.replace(".bed", "_UCSC_GenomeBrowser.bed"))


def subset_data_to_activeregions(
    replicate_1_dir, replicate_2_dir,
    cCREs_file="src/biointerpret/GRCh38-cCREs.bed",
    Meuleman_file="src/biointerpret/Meuleman.tsv", restrict_to="cCRE", locis=True):

    if locis:
        loci1, loci2 = replicate_1_dir, replicate_2_dir
    else:
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=False)

        loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)
    ##################################################################################################################

    bedloci1 = pybedtools.BedTool.from_dataframe(loci1)
    bedloci2 = pybedtools.BedTool.from_dataframe(loci2)

    if restrict_to == "cCRE":
        bed_ccre = pybedtools.BedTool(cCREs_file)

        # Get the intersection
        loci1_intersect_ccre = bedloci1.intersect(bed_ccre, wa=True, wb=False).to_dataframe()
        loci2_intersect_ccre = bedloci2.intersect(bed_ccre, wa=True, wb=False).to_dataframe()
        loci1_intersect_ccre.columns = loci1.columns
        loci2_intersect_ccre.columns = loci2.columns

        return loci1_intersect_ccre, loci2_intersect_ccre

    else:

        bed_meuleman = pd.read_csv(Meuleman_file, sep="\t")
        bed_meuleman.columns = ["chr", "start", "end", "identifier", "mean_signal", "numsamples", "summit", "core_start", "core_end", "component"]
        bed_meuleman = pybedtools.BedTool.from_dataframe(bed_meuleman)

        # Get the intersection
        loci1_intersect_meul = bedloci1.intersect(bed_meuleman, wa=True, wb=False).to_dataframe()
        loci2_intersect_meul = bedloci2.intersect(bed_meuleman, wa=True, wb=False).to_dataframe()

        loci1_intersect_meul.columns = loci1.columns
        loci2_intersect_meul.columns = loci2.columns

        return loci1_intersect_meul, loci2_intersect_meul


def get_rvals_activeregion(loci1, loci2, savedir, w=1000, restrict_to="cCRE"):
    to = 0.75
    rvalues = is_repr_posterior(
        loci1, loci2, ovr_threshold=to, window_bp=w, matching="static",
        always_include_best_match=True, return_r=True)

    rvalues.to_csv(savedir+f"/r_values_{restrict_to}.bed", sep='\t', header=True, index=False)


def ct_binned_posterior_heatmap(loci_1, loci_2, savedir, n_bins=10):
    indicator_file = "{}/binned_posterior_heatmap.txt".format(savedir)
    if os.path.exists(indicator_file):
        return

    loci_1.iloc[:,3:] = 1 / (1 + np.exp(-1 * np.array(loci_1.iloc[:,3:])))
    loci_2.iloc[:,3:] = 1 / (1 + np.exp(-1 * np.array(loci_2.iloc[:,3:])))

    num_labels = int(loci_1.shape[1]-3)
    matrix = joint_prob_with_binned_posterior(loci_1, loci_2, n_bins=n_bins, conditional=True, stratified=False)
    
    #create custom colormap
    boundaries = [x**2 for x in list(np.linspace(0, 1, 20))] + [1]# custom boundaries
    hex_colors = sns.light_palette('navy', n_colors=len(boundaries) * 2 + 2, as_cmap=False).as_hex()
    hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]
    colors=list(zip(boundaries, hex_colors))
    custom_color_map = LinearSegmentedColormap.from_list(
        name='custom_navy',
        colors=colors)
    
    p = sns.heatmap(
            matrix.astype(float), annot=False,
            linewidths=0.001,  cbar=True, 
            cmap=custom_color_map)

    sns.set(rc={'figure.figsize':(20,15)})

    yticks = list(matrix.index)
    p.set_yticks(np.arange(0, len(yticks), int(n_bins)) + 0.5)
    p.set_yticklabels([yticks[yt].split("|")[0] for yt in range(0, len(yticks), int(n_bins))], rotation=0)
    p.tick_params(axis='x', rotation=90, labelsize=9)
    p.tick_params(axis='y', rotation=0, labelsize=9)

    # plt.title('Overlap Ratio')
    plt.xlabel('Replicate 2 Labels')
    plt.ylabel("Replicate 1 Labels")
    plt.tight_layout()
    plt.savefig('{}/binned_posterior_heatmap.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/binned_posterior_heatmap.svg'.format(savedir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()


def ct_confus(loci_1, loci_2, savedir, w=1000):
    indicator_file = "{}/raw_conditional_overlap_ratio.txt".format(savedir)
    if os.path.exists(indicator_file):
        return

    """
    labels can be matched or not
    """
    #create custom colormap
    boundaries = [x**2 for x in list(np.linspace(0, 1, 20))] + [1] # custom boundaries
    hex_colors = sns.light_palette('navy', n_colors=len(boundaries) * 2 + 2, as_cmap=False).as_hex()
    hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]
    colors=list(zip(boundaries, hex_colors))
    custom_color_map = LinearSegmentedColormap.from_list(
        name='custom_navy',
        colors=colors)
    
    ####################################################################################
    confmat = IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False)

    p = sns.heatmap(
        confmat.astype(float), annot=True, fmt=".2f",
        linewidths=0.01,  cbar=True, annot_kws={"size": 8}, 
        vmin=0, vmax=1, cmap=custom_color_map)

    sns.set(rc={'figure.figsize':(20,15)})
    p.tick_params(axis='x', rotation=90, labelsize=10)
    p.tick_params(axis='y', rotation=0, labelsize=10)

    plt.title('IoU Overlap')
    plt.xlabel('Replicate 2 Labels')
    plt.ylabel("Replicate 1 Labels")
    plt.tight_layout()
    plt.savefig('{}/heatmap.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/heatmap.svg'.format(savedir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    confmat.to_csv("{}/heatmap.csv".format(savedir))

    conditional = overlap_matrix(loci_1, loci_2, type="conditional")
    with open("{}/overlap_ratio.txt".format(savedir), "w") as cf:
        overall_overlap = overall_overlap_ratio(loci_1, loci_2, w=0)
        cf.write("{} : {}\n".format("naive overall w=0", overall_overlap))

        overall_overlap = overall_overlap_ratio(loci_1, loci_2, w=w)
        cf.write("{} : {}\n".format(f"overall w={w}", overall_overlap))

        for i in conditional.index:
            cf.write("{} : {}\n".format(i, np.max(np.array(conditional.loc[i, :]))))

    ####################################################################################

    confmat = IoU_overlap(loci_1, loci_2, w=w, symmetric=False, soft=False)
    p = sns.heatmap(
        confmat.astype(float), annot=True, fmt=".2f",
        linewidths=0.01,  cbar=True, annot_kws={"size": 8}, 
        vmin=0, vmax=1, cmap=custom_color_map)

    sns.set(rc={'figure.figsize':(20,15)})
    p.tick_params(axis='x', rotation=90, labelsize=10)
    p.tick_params(axis='y', rotation=0, labelsize=10)

    plt.title('IoU Overlap | w={}'.format(w))
    plt.xlabel('Replicate 2 Labels')
    plt.ylabel("Replicate 1 Labels")
    plt.tight_layout()
    plt.savefig('{}/heatmap_w.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/heatmap_w.svg'.format(savedir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    confmat.to_csv("{}/heatmap_w.csv".format(savedir))


def ct_granul(loci_1, loci_2, savedir):
    indicator_file = savedir+"/granularity.pdf"
    if os.path.exists(indicator_file):
        return

    """
    for this function, the labels should not be matched
    """
    
    num_labels = loci_1.shape[1]-3
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    p_to_r_auc_record = {}

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[25, 16])
    label_being_plotted = 0

    for i in range(n_rows):
        for j in range(n_cols):
            if label_being_plotted < num_labels:
                c = list(loci_1.columns[3:])[label_being_plotted]
                cr, ar, ovr_rec = granularity_vs_agreement_nonsymmetric(loci_1.copy(), loci_2.copy(), k=c)

                perfect_agr = [0] + [1 for i in range(len(ar) - 1)]
                realAUC = metrics.auc(cr, ar)
                perfectAUC = metrics.auc(cr, perfect_agr)
                p_to_r_auc = float((realAUC)/(perfectAUC))

                p_to_r_auc_record[c] = p_to_r_auc

                axs[i,j].plot(cr, ar, c="yellowgreen")
                axs[i,j].set_title(
                    c + str(" | Real/Perfect AUC = {:.2f}".format(p_to_r_auc)), 
                    fontsize=15)

                axs[i,j].fill_between(cr, ar, color="yellowgreen", alpha=0.4)
                axs[i,j].fill_between(cr, perfect_agr, ar, color="palevioletred", alpha=0.4)
                axs[i,j].set_xticks(np.arange(0, 1.1, step=0.2))
                axs[i,j].set_yticks(np.arange(0, 1.1, step=0.2))
                axs[i,j].tick_params(axis='both', which='major', labelsize=15) 
                label_being_plotted+=1
        
    # fig.text(0.5, 0.02, 'Coverage' , ha='center')
    # fig.text(0.02, 0.5, 'Agreement', va='center', rotation='vertical')
    plt.tight_layout()
    plt.savefig(savedir+"/granularity.pdf", format='pdf')
    plt.savefig(savedir+"/granularity.svg", format='svg')

    plt.clf()
    sns.reset_orig
    plt.style.use('default')
    plt.close("all")

    fig, ax = plt.subplots(figsize=(12, 9))
    ax.bar(p_to_r_auc_record.keys(), p_to_r_auc_record.values(), color='black', alpha=0.5)
    ax.set_xlabel('Chromatin States')
    ax.set_ylabel('Real/Perfect AUC')
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    ax.tick_params(axis='x', rotation=30, labelsize=12)
    ax.tick_params(axis='y', rotation=0, labelsize=12)
    ax.axhline(y=0.5, color='r', linestyle='--', linewidth=3)

    with open(savedir+"/AUC_mAUC.txt", 'w') as f:
        f.write(str(p_to_r_auc_record))

    plt.savefig(savedir+"/barplot.pdf", format='pdf')
    plt.savefig(savedir+"/barplot.svg", format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()


def ct_lable_calib(loci_1, loci_2, pltsavedir):
    indicator_file = pltsavedir+"/calib_logit_enr"
    if os.path.exists(indicator_file):
        return

    """
    labels need to be matched
    """
    if os.path.exists(pltsavedir+"/calib") == False:
        os.mkdir(pltsavedir+"/calib")
    calb = posterior_calibration(
        loci_1, loci_2, window_size=1000, savedir=pltsavedir+"/calib", allow_w=False)
    calibrated_loci_1 = calb.perlabel_calibration_function(return_caliberated_matrix=False)

    # if os.path.exists(pltsavedir+"/calib_logit_enr") == False:
    #     os.mkdir(pltsavedir+"/calib_logit_enr")
    # calb = posterior_calibration(
    #     loci_1, loci_2, plot_raw=False, window_size=1000, savedir=pltsavedir+"/calib_logit_enr", allow_w=False)
    # calibrated_loci_1 = calb.perlabel_calibration_function()
    
    plt.close("all")
    plt.style.use('default')


def overall_boundary(loci_1, loci_2, savedir, match_definition="BM"):
    indicator_file = savedir+"/len_bound_{}.pdf".format("overall")
    if os.path.exists(indicator_file):
        return

    """
    a dict for match definition
    for w in range(max_distance):
        check the overall correspondence within window of size w
    
    plot like normal boundary
    """
    resolution = int(loci_1.iloc[0, 2] - loci_1.iloc[0, 1])

    max_distance = int(5000/resolution)
    num_labels = loci_1.shape[1]-3
    MAP1 = list(loci_1.iloc[:,3:].idxmax(axis=1))
    MAP2 = list(loci_2.iloc[:,3:].idxmax(axis=1))

    confmat = IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False)
    
    # define matches
    per_label_matches = {}
    for k in list(loci_1.columns[3:]):
        sorted_k_vector = confmat.loc[k,:].sort_values(ascending=False)

        good_matches = sorted_k_vector.index[0]
        per_label_matches[k] = good_matches
    #========================================================================================#
    listofws = []
    overlaprecord = []
    
    for w in range(max_distance):
        m = 0
        for i in range(len(MAP1)):
            k = MAP1[i]
            if w == 0:
                i_neighbors = [MAP2[i]]
            else:
                i_neighbors = MAP2[max(0, i-w) : min(i+w, len(MAP2)-1)]

            if per_label_matches[k] in i_neighbors:
                m += 1
        
        listofws.append(w)
        overlaprecord.append(float(m) / len(loci_1))
    
    plt.plot([x * resolution for x in listofws], overlaprecord, color="black", label="distance from boundary", linewidth=2)
    plt.xticks(np.arange(0, (max_distance+1)*resolution, step=5*resolution))
    plt.tick_params(axis='both', labelsize=12)
    plt.tick_params(axis='x', rotation=90)
    plt.yticks(np.arange(0, 1.1, step=0.1))
    plt.xlabel('bp', fontsize=10)
    plt.ylabel('Ratio Overlap',fontsize=10)
    
    plt.title("Overall", fontsize=12)
    plt.tight_layout()
    plt.savefig(savedir+"/len_bound_{}.pdf".format("overall"), format='pdf')
    plt.savefig(savedir+"/len_bound_{}.svg".format("overall"), format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()


def distance_vs_overlap(loci_1, loci_2, savedir, match_definition="BM"):
    indicator_file = savedir+"/Dist_vs_Corresp"
    if os.path.exists(indicator_file):
        return

    if os.path.exists(savedir+"/Dist_vs_Corresp") == False:
        os.mkdir(savedir+"/Dist_vs_Corresp")
    savedir = savedir+"/Dist_vs_Corresp"

    resolution = int(loci_1.iloc[0, 2] - loci_1.iloc[0, 1])

    max_distance = int(5000/resolution)
    num_labels = loci_1.shape[1]-3
    MAP1 = list(loci_1.iloc[:,3:].idxmax(axis=1))
    MAP2 = list(loci_2.iloc[:,3:].idxmax(axis=1))

    confmat = IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False)
    
    # define matches
    per_label_matches = {}
    for k in list(loci_1.columns[3:]):
        sorted_k_vector = confmat.loc[k,:].sort_values(ascending=False)

        good_matches = sorted_k_vector.index[0]
        per_label_matches[k] = good_matches
    #========================================================================================#
    distance_to_corresp = []
    for i in range(len(MAP1)):
        k = MAP1[i]

        matched = False
        for w in range(max_distance):
            if matched==False and w == 0:
                if MAP2[i] == per_label_matches[k]:
                    distance_to_corresp.append(0)
                    matched = True

            elif matched==False and w > 0:
                upst_neighbors = MAP2[i : min(i+w+1, len(MAP2))]
                downst_neighbors = MAP2[max(0, i-w) : i+1]
                if per_label_matches[k] in upst_neighbors:
                    distance_to_corresp.append(w)
                    matched = True
                
                elif per_label_matches[k] in downst_neighbors:
                    distance_to_corresp.append(-1*w)
                    matched = True

            if matched:
                break

            elif matched==False and w==(max_distance-1):
                distance_to_corresp.append(None)
    
    nonefiltered = [x * resolution for x in distance_to_corresp if x is not None] 
    matched_ratio = len(nonefiltered) / len(distance_to_corresp)
    plt.hist(nonefiltered, bins=len(set(nonefiltered)), density=True, log=True, color='black', alpha=0.6, histtype="stepfilled")
    plt.axvline(x=0, color='red', linestyle='dotted', linewidth=1.5)
    plt.yticks(np.logspace(-6, 0, 7))
    plt.xlabel("Distance (bp)")
    plt.ylabel("Matched label density (log scale)")
    plt.title("Correspondence vs. Distance -- Overall | overlap ratio = {:.2f}".format(matched_ratio))
    plt.tight_layout()
    plt.savefig(savedir+"/dist_vs_corresp_{}.pdf".format("overall"), format='pdf')
    plt.savefig(savedir+"/dist_vs_corresp_{}.svg".format("overall"), format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()


    #========================================================================================#
    perlabel = {}
    for k in loci_1.columns[3:]:
        perlabel[k] = []
    
    for i in range(len(MAP1)):
        perlabel[MAP1[i]].append(distance_to_corresp[i])

    #========================================================================================#
    for k in perlabel.keys():
        distance_to_corresp_k = perlabel[k]
        nonefiltered = [x * resolution for x in distance_to_corresp_k if x is not None] 
        matched_ratio = len(nonefiltered) / len(distance_to_corresp_k)

        plt.hist(nonefiltered, bins=len(set(nonefiltered)), density=True, log=True, color='black', alpha=0.6, histtype="stepfilled")
        plt.axvline(x=0, color='red', linestyle='dotted', linewidth=1.5)
        plt.yticks(np.logspace(-6, 0, 7))
        plt.xlabel("Distance (bp)")
        plt.ylabel("Matched label density (log scale)")
        plt.title("Correspondence vs. Distance -- {} | overlap ratio = {:.2f}".format(k, matched_ratio))
        plt.tight_layout()
        plt.savefig(savedir+"/dist_vs_corresp_{}.pdf".format(k), format='pdf')
        plt.savefig(savedir+"/dist_vs_corresp_{}.svg".format(k), format='svg')

        sns.reset_orig
        plt.close("all")
        plt.style.use('default')
        plt.clf()


    #========================================================================================#
    num_labels = loci_1.shape[1]-3
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[25, 16])
    label_being_plotted = 0
    
    for i in range(n_rows):
        for j in range(n_cols):
            if label_being_plotted < num_labels:
                k = loci_1.columns[3:][label_being_plotted]
                distance_to_corresp_k = perlabel[k]
                nonefiltered = [x * resolution for x in distance_to_corresp_k if x is not None] 
                matched_ratio = len(nonefiltered) / len(distance_to_corresp_k)

                axs[i,j].hist(
                    nonefiltered, bins=len(set(nonefiltered)), density=True, log=True, 
                    color='black', alpha=0.6, histtype="stepfilled")

                axs[i,j].axvline(x=0, color='red', linestyle='dotted', linewidth=2)
                axs[i,j].set_yticks(np.logspace(-6, 0, 7))
                axs[i,j].set_title("{} | overlap ratio = {:.2f}".format(k, matched_ratio), fontsize=15)
                axs[i,j].tick_params(axis='both', which='major', labelsize=15) 
                label_being_plotted += 1
        
    plt.tight_layout()
    plt.savefig(savedir+"/dist_vs_corresp_{}.pdf".format("subplot"), format='pdf')
    plt.savefig(savedir+"/dist_vs_corresp_{}.svg".format("subplot"), format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()


def distance_vs_overlap_3(loci_1, loci_2, savedir, match_definition="BM"):
    indicator_file = savedir+"/Dist_vs_Corresp_3"
    if os.path.exists(indicator_file):
        return
        
    if os.path.exists(savedir+"/Dist_vs_Corresp_3") == False:
        os.mkdir(savedir+"/Dist_vs_Corresp_3")
    savedir = savedir+"/Dist_vs_Corresp_3"

    resolution = int(loci_1.iloc[0, 2] - loci_1.iloc[0, 1])

    max_distance = int(5000/resolution)
    num_labels = loci_1.shape[1]-3

    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    coverage1 = {k:len(MAP1.loc[MAP1==k])/len(loci_1) for k in loci_1.columns[3:]}

    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
    coverage2 = {k:len(MAP2.loc[MAP2==k])/len(loci_2) for k in loci_2.columns[3:]}

    MAP1 = list(MAP1)
    MAP2 = list(MAP2)

    confmat = IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False)
    
    # define matches
    per_label_matches = {}
    for k in list(loci_1.columns[3:]):
        sorted_k_vector = confmat.loc[k,:].sort_values(ascending=False)

        good_matches = sorted_k_vector.index[0]
        per_label_matches[k] = good_matches
    
    # x:[[0, 0],[0, 0]] -> [[right_match, right_mismatch]], [left_match, left_mismatch]]
    label_dist_dict = {
        l : {x:[[0, 0],[0, 0]] for x in range(max_distance)} for l in loci_1.columns[3:]
    }

    for i in range(loci_1.shape[0]):
        for w in range(max_distance):
            l = MAP1[i]

            if MAP2[min([len(MAP2)-1, i + w])] == per_label_matches[l]:
                label_dist_dict[l][w][0][0] += 1
            else:
                label_dist_dict[l][w][0][1] += 1

            if MAP2[max([0, i - w])] == per_label_matches[l]:
                label_dist_dict[l][w][1][0] += 1
            else:
                label_dist_dict[l][w][1][1] += 1
    
    for l in label_dist_dict.keys():
        for w in label_dist_dict[l].keys():
            left_prob = (label_dist_dict[l][w][1][0]) / (label_dist_dict[l][w][1][0] + label_dist_dict[l][w][1][1])
            right_prob = (label_dist_dict[l][w][0][0]) / (label_dist_dict[l][w][0][0] + label_dist_dict[l][w][0][1])

            label_dist_dict[l][w] = [left_prob, right_prob]


    for l in label_dist_dict.keys():
        dist_to_corresp = {}

        for w in label_dist_dict[l].keys():
            dist_to_corresp[w] = label_dist_dict[l][w][1]
            dist_to_corresp[-1 * w] = label_dist_dict[l][w][0]

        label_dist_dict[l] = dist_to_corresp
        
    label_dist_dict = pd.DataFrame(label_dist_dict).sort_index()

    ################################################################################################
    num_labels = loci_1.shape[1]-3
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[25, 16])
    label_being_plotted = 0
    
    for i in range(n_rows):
        for j in range(n_cols):
            if label_being_plotted < num_labels:
                k = label_dist_dict.columns[label_being_plotted]

                xaxis = [x*resolution for x in label_dist_dict.index]
                yaxis = list(label_dist_dict[k])

                axs[i,j].plot(xaxis, yaxis, color="black")

                axs[i,j].fill_between(xaxis, yaxis, color="black", alpha=0.4)
                
                # axs[i,j].axhline(y=coverage2[per_label_matches[k]], color='red', linestyle='--', linewidth=1.5)
                axs[i,j].axvline(x=0, color='blue', linestyle='--', linewidth=1.5)

                axs[i,j].set_xlabel("Distance (bp)")
                axs[i,j].set_ylabel("Probability of overlap with corresponding label")
                axs[i,j].set_title(k)

                label_being_plotted += 1
        
    plt.tight_layout()
    plt.savefig(savedir+"/dist_vs_corresp_{}.pdf".format("subplot"), format='pdf')
    plt.savefig(savedir+"/dist_vs_corresp_{}.svg".format("subplot"), format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()


def ct_boundar(loci_1, loci_2, outdir, match_definition="BM", max_distance=50):
    indicator_file = outdir+"/len_bound.pdf"
    if os.path.exists(indicator_file):
        return

    """
    Here i'm gonna merge ECDF of length dist and boundary thing.
    """

    resolution = int(loci_1.iloc[0, 2] - loci_1.iloc[0, 1])

    max_distance = int(5000/resolution)
    num_labels = loci_1.shape[1]-3
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

    coverage1 = {k:len(MAP1.loc[MAP1 == k])/len(loci_1) for k in loci_1.columns[3:]}
    coverage2 = {k:len(MAP2.loc[MAP2 == k])/len(loci_2) for k in loci_2.columns[3:]}

    with open(outdir+"/coverages1.txt", 'w') as coveragesfile:
        coveragesfile.write(str(coverage1))

    with open(outdir+"/coverages2.txt", 'w') as coveragesfile:
        coveragesfile.write(str(coverage2))
    
    """
    we need to define what a "good match" is.
    by default, we will consider any label with log(OO/EO)>0 
    to be a match. this definition can be refined later.
    """

    confmat = overlap_matrix(loci_1, loci_2, type="IoU")
    
    # define matches
    per_label_matches = {}
    for k in list(loci_1.columns[3:]):
        sorted_k_vector = confmat.loc[k,:].sort_values(ascending=False)

        if match_definition=="EO" or match_definition=="enrichment_of_overlap":
            good_matches = [sorted_k_vector.index[j] for j in range(len(sorted_k_vector)) if sorted_k_vector[j]>0]
            per_label_matches[k] = good_matches

        elif match_definition=="BM" or match_definition=="best_match":
            good_matches = [sorted_k_vector.index[0]]
            per_label_matches[k] = good_matches

    #========================================================================================#
    # now define a boundary
    boundary_distances = {}
    count_unmatched_boundaries = {}
    for k in list(loci_1.columns[3:]):
        boundary_distances[k] = []
        count_unmatched_boundaries[k] = 0

    for i in range(len(MAP1)-1):
        # if MAP1[i] != MAP1[i+1]: 
        #this means transition
        """
        now check how far should we go from point i in the other replicate to see a good match.
        """
        if MAP2[i] not in per_label_matches[MAP1[i]]:
            matched = False
            dist = 0
            while matched == False and dist<max_distance:
                if i + (-1*dist) >= 0:
                    if MAP2[i + (-1*dist)] in per_label_matches[MAP1[i]]:
                        matched = True

                if i + dist < len(MAP2):
                    if MAP2[i + dist] in per_label_matches[MAP1[i]]:
                        matched = True

                if matched==False:
                    dist += 1

            if matched==True:
                boundary_distances[MAP1[i]].append(dist)

            else:  
                count_unmatched_boundaries[MAP1[i]] += 1

        else:
            boundary_distances[MAP1[i]].append(0)

    #========================================================================================#
    # plot histograms

    num_labels = loci_1.shape[1]-3
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[16, 12])
    label_being_plotted = 0
    
    for i in range(n_rows):
        for j in range(n_cols):
            if label_being_plotted < num_labels:
                k = loci_1.columns[3:][label_being_plotted]

                matched_hist = np.histogram(boundary_distances[k], bins=max_distance, range=(0, max_distance)) #[0]is the bin size [1] is the bin edge

                cdf = []
                for jb in range(len(matched_hist[0])):
                    if len(cdf) == 0:
                        cdf.append(matched_hist[0][jb])
                    else:
                        cdf.append((cdf[-1] + matched_hist[0][jb]))

                cdf = np.array(cdf) / (np.sum(matched_hist[0]) + count_unmatched_boundaries[k])
                
                
                axs[i, j].plot(list(matched_hist[1][:-1]*resolution), list(cdf), color="black", label="distance from boundary")
                axs[i, j].fill_between(list(matched_hist[1][:-1]*resolution), list(cdf), color="black", alpha=0.35)
                
                # lendist = get_ecdf_for_label(loci_1, k, max=(max_distance+1)*resolution)
                # if len(lendist)>0:
                #     sns.ecdfplot(lendist, ax=axs[i, j], color="red")

                axs[i, j].set_xticks(np.arange(0, (max_distance+1)*resolution, step=5*resolution))
                axs[i, j].tick_params(axis='both', labelsize=11)
                axs[i, j].tick_params(axis='x', rotation=90)
                axs[i, j].set_yticks(np.arange(0, 1.1, step=0.2))
                axs[i, j].set_xlabel('bp', fontsize=13)
                axs[i, j].set_ylabel('overlap ratio',fontsize=13)
                
                
                axs[i, j].set_title(k, fontsize=14)

                label_being_plotted += 1
            
    plt.tight_layout()
    plt.savefig(outdir+"/len_bound.pdf", format='pdf')
    plt.savefig(outdir+"/len_bound.svg", format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()


def get_all_ct(replicate_1_dir, replicate_2_dir, savedir, locis=False, w=1000):
    if locis:
        loci1, loci2 = replicate_1_dir, replicate_2_dir
    else:
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=True)

        loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    ct_binned_posterior_heatmap(loci1.copy(), loci2.copy(), savedir)

    ct_confus(loci1.copy(), loci2.copy(), savedir, w=w)

    ct_lable_calib(loci1.copy(), loci2.copy(), savedir)

    ct_granul(loci1.copy(), loci2.copy(), savedir)

    ct_boundar(loci1.copy(), loci2.copy(), savedir, match_definition="BM", max_distance=50)

    overall_boundary(loci1.copy(), loci2.copy(), savedir, match_definition="BM")

    distance_vs_overlap(loci1.copy(), loci2.copy(), savedir, match_definition="BM")
    
    distance_vs_overlap_3(loci1.copy(), loci2.copy(), savedir, match_definition="BM")


def get_overalls(replicate_1_dir, replicate_2_dir, savedir, locis=False, w=1000, to=0.75, tr=0.9):
    indicator_file = savedir+"/boolean_reproducibility_report_POSTERIOR.txt"
    if os.path.exists(indicator_file):
        return

    if locis:
        loci1, loci2 = replicate_1_dir, replicate_2_dir
    else: 
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=True)

        loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)
    
    ##################################################################################################################

    bool_reprod_report = single_point_repr(
        loci1, loci2, ovr_threshold=to, window_bp=w, posterior=True, reproducibility_threshold=tr)

    ##################################################################################################################

    reproduced_loci1 = keep_reproducible_annotations(loci1, bool_reprod_report)
    write_MAPloci_in_BED(reproduced_loci1, savedir)
    # reproduced_loci1.to_csv(savedir + "/confident_posteriors.bed", sep='\t', header=True, index=False)

    ##################################################################################################################

    bool_reprod_report = pd.concat(
        [loci1["chr"], loci1["start"], loci1["end"], pd.Series(bool_reprod_report)], axis=1)
    bool_reprod_report.columns = ["chr", "start", "end", "is_repr"]
        
    # bool_reprod_report.to_csv(savedir+"/boolean_reproducibility_report_POSTERIOR.csv")

    general_rep_score = len(bool_reprod_report.loc[bool_reprod_report["is_repr"]==True]) / len(bool_reprod_report)
    perlabel_rec = {}
    lab_rep = perlabel_is_reproduced(bool_reprod_report, loci1)
    for k, v in lab_rep.items():
        perlabel_rec[k] = v[0]/ (v[0]+v[1])

    with open(savedir+"/ratio_robust.txt", "w") as scorefile:
        scorefile.write("general ratio robust = {}\n".format(general_rep_score))
        for k, v in perlabel_rec.items():
            scorefile.write("{} ratio robust = {}\n".format(k, str(v)))

    ##################################################################################################################
    rvalues = is_repr_posterior(
            loci1, loci2, ovr_threshold=to, window_bp=w, matching="static",
            always_include_best_match=True, return_r=True)
    
    rvalues.to_csv(savedir+"/r_values.bed", sep='\t', header=True, index=False)
    convert_to_GenomeBrowser_viewable_BED(savedir+"/r_values.bed")

    avg_r = np.mean(np.array(rvalues["r_value"]))
    perlabel_r = {}

    labels = rvalues.MAP.unique()

    for l in range(len(labels)):
        r_l = rvalues.loc[rvalues["MAP"] == labels[l], "r_value"]
        perlabel_r[labels[l]] = np.mean(np.array(r_l))
    
    with open(savedir+"/r_values_report.txt", "w") as scorefile:
        scorefile.write("GW average r_value = {}\n".format(avg_r))
        for k, v in perlabel_r.items():
            scorefile.write("{} average r_value = {}\n".format(k, str(v)))

    ##################################################################################################################


    MAP_NMI =  NMI_from_matrix(joint_overlap_prob(loci1, loci2, w=0, symmetric=True))
    POST_NMI = NMI_from_matrix(joint_prob_MAP_with_posterior(loci1, loci2, n_bins=200, conditional=False, stratified=True), posterior=True)

    with open(savedir+"/NMI.txt", "w") as scorefile:
        scorefile.write("NMI with MAP = {}\n".format(MAP_NMI))
        scorefile.write("NMI with binned posterior of R1 (n_bins=200) = {}\n".format(POST_NMI))


def post_clustering_keep_k_states(replicate_1_dir, replicate_2_dir, savedir, k, w, locis=False, write_csv=True, r_val=True):
    if locis:
        loci_1, loci_2 = replicate_1_dir, replicate_2_dir

    else:
        loci_1, loci_2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=False)

        loci_1, loci_2 = process_data(loci_1, loci_2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    joint = joint_overlap_prob(loci_1, loci_2, w=0, symmetric=False)

    ###################################################################################
    while loci_1.shape[1]-3 > 1:
        joint = joint_overlap_prob(loci_1, loci_2, w=0, symmetric=True)
        
        loci_1, loci_2 = merge_clusters(joint, loci_1, loci_2, r1=True)
        loci_1, loci_2 = merge_clusters(joint, loci_1, loci_2, r1=False)

        if loci_1.shape[1]-3 == k:
            if write_csv:
                loci_1.to_csv(savedir + f"/{str(k)}_states_post_clustered_posterior.csv")
            loci_1.to_csv(savedir + f"/{str(k)}_states_post_clustered_posterior.bed", sep='\t', header=True, index=False)

            MAP = loci_1.iloc[:,3:].idxmax(axis=1)
            coordMAP = pd.concat([loci_1.iloc[:, :3], MAP], axis=1)
            coordMAP.columns = ["chr", "start", "end", "MAP"]
            denseMAP = condense_segments(coordMAP)
            if write_csv:
                denseMAP.to_csv(savedir + f"/{str(k)}_states_confident_segments_dense.csv")
            denseMAP.to_csv(savedir + f"/{str(k)}_states_confident_segments_dense.bed", sep='\t', header=True, index=False)

            to = 0.75
            rvalues = is_repr_posterior(
                loci_1, loci_2, ovr_threshold=to, window_bp=w, matching="static",
                always_include_best_match=True, return_r=True)

            rvalues.to_csv(savedir+f"/r_values_{k}_states.bed", sep='\t', header=True, index=False)
            convert_to_GenomeBrowser_viewable_BED(savedir+f"/r_values_{k}_states.bed")
            return


def quick_report(replicate_1_dir, replicate_2_dir, savedir, locis=False, w=1000, to=0.75, tr=0.9):
    if locis:
        loci1, loci2 = replicate_1_dir, replicate_2_dir
    else:
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=True)

        loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)
    
    ct_confus(loci1.copy(), loci2.copy(), savedir, w=w)
    ct_granul(loci1.copy(), loci2.copy(), savedir)
    ct_boundar(loci1.copy(), loci2.copy(), savedir, match_definition="BM", max_distance=50)
    ##################################################################################################################

    bool_reprod_report = single_point_repr(
        loci1, loci2, ovr_threshold=to, window_bp=w, posterior=True, reproducibility_threshold=tr)

    ##################################################################################################################

    bool_reprod_report = pd.concat(
        [loci1["chr"], loci1["start"], loci1["end"], pd.Series(bool_reprod_report)], axis=1)
    bool_reprod_report.columns = ["chr", "start", "end", "is_repr"]

    general_rep_score = len(bool_reprod_report.loc[bool_reprod_report["is_repr"]==True]) / len(bool_reprod_report)
    perlabel_rec = {}
    lab_rep = perlabel_is_reproduced(bool_reprod_report, loci1)
    for k, v in lab_rep.items():
        perlabel_rec[k] = v[0]/ (v[0]+v[1])

    with open(savedir+"/ratio_robust.txt", "w") as scorefile:
        scorefile.write("general ratio robust = {}\n".format(general_rep_score))
        for k, v in perlabel_rec.items():
            scorefile.write("{} ratio robust = {}\n".format(k, str(v)))

    ##################################################################################################################
    rvalues = is_repr_posterior(
            loci1, loci2, ovr_threshold=to, window_bp=w, matching="static",
            always_include_best_match=True, return_r=True)

    avg_r = np.mean(np.array(rvalues["r_value"]))
    perlabel_r = {}

    labels = rvalues.MAP.unique()

    for l in range(len(labels)):
        r_l = rvalues.loc[rvalues["MAP"] == labels[l], "r_value"]
        perlabel_r[labels[l]] = np.mean(np.array(r_l))
    
    with open(savedir+"/r_values_report.txt", "w") as scorefile:
        scorefile.write("GW average r_value = {}\n".format(avg_r))
        for k, v in perlabel_r.items():
            scorefile.write("{} average r_value = {}\n".format(k, str(v)))

    ##################################################################################################################

    MAP_NMI =  NMI_from_matrix(joint_overlap_prob(loci1, loci2, w=0, symmetric=True))
    POST_NMI = NMI_from_matrix(joint_prob_MAP_with_posterior(loci1, loci2, n_bins=200, conditional=False, stratified=True), posterior=True)

    with open(savedir+"/NMI.txt", "w") as scorefile:
        scorefile.write("NMI with MAP = {}\n".format(MAP_NMI))
        scorefile.write("NMI with binned posterior of R1 (n_bins=200) = {}\n".format(POST_NMI))


def overlap_vs_segment_length(replicate_1_dir, replicate_2_dir, savedir, locis=True, custom_bin=True):
    if locis:
        loci1, loci2 = replicate_1_dir, replicate_2_dir
    else:
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=False)

        loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)
    ##################################################################################################################
    if os.path.exists(f"{savedir}/overlap_vs_segmentlength/") == False:
        os.mkdir(f"{savedir}/overlap_vs_segmentlength/")
    savedir = f"{savedir}/overlap_vs_segmentlength/"

    MAP1 = loci1.iloc[:,3:].idxmax(axis=1)
    MAP1 = pd.concat([loci1.chr, loci1.start, loci1.end, MAP1], axis=1)
    MAP1.columns = ["chr", "start", "end", "MAP"]

    MAP2 = loci2.iloc[:,3:].idxmax(axis=1)
    MAP2 = pd.concat([loci2.chr, loci2.start, loci2.end, MAP2], axis=1)
    MAP2.columns = ["chr", "start", "end", "MAP"]
    resolution = MAP1["end"][0] - MAP1["start"][0]
    #############################################ESTABLISH CORRESPONDENCE#############################################
    confmat = IoU_overlap(loci1, loci2, w=0, symmetric=True, soft=False)
    
    # define matches
    per_label_matches = {}
    for k in list(loci1.columns[3:]):
        sorted_k_vector = confmat.loc[k, :].sort_values(ascending=False)
        per_label_matches[k] = sorted_k_vector.index[0]
    ##################################################################################################################
    interpretation_terms = ["Prom", "Prom_fla", "Enha", "Enha_low", "Biva", "Tran", "Cons", "Facu", "K9K3", "Quie", "Unkn"]

    MAP1 = MAP1.to_numpy()
    MAP2 = MAP2.to_numpy()

    all_segs = {term: [] for term in interpretation_terms}
    # all_segs = {term: [] for term in loci1.columns[3:]}
    ##################################################################################################################
    i = 0
    current_map = MAP1[i, 3]
    current_start = MAP1[i, 1]
    current_seg = []

    is_matched = bool(MAP2[i, 3] == per_label_matches[MAP1[i, 3]])
    #middle of the segment
    if is_matched:
        current_seg.append(
            [MAP1[i, 1] - current_start, 
            1])
    else:
        current_seg.append(
            [MAP1[i, 1] - current_start, 
            0])

    for i in range(1, len(MAP1)):
        if MAP1[i, 0] == MAP1[i-1, 0] and MAP1[i, 3] == current_map:
            is_matched = bool(MAP2[i, 3] == per_label_matches[MAP1[i, 3]])
            #middle of the segment
            if is_matched:
                current_seg.append(
                    [MAP1[i, 1] - current_start, 
                    1])
            else:
                current_seg.append(
                    [MAP1[i, 1] - current_start, 
                    0])
            
        else:
            #last_segment_ends
            seg_length = MAP1[i-1, 2] - current_start
            current_seg = np.reshape(np.array(current_seg).astype(float), (-1, 2))

            current_seg[:, 0] = (current_seg[:, 0] + (resolution/2)) / seg_length

            if len(current_seg) < 1:
                print(current_seg)

            translate_to_term = max([x for x in interpretation_terms if x in current_map], key=len)
            all_segs[translate_to_term].append(current_seg)

            # all_segs[current_map].append(current_seg)

            #new_segment
            current_map = MAP1[i, 3]
            current_start = MAP1[i, 1]
            current_seg = []

            is_matched = bool(MAP2[i, 3] == per_label_matches[MAP1[i, 3]])
            #middle of the segment
            if is_matched:
                current_seg.append(
                    [MAP1[i, 1] - current_start, 
                    1])
            else:
                current_seg.append(
                    [MAP1[i, 1] - current_start, 
                    0])

    ##################################################################################################################
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
            ax0.set_ylabel('naive overlap')
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
            ax1.set_ylabel('naive overlap')
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
            ax2.set_ylabel('naive overlap')
        except:
            pass

    # Show the plot
    # Create a separate subplot for the legend at the top
    ax_legend = plt.subplot(gs[0, :])
    ax_legend.axis('off')  # Hide the axes

    # Show the legend in this subplot
    fig.legend(lines, labels, loc='center', ncol=len(labels), bbox_to_anchor=(0.5, 0.5), bbox_transform=ax_legend.transAxes)
    plt.tight_layout()
    plt.savefig(f"{savedir}/naive_overlap_v_segment_length.pdf", format='pdf')
    plt.savefig(f"{savedir}/naive_overlap_v_segment_length.svg", format='svg')


if __name__=="__main__":  
    
    test_new_functions(
        replicate_1_dir="tests/cedar_runs/chmm/GM12878_R1/", 
        replicate_2_dir="tests/cedar_runs/chmm/GM12878_R2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        savedir="tests/cedar_runs/chmm/GM12878_R1/")
    GET_ALL(
        replicate_1_dir="tests/cedar_runs/segway/GM12878_R1/", 
        replicate_2_dir="tests/cedar_runs/segway/GM12878_R2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv", 
        savedir="tests/cedar_runs/segway/GM12878_R1/", contour=False)
    
    GET_ALL(
        replicate_1_dir="tests/cedar_runs/chmm/GM12878_R1/", 
        replicate_2_dir="tests/cedar_runs/chmm/GM12878_R2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv", 
        savedir="tests/cedar_runs/chmm/GM12878_R1/", contour=False)

    exit()

    # test_new_functions(
    #     replicate_1_dir="tests/cedar_runs/chmm/GM12878_R1/", 
    #     replicate_2_dir="tests/cedar_runs/chmm/GM12878_R2/", 
    #     genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
    #     savedir="tests/cedar_runs/chmm/GM12878_R1/")
    
    # test_new_functions(
    #     replicate_1_dir="tests/cedar_runs/segway/GM12878_R1/", 
    #     replicate_2_dir="tests/cedar_runs/segway/GM12878_R2/", 
    #     genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
    #     savedir="tests/cedar_runs/segway/GM12878_R1/")

    # test_new_functions(
    #     replicate_1_dir="tests/cedar_runs/segway/GM12878_R1/", 
    #     replicate_2_dir="tests/cedar_runs/segway/GM12878_R2/", 
    #     genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
    #     savedir="tests/cedar_runs/segway/GM12878_R1/")
    # exit()
    # test_new_functions(
    #     replicate_1_dir="tests/cedar_runs/segway/GM12878_R1/", 
    #     replicate_2_dir="tests/cedar_runs/segway/GM12878_R2/", 
    #     genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
    #     savedir="tests/cedar_runs/segway/GM12878_R1/")

    # exit()
    # post_clustering(
    #     replicate_1_dir="tests/cedar_runs/chmm/GM12878_R1/", 
    #     replicate_2_dir="tests/cedar_runs/chmm/GM12878_R2/", 
    #     savedir="tests/cedar_runs/chmm/GM12878_R1/")

    # post_clustering(
    #     replicate_1_dir="tests/cedar_runs/segway/GM12878_R1/", 
    #     replicate_2_dir="tests/cedar_runs/segway/GM12878_R2/", 
    #     savedir="tests/cedar_runs/segway/GM12878_R1/")

    # exit()

    # GET_ALL(
    #     replicate_1_dir="tests/cedar_runs/segway_concat/GM12878_concat_rep1/", 
    #     replicate_2_dir="tests/cedar_runs/segway_concat/GM12878_concat_rep2/", 
    #     genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
    #     rnaseq="biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv", 
    #     savedir="tests/cedar_runs/segway_concat/GM12878_concat_rep1/", contour=False)

    # GET_ALL(
    #     replicate_1_dir="tests/cedar_runs/segway_concat/K562_concat_rep1/", 
    #     replicate_2_dir="tests/cedar_runs/segway_concat/K562_concat_rep2/", 
    #     genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
    #     rnaseq="biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv", 
    #     savedir="tests/cedar_runs/segway_concat/K562_concat_rep1/", contour=False)

    
    exit()

    GET_ALL(
        replicate_1_dir="tests/cedar_runs/segway/GM12878_R1/", 
        replicate_2_dir="tests/cedar_runs/segway/GM12878_R2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv", 
        savedir="tests/cedar_runs/segway/GM12878_R1/", contour=True)

    exit()
    
    GET_ALL(
        replicate_1_dir="tests/cedar_runs/chmm/GM12878_R2/", 
        replicate_2_dir="tests/cedar_runs/chmm/GM12878_R1/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv", 
        savedir="tests/cedar_runs/chmm/GM12878_R2/")

    GET_ALL(
        replicate_1_dir="tests/cedar_runs/segway/GM12878_R2/", 
        replicate_2_dir="tests/cedar_runs/segway/GM12878_R1/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv", 
        savedir="tests/cedar_runs/segway/GM12878_R2/")
