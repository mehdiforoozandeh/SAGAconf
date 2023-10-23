from src.run import *
from src.indepth import *
from src.granul import *
from src.bio_valid import *
from src.overall import *

from matplotlib.colors import LinearSegmentedColormap
import ast, pybedtools
import matplotlib.gridspec as gridspec
from scipy.interpolate import UnivariateSpline
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm

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

#######################################################################################################################

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

def plot_raw_and_iou(replicate_1_dir, replicate_2_dir):
    loci_1, loci_2 = load_data(
        replicate_1_dir+"/parsed_posterior.csv",
        replicate_2_dir+"/parsed_posterior.csv",
        subset=True, logit_transform=False)

    loci_1, loci_2 = process_data(loci_1, loci_2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False, custom_order=True) 
    savedir = replicate_1_dir
    boundaries = [x for x in list(np.linspace(0, 1, 16))]  # custom boundaries
    hex_colors = sns.light_palette('navy', n_colors=len(boundaries) * 2 + 2, as_cmap=False).as_hex()
    hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]
    colors=list(zip(boundaries, hex_colors))
    custom_color_map = LinearSegmentedColormap.from_list(
        name='custom_navy',
        colors=colors)

    # ####################################################################################
    confmat = joint_overlap_prob(loci_1, loci_2, w=0, symmetric=True)
    confmat = confmat + 1e-6

    vmin, vmax = 1e-6, confmat.max().max()
    # Create an instance of the LogNorm class
    norm = LogNorm(vmin=vmin, vmax=vmax)
    # sns.set(rc={'figure.figsize':(15,15)})

    p = sns.heatmap(
    confmat.astype(float), annot=True, fmt=".0e",
        linewidths=0.01, cbar=True, annot_kws={"size": 6, "rotation":-45},
        vmin=vmin, vmax=vmax, cmap=custom_color_map, norm=norm)

    
    p.tick_params(axis='x', rotation=90, labelsize=10)
    p.tick_params(axis='y', rotation=0, labelsize=10)

    plt.xlabel('Replicate 2 Labels')
    plt.ylabel("Replicate 1 Labels")
    plt.tight_layout()
    plt.savefig('{}/raw.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/raw.svg'.format(savedir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    # ####################################################################################
    confmat = IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False, overlap_coeff=False)
    confmat = confmat + 1e-6

    vmin, vmax = 1e-6, confmat.max().max()
    # Create an instance of the LogNorm class
    norm = LogNorm(vmin=vmin, vmax=vmax)
    
    p = sns.heatmap(
    confmat.astype(float), annot=True, fmt=".2f",
        linewidths=0.01, cbar=True, annot_kws={"size": 6},
        vmin=vmin, vmax=vmax, cmap=custom_color_map, norm=norm)

    # sns.set(rc={'figure.figsize':(15,15)})
    p.tick_params(axis='x', rotation=90, labelsize=10)
    p.tick_params(axis='y', rotation=0, labelsize=10)
    # Iterate over the rows of the heatmap data
    for i, row in enumerate(confmat.values):
        # Find the column index of the maximum value in this row
        j = row.argmax()
        
        # Add a green border to the cell at location (i,j)
        p.add_patch(Rectangle((j, i), 1, 1, edgecolor='red', fill=False, lw=3))

    plt.title('IoU Overlap')
    plt.xlabel('Replicate 2 Labels')
    plt.ylabel("Replicate 1 Labels")
    plt.tight_layout()
    plt.savefig('{}/iou_logscale.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/iou_logscale.svg'.format(savedir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    # ####################################################################################
    boundaries = [x**2 for x in list(np.linspace(0, 1, 16))] + [1] # custom boundaries
    hex_colors = sns.light_palette('navy', n_colors=len(boundaries) * 2 + 2, as_cmap=False).as_hex()
    hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]
    colors=list(zip(boundaries, hex_colors))
    custom_color_map = LinearSegmentedColormap.from_list(
        name='custom_navy',
        colors=colors)

    # ####################################################################################

    confmat = IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False, overlap_coeff=False)
    # confmat = confmat + 1e-6

    # vmin, vmax = 1e-6, confmat.max().max()
    # Create an instance of the LogNorm class
    # norm = LogNorm(vmin=vmin, vmax=vmax)
    
    p = sns.heatmap(
    confmat.astype(float), annot=True, fmt=".2f",
        linewidths=0.01, cbar=True, annot_kws={"size": 6},
        vmin=vmin, vmax=vmax, cmap=custom_color_map)#, norm=norm)

    # sns.set(rc={'figure.figsize':(15,15)})
    p.tick_params(axis='x', rotation=90, labelsize=10)
    p.tick_params(axis='y', rotation=0, labelsize=10)
    # Iterate over the rows of the heatmap data
    for i, row in enumerate(confmat.values):
        # Find the column index of the maximum value in this row
        j = row.argmax()
        
        # Add a green border to the cell at location (i,j)
        p.add_patch(Rectangle((j, i), 1, 1, edgecolor='red', fill=False, lw=3))

    plt.title('IoU Overlap')
    plt.xlabel('Replicate 2 Labels')
    plt.ylabel("Replicate 1 Labels")
    plt.tight_layout()
    plt.savefig('{}/iou.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/iou.svg'.format(savedir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

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

def distance_vs_overlap_2(loci_1, loci_2, savedir, match_definition="BM"):
    indicator_file = savedir+"/Dist_vs_Corresp_2"
    if os.path.exists(indicator_file):
        return
        
    if os.path.exists(savedir+"/Dist_vs_Corresp_2") == False:
        os.mkdir(savedir+"/Dist_vs_Corresp_2")
    savedir = savedir+"/Dist_vs_Corresp_2"

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
    #========================================================================================#
    
    """
    for each position i
        for each distance d
            what fraction of i+d is the corresponding label?
    """

    distance_to_corresp = {}
    for d in range(max_distance):
        distance_to_corresp[d] = [0, 0] #[num_correspond, num_all]

    for i in range(len(loci_1)):
        for d in range(max_distance):
            if i+d < len(loci_1):
                distance_to_corresp[d][1] +=1

                if MAP2[i+d] == per_label_matches[MAP1[i]]:
                    distance_to_corresp[d][0] +=1

            elif i-d >= 0:
                distance_to_corresp[d][1] +=1

                if MAP2[i-d] == per_label_matches[MAP1[i]]:
                    distance_to_corresp[d][0] +=1
    
    for d in distance_to_corresp.keys():
        distance_to_corresp[d] =  float(distance_to_corresp[d][0] / distance_to_corresp[d][1])

    xaxis = [x*resolution for x in distance_to_corresp.keys()]
    plt.plot(xaxis, distance_to_corresp.values(), color="black")
    plt.fill_between(xaxis, distance_to_corresp.values(), color="black", alpha=0.4)
    plt.xlabel("Distance (bp)")
    plt.ylabel("Probability of overlap with corresponding label")
    plt.title("Correspondence vs. Distance - Overall")

    plt.savefig(savedir+"/dist_vs_corresp_{}.pdf".format("overall"), format='pdf')
    plt.savefig(savedir+"/dist_vs_corresp_{}.svg".format("overall"), format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    
    #========================================================================================#
    
    distance_to_corresp = {}
    for k in loci_1.columns[3:]:
        distance_to_corresp[k] = {}
        for d in range(max_distance):
            distance_to_corresp[k][d] = [0, coverage1[k]*len(loci_1)] #[num_correspond, num_all]

    for i in range(len(loci_1)):
        for d in range(max_distance):
            if i+d < len(loci_1):
                distance_to_corresp[MAP1[i]][d][1] +=1
                if MAP2[i+d] == per_label_matches[MAP1[i]]:
                    distance_to_corresp[MAP1[i]][d][0] += 1

            if i-d >= 0:
                distance_to_corresp[MAP1[i]][d][1] +=1
                if MAP2[i-d] == per_label_matches[MAP1[i]]:
                    distance_to_corresp[MAP1[i]][d][0] += 1
    
    for k in distance_to_corresp.keys():
        for d in distance_to_corresp[k].keys():
            distance_to_corresp[k][d] =  float(distance_to_corresp[k][d][0] / distance_to_corresp[k][d][1])

    distance_to_corresp = pd.DataFrame(distance_to_corresp)
    distance_to_corresp.index = [x*resolution for x in distance_to_corresp.index]

    print(distance_to_corresp)

    for k in distance_to_corresp.columns:
        plt.plot(
            distance_to_corresp.index, distance_to_corresp.loc[:, k], color="black")

        plt.fill_between(
            distance_to_corresp.index, distance_to_corresp.loc[:, k], color="black", alpha=0.4)

        plt.axhline(y=coverage2[per_label_matches[k]], color='red', linestyle='--', linewidth=2)

        plt.xlabel("Distance (bp)")
        plt.ylabel("Probability of overlap with corresponding label")
        plt.title("Correspondence vs. Distance - {}".format(k))

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
            k = loci_1.columns[3:][label_being_plotted]
            axs[i,j].plot(
                distance_to_corresp.index, distance_to_corresp.loc[:, k], color="black")

            axs[i,j].fill_between(
                distance_to_corresp.index, distance_to_corresp.loc[:, k], color="black", alpha=0.4)
            
            axs[i,j].axhline(y=coverage2[per_label_matches[k]], color='red', linestyle='--', linewidth=1.5)

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

def label_granularity(loci_1, loci_2, savedir):
    indicator_file = savedir+"/granul/"
    if os.path.exists(indicator_file):
        return

    savedir = savedir+"/granul/"
    if os.path.exists(savedir)==False:
        os.mkdir(savedir)
    
    for c in list(loci_1.columns[3:]):
        cr, ar, ovr_rec = granularity_vs_agreement_nonsymmetric(loci_1.copy(), loci_2.copy(), k=c)

        perfect_agr = [0] + [1 for i in range(len(ar) - 1)]
        realAUC = metrics.auc(cr, ar)
        perfectAUC = metrics.auc(cr, perfect_agr)
        p_to_r_auc = float(realAUC/perfectAUC)

        plt.plot(cr, ar, c="yellowgreen")
        plt.title(
            c + str(" | Real/Perfect AUC = {:.2f}".format(p_to_r_auc)), 
            fontsize=13)

        plt.fill_between(cr, ar, color="yellowgreen", alpha=0.4)
        plt.fill_between(cr, perfect_agr, ar, color="palevioletred", alpha=0.4)
        plt.xticks(np.arange(0, 1.1, step=0.2))
        plt.yticks(np.arange(0, 1.1, step=0.2))
    
        plt.xlabel('Coverage' , ha='center')
        plt.ylabel('Agreement', va='center', rotation='vertical')
        plt.tight_layout()
        plt.savefig(savedir+"/granularity_{}.pdf".format(c), format='pdf')
        plt.savefig(savedir+"/granularity_{}.svg".format(c), format='svg')

        sns.reset_orig
        plt.close("all")
        plt.style.use('default')
        plt.clf()

def label_boundary(loci_1, loci_2, savedir, match_definition="BM", max_distance=50):
    indicator_file = savedir+"/len_bound/"
    if os.path.exists(indicator_file):
        return

    savedir = savedir+"/len_bound/"
    if os.path.exists(savedir)==False:
        os.mkdir(savedir)

    """
    Here i'm gonna merge ECDF of length dist and boundary thing.
    """

    resolution = int(loci_1.iloc[0, 2] - loci_1.iloc[0, 1])

    max_distance = int(5000/resolution)
    num_labels = loci_1.shape[1]-3
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
    
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

    num_labels = len(boundary_distances.keys())

    for k in loci_1.columns[3:]:
        matched_hist = np.histogram(boundary_distances[k], bins=max_distance, range=(0, max_distance)) #[0]is the bin size [1] is the bin edge

        cdf = []
        for jb in range(len(matched_hist[0])):
            if len(cdf) == 0:
                cdf.append(matched_hist[0][jb])
            else:
                cdf.append((cdf[-1] + matched_hist[0][jb]))

        cdf = np.array(cdf) / (np.sum(matched_hist[0]) + count_unmatched_boundaries[k])
        
        plt.plot(list(matched_hist[1][:-1]*resolution), list(cdf), color="black", label="distance from boundary")
        
        # lendist = get_ecdf_for_label(loci_1, k, max=(max_distance+1)*resolution)
        # if len(lendist)>0:
        #     sns.ecdfplot(lendist, color="red")

        plt.xticks(np.arange(0, (max_distance+1)*resolution, step=5*resolution))
        plt.tick_params(axis='both', labelsize=12)
        plt.tick_params(axis='x', rotation=90)
        plt.yticks(np.arange(0, 1.1, step=0.2))
        plt.xlabel('bp', fontsize=12)
        plt.ylabel('Ratio Overlap',fontsize=12)
        
        plt.title(k, fontsize=12)
            
        plt.tight_layout()
        plt.savefig(savedir+"/len_bound_{}.pdf".format(k), format='pdf')
        plt.savefig(savedir+"/len_bound_{}.svg".format(k), format='svg')

        sns.reset_orig
        plt.close("all")
        plt.style.use('default')
        plt.clf()

def label_merging_progression(loci_1, loci_2, savedir):
    indicator_file = savedir+"/prog/"
    if os.path.exists(indicator_file):
        return

    savedir = savedir+"/prog/"
    if os.path.exists(savedir)==False:
        os.mkdir(savedir)

    for c in list(loci_1.columns[3:]):
        cr, ar, ovr_rec = granularity_vs_agreement_nonsymmetric(loci_1.copy(), loci_2.copy(), k=c)
        plot_progression(ar, cr, ovr_rec, c, savedir)

        sns.reset_orig
        plt.close("all")
        plt.style.use('default')

#######################################################################################################################

def get_all_labels(replicate_1_dir, replicate_2_dir, savedir, locis=False):
    """
    -single granularity
    -single prog
    -single boundar
    """
    if locis:
        loci1, loci2 = replicate_1_dir, replicate_2_dir
    else:
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True)

        loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    label_granularity(loci1, loci2, savedir)
    label_merging_progression(loci1, loci2, savedir)
    label_boundary(loci1, loci2, savedir, match_definition="BM", max_distance=50)

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

def gather_labels(original_ct_dir, savedir, contour=True):
    
    print("loading mnemonics...")
    if os.path.exists(
        "/".join(original_ct_dir.split("/")[:-1])+"/mnemonics_rep1.txt"):

        print("reading concat mnemonics")
        loci_1_mnemon = read_mnemonics("/".join(original_ct_dir.split("/")[:-1])+"/mnemonics_rep1.txt")
    else:
        print("reading mnemonics")
        loci_1_mnemon = read_mnemonics("/".join(original_ct_dir.split("/")[:-1])+"/mnemonics.txt")

    mnemon1_dict = {}
    for i in loci_1_mnemon:
        mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]

    print(mnemon1_dict)

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
    if os.path.exists(savedir+"/labels/")==False:
        os.mkdir(savedir+"/labels/")

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

    for k in mnemon1_dict.values():
        if os.path.exists("{}/labels/{}".format(savedir, k)) == False:
            os.mkdir("{}/labels/{}".format(savedir, k))
        
        for f in os.listdir(savedir+"/len_bound/"):
            if "len_bound_"+k in f:
                os.system(
                    "cp {} {}".format(savedir+"/len_bound/"+f, savedir+"/labels/"+k)
                    )
        
        for f in os.listdir(savedir+"/calib/"):
            if "caliberation_"+k in f:
                os.system(
                    "cp {} {}".format(savedir+"/calib/"+f, savedir+"/labels/"+k)
                    )

        for f in os.listdir(savedir+"/granul/"):
            if "granularity_"+k in f:
                os.system(
                    "cp {} {}".format(savedir+"/granul/"+f, savedir+"/labels/"+k)
                    )

        for f in os.listdir(savedir+"/prog/"):
            if "prog_"+k in f:
                os.system(
                    "cp {} {}".format(savedir+"/prog/"+f, savedir+"/labels/"+k)
                    )
        if contour:
            for f in os.listdir(savedir+"/contours_OvrWind/"):
                if "reprod_lines_"+k in f:
                    os.system(
                        "cp {} {}".format(savedir+"/contours_OvrWind/"+f, savedir+"/labels/"+k)
                        )
                
                if "reprod_contour_"+k in f:
                    os.system(
                        "cp {} {}".format(savedir+"/contours_OvrWind/"+f, savedir+"/labels/"+k)
                        )
                    
            for f in os.listdir(savedir+"/contours_ReprThresWind/"):
                if "reprod_lines_"+k in f:
                    os.system(
                        "cp {} {}".format(savedir+"/contours_ReprThresWind/"+f, savedir+"/labels/"+k)
                        )
                    
                if "reprod_contour_"+k in f:
                    os.system(
                        "cp {} {}".format(savedir+"/contours_ReprThresWind/"+f, savedir+"/labels/"+k)
                        )
                
def get_all_bioval(replicate_1_dir, replicate_2_dir, savedir, genecode_dir, rnaseq=None, locis=False):
    indicator_file = savedir+"/tss_enr"
    if os.path.exists(indicator_file):
        return

    if locis:
        loci1, loci2 = replicate_1_dir, replicate_2_dir
    else:
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=False)

        loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    gene_coords = load_gene_coords(genecode_dir)
    if rnaseq != None:
        trans_data = load_transcription_data(rnaseq, gene_coords)

        trans_data = trans_data.drop(trans_data[trans_data.TPM==0].index).reset_index(drop=True)

        posterior_transcription_enrichment(loci1, trans_data, savedir+"/trans_post_enr")
        posterior_transcription_correlation(loci1, trans_data, savedir=savedir+"/trans_post_correl")

        posterior_transcription_enrichment(loci1, trans_data, TSS=True, savedir=savedir+"/trans_post_enr_aroundTSS")

    overal_TSS_enrichment(loci1, savedir+"/tss_enr")

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

def get_contour(replicate_1_dir, replicate_2_dir, savedir, locis=False):
    indicator_file = savedir+"/contours_OvrWind/"
    if os.path.exists(indicator_file):
        return

    if locis:
        loci1, loci2 = replicate_1_dir, replicate_2_dir
    else:
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=False)

        loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    print("getting contours 1 : overlapT-window-repr")
    OvrWind_contour(
        loci1, loci2, savedir, w_range=[0, 2600, 500], t_range=[0, 11, 2], posterior=True, repr_threshold=0.75)
    
    # print("getting contours 2 : reprT-window-repr")
    # ReprThresWind_contour(
    #     loci1, loci2, savedir, w_range=[0, 4000, 600], t_range=[50, 105, 5], posterior=True, matching="static", static_thres=0.75)
    
    # print("getting contours 3 : overlapT-window-deltaNMI")
    # OvrWind_delta_NMI_contour(
    #     loci1, loci2, savedir, w_range=[0, 3000, 500], t_range=[0, 11, 2], posterior=True, repr_threshold=0.75)
    
    # print("getting contours 4 : reprT-window-deltaNMI")
    # ReprThresWind_delta_NMI_contour(
    #     loci1, loci2, savedir, w_range=[0, 3000, 500], t_range=[50, 100, 15], posterior=True, matching="static")

def after_SAGAconf_metrics(replicate_1_dir, replicate_2_dir, genecode_dir, savedir, w=1000, rnaseq=None, intersect_r1r2=False, locis=False, to=0.75, tr=0.9):
    indicator_file = savedir+"/after_SAGAconf/NMI.txt"
    if os.path.exists(indicator_file):
        return

    if locis:
        loci_1, loci_2 = replicate_1_dir, replicate_2_dir
    else:
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=True)
        loci_1, loci_2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    savedir = savedir+"/after_SAGAconf"
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)
    
    isrep1 = is_repr_posterior(
        loci_1, loci_2, ovr_threshold=to, window_bp=w, reproducibility_threshold=tr, matching="static")
    
    isrep2 = is_repr_posterior(
            loci_2, loci_1, ovr_threshold=to, window_bp=w, reproducibility_threshold=tr, matching="static")
    
    is_rep_intersec = []
    conf = np.zeros((2,2))
    for i in range(len(isrep1)):
        if isrep1[i]==isrep2[i]==1:
            is_rep_intersec.append(True)
            conf[1,1] += 1

        else:
            is_rep_intersec.append(False)
            if isrep1[i]==1 and isrep2[i]==0:
                conf[1,0] += 1

            elif isrep1[i]==0 and isrep2[i]==1:
                conf[0,1] += 1

            elif isrep1[i]==0 and isrep2[i]==0:
                conf[0,0] += 1

    try:
        with open(savedir+"/sym_r1r2_repr.txt", "w") as scorefile:
            scorefile.write("reproduced in both = {}\n".format(float(conf[1,1] / len(isrep1))))
            scorefile.write("not-reproduced in both = {}\n".format(float(conf[0,0] / len(isrep1))))
            scorefile.write("fraction similar = {}\n".format(float((conf[1,1] + conf[0,0]) / len(isrep1))))

    except:
        pass

    if intersect_r1r2:
        # at this point we can get the intersection of isrep1 and isrep2 as the positions that are considered 
        calibrated_loci1 = keep_reproducible_annotations(loci_1, is_rep_intersec) 
        calibrated_loci2 = keep_reproducible_annotations(loci_2, is_rep_intersec)

    else:
        #do not intersect and just use isrep1
        calibrated_loci1 = keep_reproducible_annotations(loci_1, isrep1) 
        calibrated_loci2 = keep_reproducible_annotations(loci_2, isrep1)

    # ========================================================================================================= #
    try:
        MAP_NMI =  NMI_from_matrix(
            joint_overlap_prob(calibrated_loci1, calibrated_loci2, w=0, symmetric=True))
        POST_NMI = NMI_from_matrix(
            joint_prob_with_binned_posterior(
                calibrated_loci1, calibrated_loci2, n_bins=200, conditional=False, stratified=True), 
                posterior=True)

        with open(savedir+"/NMI.txt", "w") as scorefile:
            scorefile.write("NMI with MAP = {}\n".format(MAP_NMI))
            scorefile.write("NMI with binned posterior of R1 (n_bins=200) = {}\n".format(POST_NMI))
    except:
        print("FAILED. EXCEPTION...")
    
    try:
        ct_granul(calibrated_loci1, calibrated_loci2, savedir)
    except:
        print("FAILED. EXCEPTION...")

    try:
        conditional = overlap_matrix(calibrated_loci1, calibrated_loci2, type="conditional")
        with open("{}/raw_conditional_overlap_ratio.txt".format(savedir), "w") as cf:
            for i in conditional.index:
                cf.write("{} : {}\n".format(i, np.max(np.array(conditional.loc[i, :]))))
    except:
        print("FAILED. EXCEPTION...")

    try:
        gene_coords = load_gene_coords(genecode_dir)
        if rnaseq != None:
            trans_data = load_transcription_data(rnaseq, gene_coords)

            # trans_data = trans_data.drop(trans_data[trans_data.TPM==0].index).reset_index(drop=True)

            posterior_transcription_enrichment(calibrated_loci1, trans_data, savedir+"/trans_post_enr")
            posterior_transcription_correlation(calibrated_loci1, trans_data, savedir=savedir+"/trans_post_correl")

            posterior_transcription_enrichment(calibrated_loci1, trans_data, TSS=True, savedir=savedir+"/trans_post_enr_aroundTSS")

        overal_TSS_enrichment(loci1, savedir+"/tss_enr")
    except:
        print("FAILED. EXCEPTION...")

def before_after_saga(savedir):
    indicator_file = '{}/AUC_before_after_SAGAconf.pdf'.format(savedir)
    if os.path.exists(indicator_file):
        return

    auc_before = ast.literal_eval(open(savedir+"/AUC_mAUC.txt", "r").read())
    auc_after = ast.literal_eval(open(savedir+"/after_SAGAconf/AUC_mAUC.txt", "r").read())

    auc_df = []
    for k in auc_before.keys():
        auc_df.append([k, "Before_SAGAconf", auc_before[k]])
        auc_df.append([k, "After_SAGAconf", auc_after[k]])

    auc_df = pd.DataFrame(auc_df, columns = ["Chromatin State", "b/a", "AUC/mAUC"])
    sns.set_theme(style="whitegrid")
    sns.reset_orig
    plt.style.use('default')

    sns.set(rc={'figure.figsize':(10,8)})
    sns.set_palette([(0, 0, 0, 0.7), (0, 0.7, 0.3, 0.7)])

    sns.barplot(x="Chromatin State", y="AUC/mAUC", hue="b/a", data=auc_df, alpha=0.7)
    plt.axhline(y=0.5, linestyle='--', color='red', alpha=0.7)
    plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
    plt.ylabel("AUC/mAUC")
    plt.xticks(rotation=45)
    plt.yticks(np.arange(0, 1.1, step=0.1))
    plt.tight_layout()

    plt.savefig(savedir+"/AUC_before_after_SAGAconf.pdf", format="pdf")
    plt.savefig(savedir+"/AUC_before_after_SAGAconf.svg", format="svg")
    
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    sns.set_theme(style="whitegrid")
    sns.reset_orig
    plt.style.use('default')

    #########################################################################################################
    before_nmi = savedir + "/NMI.txt"
    after_nmi =  savedir + "/after_SAGAconf/NMI.txt"

    with open(before_nmi, "r") as nmif:
        nmilines = nmif.readlines()
        before = float(nmilines[0].split("=")[1][:7])
    
    with open(after_nmi, "r") as nmif:
        nmilines = nmif.readlines()
        after = float(nmilines[0].split("=")[1][:7])
    
    plt.bar(["Before SAGAconf", "After SAGAconf"], [before, after], color=[(0, 0, 0, 0.7), (0, 0.7, 0.3, 0.7)], alpha=0.7)
    plt.ylabel("NMI")
    plt.yticks(np.arange(0, 1.1, step=0.1))
    plt.tight_layout()

    plt.savefig(savedir+"/NMI_before_after_SAGAconf.pdf", format="pdf")
    plt.savefig(savedir+"/NMI_before_after_SAGAconf.svg", format="svg")
    
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()


    sns.set_theme(style="whitegrid")
    sns.reset_orig
    plt.style.use('default')

def post_clustering(replicate_1_dir, replicate_2_dir, savedir, locis=False, to=0.75, tr=0.9):
    """
    get the IoU
    make it symmetric
    get two separate dendrograms
        one on rows 
        one on columns

    while num_rows > 1:
        merge best pair of R1
        merge best pair of R2
        generate new dendrogram
    """
    indicator_file = '{}/MI_post_clustering_progress.txt'.format(savedir)
    if os.path.exists(indicator_file):
        return

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
    robust_rec = {}
    robust_rec_entropy = {}
    avgr_rec = {} # avgr vs num_labels
    avgr_entropy = {} # avgr vs enrtopy
    nmi_rec = {}
    NMI_rec_entropy = {}

    MI_rec_entropy = {}
    MI_rec = {}

    naive_rec = {}
    naive_ent = {}
    
    while loci_1.shape[1]-3 > 1 :
        bool_reprod_report, avg_r = single_point_repr(
            loci_1, loci_2, ovr_threshold=to, window_bp=1000, posterior=True, 
            reproducibility_threshold=tr, return_mean_r=True)

        bool_reprod_report = pd.concat(
            [loci_1["chr"], loci_1["start"], loci_1["end"], pd.Series(bool_reprod_report)], axis=1)
        bool_reprod_report.columns = ["chr", "start", "end", "is_repr"]

        general_rep_score = len(bool_reprod_report.loc[bool_reprod_report["is_repr"]==True]) / len(bool_reprod_report)
        perlabel_rec = {}

        lab_rep = perlabel_is_reproduced(bool_reprod_report, loci_1)

        for k, v in lab_rep.items():
            perlabel_rec[k] = v[0]/ (v[0]+v[1])

        joint = joint_overlap_prob(loci_1, loci_2, w=0, symmetric=True)
        NMI = NMI_from_matrix(joint)
        MI = NMI_from_matrix(joint, return_MI=True)
        loci1_entropy = calc_entropy(joint, rows=True)
        naive_overlap = overall_overlap_ratio(loci_1, loci_2, w=0)

        nmi_rec["{}".format(loci_1.shape[1]-3)] = NMI
        robust_rec["{}".format(loci_1.shape[1]-3)] = general_rep_score
        avgr_rec["{}".format(loci_1.shape[1]-3)] = avg_r
        MI_rec["{}".format(loci_1.shape[1]-3)] = MI
        naive_rec["{}".format(loci_1.shape[1]-3)] = naive_overlap

        MI_rec_entropy[loci1_entropy] = MI
        robust_rec_entropy[loci1_entropy] = general_rep_score
        NMI_rec_entropy[loci1_entropy] = NMI
        avgr_entropy[loci1_entropy] = avg_r
        naive_ent[loci1_entropy] = naive_overlap

        print(f"entropy: {loci1_entropy} | num_labels: {loci_1.shape[1]-3} | repr_score: {general_rep_score}\n\n")
        
        if loci_1.shape[1]-3 == 2:
            break

        loci_1, loci_2 = merge_clusters(joint, loci_1, loci_2, r1=True)
        loci_1, loci_2 = merge_clusters(joint, loci_1, loci_2, r1=False)
    
    with open('{}/naive_post_clustering_progress.txt'.format(savedir), 'w') as savefile:
        savefile.write(str(naive_rec))
        savefile.write("\n")
        savefile.write(str(naive_ent))
        savefile.write("\n")
    
    with open('{}/MI_post_clustering_progress.txt'.format(savedir), 'w') as savefile:
        savefile.write(str(MI_rec))
        savefile.write("\n")
        savefile.write(str(MI_rec_entropy))
        savefile.write("\n")

    with open('{}/post_clustering_progress.txt'.format(savedir), 'w') as savefile:
        savefile.write(str(nmi_rec))
        savefile.write("\n")
        savefile.write(str(robust_rec))
        savefile.write("\n")
        savefile.write(str(avgr_rec))
        savefile.write("\n")


        savefile.write(str(robust_rec_entropy))
        savefile.write("\n")
        savefile.write(str(NMI_rec_entropy))
        savefile.write("\n")
        savefile.write(str(avgr_entropy))

    ####################################################################################################

    nl = list(robust_rec.keys())
    ys =[robust_rec[k] for k in nl]
    plt.bar(list(nl), list(ys), color="grey")
    plt.xlabel("Number of Labels")
    plt.ylabel("ratio confident")
    plt.yticks(np.arange(0,1.05,0.1))
    plt.xticks(rotation=90)
    plt.savefig('{}/conf_progress.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/conf_progress.svg'.format(savedir), format='svg')
    
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    ####################################################################################################

    nl = list(robust_rec_entropy.keys())
    ys = [robust_rec_entropy[k] for k in nl]
    plt.scatter(list(nl), list(ys), color="Black")
    plt.plot(list(nl), list(ys), "--", color="red", linewidth=3)

    keyss = list(robust_rec.keys())
    for i in range(len(nl)):
        plt.annotate(f"{keyss[i]}", (nl[i], ys[i] + 0.1), rotation=90, fontsize=6)
    
    plt.xlabel("Entropy of Base Annotation")
    plt.ylabel("ratio confident")
    plt.yticks(np.arange(0,1.05,0.1))
    plt.xticks(rotation=90)
    plt.savefig('{}/conf_progress_entropy.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/conf_progress_entropy.svg'.format(savedir), format='svg')
    
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    ####################################################################################################

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
            return

def compare_corresp_methods(replicate_1_dir, replicate_2_dir, outdir, saga="chmm"):
    indicator_file = outdir+"/corresp.pdf"
    if os.path.exists(indicator_file):
        return

    cos_mat = correspondence_based_on_emission(replicate_1_dir, replicate_2_dir, outdir=outdir, saga=saga, metric="cosine")

    ####################################################################################
    loci1, loci2 = load_data(
        f"{replicate_1_dir}/parsed_posterior.csv",
        f"{replicate_2_dir}/parsed_posterior.csv",
        subset=True, logit_transform=False)
    
    loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False, custom_order=False)

    cos_mat = pd.DataFrame(cos_mat, columns=iou.columns, index=iou.index)
    a = Hungarian_algorithm(cos_mat)
    a = [(iou.index[i], iou.columns[j]) for i, j in a]

    iou = IoU_overlap(loci1, loci2)
    c = Hungarian_algorithm(iou)
    c = [(iou.index[i], iou.columns[j]) for i, j in c]
    ####################################################################################
    
    boundaries = [x for x in list(np.linspace(0, 1, 20))] + [1] # custom boundaries
    hex_colors = sns.light_palette('navy', n_colors=len(boundaries) * 2 + 2, as_cmap=False).as_hex()
    hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]
    colors=list(zip(boundaries, hex_colors))
    custom_color_map = LinearSegmentedColormap.from_list(
        name='custom_navy',
        colors=colors)
    
    ####################################################################################
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    sns.heatmap(
        cos_mat.astype(float), annot=False, fmt=".1f",
        linewidths=0.01,  cbar=False, annot_kws={"size": 6}, 
        cmap=custom_color_map, ax=ax1)
    
    sns.heatmap(
        iou.astype(float), annot=False, fmt=".1f",
        linewidths=0.01,  cbar=True, annot_kws={"size": 6}, 
        cmap=custom_color_map, ax=ax2)

    ax1.set_title("Cosine Similarity")
    ax2.set_title("IoU Overlap")

    for i, j in c:
        column_index = cos_mat.columns.get_loc(j)
        row_index = cos_mat.index.get_loc(i)

        # if (i, j) in c:
        ax1.annotate('X', xy=(column_index+0.5, row_index+0.5), ha='center', va='center')
        ax2.annotate('X', xy=(column_index+0.5, row_index+0.5), ha='center', va='center')
        # else:
        #     ax1.annotate('O', xy=(column_index+0.5, row_index+0.5), ha='center', va='center')


    # for i, j in c:
    #     column_index = iou.columns.get_loc(j)
    #     row_index = iou.index.get_loc(i)

    #     # if (i, j) in a:
    #     ax2.annotate('X', xy=(column_index+0.5, row_index+0.5), ha='center', va='center')
    #     # else:
    #     #     ax2.annotate('O', xy=(column_index+0.5, row_index+0.5), ha='center', va='center')

    plt.tight_layout()
    plt.savefig('{}/corresp.pdf'.format(outdir), format='pdf')
    plt.savefig('{}/corresp.svg'.format(outdir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

def r_value_hists(replicate_1_dir, replicate_2_dir, savedir, locis=False, w=1000, to=0.75):
    indicator_file = savedir+"/rval_hist.pdf"
    if os.path.exists(indicator_file):
        return

    if locis:
        loci1, loci2 = replicate_1_dir, replicate_2_dir
    else: 
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=False)

        loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)
    
    ###################################################################################################################
    def ecdf(data):
        """Compute ECDF for a one-dimensional array of measurements."""
        n = len(data)
        x = np.sort(data)
        y = np.arange(1, n+1) / n
        return x, y

    rvalues = is_repr_posterior(
            loci1, loci2, ovr_threshold=to, window_bp=w, matching="static",
            always_include_best_match=True, return_r=True)
        
    labels = rvalues.MAP.unique()
    fig, ax = plt.subplots()

    for l in range(len(labels)):
        data = rvalues.loc[rvalues["MAP"] == labels[l], "r_value"]
        x, y = ecdf(data)
        ax.plot(x, y, label=labels[l])

    ax.set_xlabel("r_value")
    ax.set_ylabel("ECDF")
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig('{}/rval_ecdf.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/rval_ecdf.svg'.format(savedir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

    
    ###################################################################################################################
    rvalues = is_repr_posterior(
        loci1, loci2, ovr_threshold=to, window_bp=w, matching="static",
         always_include_best_match=True, return_r=True)
    
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
        # axs[l].set_xlabel("r_value")
        # axs[l].set_ylabel("Proportion")
        axs[l].tick_params(axis='y', labelsize=8)
        axs[l].set_xlim([0, 1])

    plt.tight_layout()
    plt.savefig('{}/rval_hist.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/rval_hist.svg'.format(savedir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
    plt.clf()

#######################################################################################################################

def GET_ALL(replicate_1_dir, replicate_2_dir, genecode_dir, savedir, rnaseq=None, contour=False):
    print(replicate_1_dir, replicate_2_dir, genecode_dir, savedir)
    if os.path.exists(savedir)==False:
        os.mkdir(savedir)

    try:
        get_all_ct(replicate_1_dir, replicate_2_dir, savedir)
    except:
        pass

    try:
        get_all_labels(replicate_1_dir, replicate_2_dir, savedir)
    except:
        pass

    try:
        get_overalls(replicate_1_dir, replicate_2_dir, savedir, tr=0.9)
    except:
        pass

    try:
        if "chmm" in replicate_1_dir.lower() or "chromhmm_runs" in replicate_1_dir.lower():
            compare_corresp_methods(replicate_1_dir, replicate_2_dir, savedir, saga="chmm")

        elif "segway" in replicate_1_dir.lower():
            compare_corresp_methods(replicate_1_dir, replicate_2_dir, savedir, saga="segway")
            
    except:
        pass

    try:
        r_value_hists(replicate_1_dir, replicate_2_dir, savedir)
    except:
        pass

    try:
        post_clustering(replicate_1_dir, replicate_2_dir, savedir, tr=0.9)

    except:
        pass

    try:
        after_SAGAconf_metrics(replicate_1_dir, replicate_2_dir, genecode_dir, savedir, rnaseq=None, tr=0.9)
        before_after_saga(savedir)

    except:
        pass
    
    ################################################################################################################

    try:
        get_all_bioval(
            replicate_1_dir, replicate_2_dir, 
            savedir,
            genecode_dir=genecode_dir, 
            rnaseq=rnaseq)
    except:
        pass

    if contour:
        try:
            get_contour(replicate_1_dir, replicate_2_dir, savedir)
        except:
            pass

    try:
        gather_labels(replicate_1_dir, savedir, contour=contour)

    except:
        pass

def test_new_functions(replicate_1_dir, replicate_2_dir, genecode_dir, savedir):
    # plot_raw_and_iou(replicate_1_dir, replicate_2_dir)
    post_clustering(replicate_1_dir, replicate_2_dir, savedir, locis=False, to=0.75, tr=0.9)
    exit()
    loci_1, loci_2 = load_data(
        replicate_1_dir+"/parsed_posterior.csv",
        replicate_2_dir+"/parsed_posterior.csv",
        subset=True, logit_transform=False)

    loci_1, loci_2 = process_data(loci_1, loci_2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False, custom_order=True) 
   
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

def heatmaps_on_active_regions(replicate_1_dir, replicate_2_dir, cCREs_file, Meuleman_file, savedir, locis=True, w=1000):
    if locis:
        loci1, loci2 = replicate_1_dir, replicate_2_dir
    else:
        loci1, loci2 = load_data(
            replicate_1_dir+"/parsed_posterior.csv",
            replicate_2_dir+"/parsed_posterior.csv",
            subset=True, logit_transform=False)

        loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)
    ##################################################################################################################

    if os.path.exists(f"{savedir}/activeregions/") == False:
        os.mkdir(f"{savedir}/activeregions/")
    savedir = f"{savedir}/activeregions/"

    MAP1 = loci1.iloc[:,3:].idxmax(axis=1)
    MAP1 = pd.concat([loci1.chr, loci1.start, loci1.end, MAP1], axis=1)
    MAP1.columns = ["chr", "start", "end", "MAP"]

    MAP2 = loci2.iloc[:,3:].idxmax(axis=1)
    MAP2 = pd.concat([loci2.chr, loci2.start, loci2.end, MAP2], axis=1)
    MAP2.columns = ["chr", "start", "end", "MAP"]

    bed_meuleman = pd.read_csv(Meuleman_file, sep="\t")
    bed_meuleman.columns = ["chr", "start", "end", "identifier", "mean_signal", "numsamples", "summit", "core_start", "core_end", "component"]

    bed_ccre = pybedtools.BedTool(cCREs_file)

    bed_meuleman = pybedtools.BedTool.from_dataframe(bed_meuleman)
    bed_MAP1 = pybedtools.BedTool.from_dataframe(MAP1)
    bed_MAP2 = pybedtools.BedTool.from_dataframe(MAP2)

    # Get the intersection
    MAP1_intersect_ccre = bed_MAP1.intersect(bed_ccre, wa=True, wb=True).to_dataframe()
    MAP1_intersect_meul = bed_MAP1.intersect(bed_meuleman, wa=True, wb=True).to_dataframe()
    MAP1_intersect_meul.columns = ["chr","start","end","MAP","_chr", "_start", "_end", "identifier", "mean_signal", "numsamples", "summit", "core_start", "core_end", "component"]

    MAP2_intersect_ccre = bed_MAP2.intersect(bed_ccre, wa=True, wb=True).to_dataframe()
    MAP2_intersect_meul = bed_MAP2.intersect(bed_meuleman, wa=True, wb=True).to_dataframe()
    MAP2_intersect_meul.columns = ["chr","start","end","MAP","_chr", "_start", "_end", "identifier", "mean_signal", "numsamples", "summit", "core_start", "core_end", "component"]

    MAP1_intersect_ccre = MAP1_intersect_ccre[["chrom", "start", "end", "name"]]
    MAP1_intersect_ccre.columns = ["chr", "start", "end", "MAP"]
    
    MAP1_intersect_meul = MAP1_intersect_meul[["chr", "start", "end", "MAP"]]
    
    MAP2_intersect_ccre = MAP2_intersect_ccre[["chrom", "start", "end", "name"]]
    MAP2_intersect_ccre.columns = ["chr", "start", "end", "MAP"]
    
    MAP2_intersect_meul = MAP2_intersect_meul[["chr", "start", "end", "MAP"]]

    def IoU_from_MAP(MAP1_df, MAP2_df, w):
        num_labels1 = len(MAP1_df["MAP"].unique())
        num_labels2 = len(MAP2_df["MAP"].unique())

        states1 = list(MAP1_df["MAP"].unique())
        states2 = list(MAP2_df["MAP"].unique())

        resolution = MAP1_df["end"][0] - MAP1_df["start"][0] 
        w = math.ceil(w/resolution)

        observed_overlap = {} 

        for k in states1:
            for j in states2:
                observed_overlap[str(k + "|" + j)] = 0

        oo_mat = pd.DataFrame(np.zeros((num_labels1, num_labels2)), columns=states2, index=states1)
        IoU = pd.DataFrame(np.zeros((num_labels1, num_labels2)), columns=states2, index=states1)

        def is_conseq(MAP1_df, i, w, resolution, upstream=False):
            statement1 = bool(MAP1_df["chr"][i] == MAP1_df["chr"][max(0, i-w)] == MAP1_df["chr"][min(i+w_i_d, len(MAP1_df)-1)])

            if upstream:
                statement2 = bool( int(MAP1_df["start"][i] - MAP1_df["start"][max(0, i-w)]) == int(w * resolution) )

            else:
                statement2 = bool( int(MAP1_df["start"][min(i+w, len(MAP1_df)-1)] - MAP1_df["start"][i]) == int(w * resolution))

            if statement1 and statement2:
                return True
            else:
                return False

        MAP1 = list(MAP1_df["MAP"])
        MAP2 = list(MAP2_df["MAP"])
                    
        if w == 0:
            for i in range(len(MAP1)):
                k = MAP1[i]
                j = MAP2[i]
                observed_overlap[str(k + "|" + j)] += 1
            
        else:
            for i in range(len(MAP1)):
                k = MAP1[i]
                w_i_u = w
                w_i_d = w
                
                while not is_conseq(MAP1_df, i, w_i_u, resolution, upstream=True):
                    w_i_u -= 1
                    
            
                while not is_conseq(MAP1_df, i, w_i_d, resolution, upstream=False):
                    w_i_d -= 1

                i_neighbors = MAP2[max(0, i-w_i_u) : min(i+w_i_d, len(MAP2)-1)]
                            
                for j in set(i_neighbors):
                    observed_overlap[str(k + "|" + j)] += 1


        for p in observed_overlap.keys():
            oo_mat.loc[p.split("|")[0], p.split("|")[1]] = observed_overlap[p]
        
        joint = oo_mat / len(MAP1)
        MAP1 = MAP1_df["MAP"]
        MAP2 = MAP2_df["MAP"]

        coverage1 = {k:len(MAP1.loc[MAP1 == k])/len(MAP1) for k in states1}
        coverage2 = {k:len(MAP2.loc[MAP2 == k])/len(MAP2) for k in states2}

        for A in states1:
            for B in states2:
                    IoU.loc[A, B] = (joint.loc[A, B])/(coverage1[A] + coverage2[B] - ((joint.loc[A, B])))

        return IoU

    iou_0_ccre = IoU_from_MAP(MAP1_intersect_ccre, MAP2_intersect_ccre, w=0)
    iou_1k_ccre = IoU_from_MAP(MAP1_intersect_ccre, MAP2_intersect_ccre, w=1000)

    iou_0_meul = IoU_from_MAP(MAP1_intersect_meul, MAP2_intersect_meul, w=0)
    iou_1k_meul = IoU_from_MAP(MAP1_intersect_meul, MAP2_intersect_meul, w=1000)

    iou_0 = IoU_from_MAP(MAP1, MAP2, w=0)
    iou_1k = IoU_from_MAP(MAP1, MAP2, w=1000)

    def plot_heatmap(confmat1, confmat2, savedir, sufix, w1, w2):
        boundaries = [x**2 for x in list(np.linspace(0, 1, 16))] + [1] # custom boundaries
        hex_colors = sns.light_palette('navy', n_colors=len(boundaries) * 2 + 2, as_cmap=False).as_hex()
        hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]
        colors=list(zip(boundaries, hex_colors))
        custom_color_map = LinearSegmentedColormap.from_list(
            name='custom_navy',
            colors=colors)

        # Calculate the Mean Absolute Error (MAE) between the two matrices
        mae = np.mean(np.abs(confmat1 - confmat2))

        # Adjust the figure size here
        fig, axs = plt.subplots(1, 2, figsize=(12, 6)) 
        sns.set()

        for ax, confmat, w in zip(axs, [confmat1, confmat2], [w1, w2]):
            p = sns.heatmap(
                confmat.astype(float), annot=True, fmt=".2f",
                linewidths=0.01, cbar=True, annot_kws={"size": 6}, cmap=custom_color_map,
                ax=ax)

            ax.tick_params(axis='x', rotation=90, labelsize=10)
            ax.tick_params(axis='y', rotation=0, labelsize=10)
            ax.set_title(f'IoU Overlap {sufix} | w = {w} | MAE = {mae:.2f}')
            ax.set_xlabel('Replicate 2 Labels')
            ax.set_ylabel("Replicate 1 Labels")

        plt.tight_layout()
        plt.savefig('{}/IoU_{}_w{}_w{}.svg'.format(savedir, sufix, str(w1), str(w2)), format='svg')
        plt.savefig('{}/IoU_{}_w{}_w{}.pdf'.format(savedir, sufix, str(w1), str(w2)), format='pdf')
        sns.reset_orig()
        plt.close("all")
        plt.style.use('default')
        plt.clf()

    # Now you can call the function like this:
    plot_heatmap(iou_0_ccre, iou_1k_ccre, savedir, sufix="cCRE", w1=0, w2=1000)
    plot_heatmap(iou_0_meul, iou_1k_meul, savedir, sufix="Meuleman", w1=0, w2=1000)
    plot_heatmap(iou_0, iou_1k, savedir, sufix="WG", w1=0, w2=1000)

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
    
