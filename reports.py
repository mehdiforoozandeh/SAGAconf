from run import *
from indepth import *
from granul import *
from bio_valid import *
from overall import *
from matplotlib.colors import LinearSegmentedColormap
import ast
from scipy.interpolate import UnivariateSpline


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
    with open("{}/raw_conditional_overlap_ratio.txt".format(savedir), "w") as cf:
        overall_overlap = overall_overlap_ratio(loci_1, loci_2, w=0)
        cf.write("{} : {}\n".format("overall w=0", overall_overlap))

        overall_overlap = overall_overlap_ratio(loci_1, loci_2, w=1000)
        cf.write("{} : {}\n".format("overall w=1000", overall_overlap))

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
    calibrated_loci_1 = calb.perlabel_calibration_function()

    if os.path.exists(pltsavedir+"/calib_logit_enr") == False:
        os.mkdir(pltsavedir+"/calib_logit_enr")
    calb = posterior_calibration(
        loci_1, loci_2, plot_raw=False, window_size=1000, savedir=pltsavedir+"/calib_logit_enr", allow_w=False)
    calibrated_loci_1 = calb.perlabel_calibration_function()
    
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
    reproduced_loci1.to_csv(savedir + "/confident_posteriors.bed", sep='\t', header=True, index=False)

    ##################################################################################################################

    bool_reprod_report = pd.concat(
        [loci1["chr"], loci1["start"], loci1["end"], pd.Series(bool_reprod_report)], axis=1)
    bool_reprod_report.columns = ["chr", "start", "end", "is_repr"]
        
    bool_reprod_report.to_csv(savedir+"/boolean_reproducibility_report_POSTERIOR.csv")

    general_rep_score = len(bool_reprod_report.loc[bool_reprod_report["is_repr"]==True]) / len(bool_reprod_report)
    perlabel_rec = {}
    lab_rep = perlabel_is_reproduced(bool_reprod_report, loci1)
    for k, v in lab_rep.items():
        perlabel_rec[k] = v[0]/ (v[0]+v[1])

    with open(savedir+"/general_reproducibility_score_POSTERIOR.txt", "w") as scorefile:
        scorefile.write("general reprod score = {}\n".format(general_rep_score))
        for k, v in perlabel_rec.items():
            scorefile.write("{} reprod score = {}\n".format(k, str(v)))

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
    
    with open(savedir+"/r_values.txt", "w") as scorefile:
        scorefile.write("average r_value = {}\n".format(avg_r))
        for k, v in perlabel_r.items():
            scorefile.write("{} r_value = {}\n".format(k, str(v)))

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
    indicator_file = '{}/conf_rogress.pdf'.format(savedir)
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
        loci1_entropy = calc_entropy(joint, rows=True)

        nmi_rec["{}".format(loci_1.shape[1]-3)] = NMI
        robust_rec["{}".format(loci_1.shape[1]-3)] = general_rep_score
        avgr_rec["{}".format(loci_1.shape[1]-3)] = avg_r

        robust_rec_entropy[loci1_entropy] = general_rep_score
        NMI_rec_entropy[loci1_entropy] = NMI
        avgr_entropy[loci1_entropy] = avg_r

        print(f"entropy: {loci1_entropy} | num_labels: {loci_1.shape[1]-3} | repr_score: {general_rep_score}\n\n")
        
        if loci_1.shape[1]-3 == 2:
            break

        loci_1, loci_2 = merge_clusters(joint, loci_1, loci_2, r1=True)
        loci_1, loci_2 = merge_clusters(joint, loci_1, loci_2, r1=False)

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

def post_clustering_keep_k_states(replicate_1_dir, replicate_2_dir, savedir, k, locis=False, write_csv=True):

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
                loci_1.to_csv(savedir + f"/{str(k)}states_post_clustered_posterior.csv")
            loci_1.to_csv(savedir + f"/{str(k)}states_post_clustered_posterior.bed", sep='\t', header=True, index=False)

            MAP = loci_1.iloc[:,3:].idxmax(axis=1)
            coordMAP = pd.concat([loci_1.iloc[:, :3], MAP], axis=1)
            coordMAP.columns = ["chr", "start", "end", "MAP"]
            denseMAP = condense_segments(coordMAP)
            if write_csv:
                denseMAP.to_csv(savedir + f"/{str(k)}_states_confident_segments_dense.csv")
            denseMAP.to_csv(savedir + f"/{str(k)}_states_confident_segments_dense.bed", sep='\t', header=True, index=False)
            return

def compare_corresp_methods(replicate_1_dir, replicate_2_dir, outdir, saga="chmm"):
    loci1, loci2 = load_data(
        f"{replicate_1_dir}/parsed_posterior.csv",
        f"{replicate_2_dir}/parsed_posterior.csv",
        subset=True, logit_transform=False)
    
    loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False, custom_order=False)

    cos_mat = correspondence_based_on_emission(replicate_1_dir, replicate_2_dir, outdir=outdir, saga=saga, metric="cosine")

    iou = IoU_overlap(loci1, loci2)
    cos_mat = pd.DataFrame(cos_mat, columns=iou.columns, index=iou.index)
    
    a = Hungarian_algorithm(cos_mat)
    a = [(iou.index[i], iou.columns[j]) for i, j in a]
        
    c = Hungarian_algorithm(iou)
    c = [(iou.index[i], iou.columns[j]) for i, j in c]

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
        if "chmm" in replicate_1_dir.lower() or "chromhmm" in replicate_1_dir.lower():
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

    post_clustering(replicate_1_dir, replicate_2_dir, savedir, locis=False, to=0.75, tr=0.9)
    exit()
    
    # r_value_hists(replicate_1_dir, replicate_2_dir, savedir, locis=False, w=1000, to=0.75)
    if "chmm" in replicate_1_dir:
        compare_corresp_methods(replicate_1_dir, replicate_2_dir, savedir, saga="chmm")
    elif "segway" in replicate_1_dir:
        compare_corresp_methods(replicate_1_dir, replicate_2_dir, savedir, saga="segway")

    # # post_clustering(replicate_1_dir, replicate_2_dir, savedir, locis=False, to=0.75, tr=0.8)
    # # get_overalls(replicate_1_dir, replicate_2_dir, savedir, locis=False, w=1000)

    # loci_1, loci_2 = load_data(
    #     replicate_1_dir+"/parsed_posterior.csv",
    #     replicate_2_dir+"/parsed_posterior.csv",
    #     subset=True, logit_transform=False)

    # loci_1, loci_2 = process_data(loci_1, loci_2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)
    # distance_vs_overlap_3(loci_1, loci_2, savedir, match_definition="BM")

    # print(loci1, loci2)

    # post_clustering_keep_k_states(replicate_1_dir, replicate_2_dir, replicate_1_dir, k=5, locis=False)

    # print(NMI_from_matrix(joint_prob_with_binned_posterior(loci1, loci2, n_bins=50, conditional=False, stratified=True)))
    # print(NMI_from_matrix(joint_prob_MAP_with_posterior(loci1, loci2, n_bins=50, conditional=False, stratified=True)))

    # after_SAGAconf_metrics(replicate_1_dir, replicate_2_dir, genecode_dir, savedir, rnaseq=None)
    # before_after_saga(savedir)

    # boundaries = [x for x in list(np.linspace(0, 1, 20))] + [1] # custom boundaries
    # hex_colors = sns.light_palette('navy', n_colors=len(boundaries) * 2 + 2, as_cmap=False).as_hex()
    # hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]
    # colors=list(zip(boundaries, hex_colors))
    # custom_color_map = LinearSegmentedColormap.from_list(
    #     name='custom_navy',
    #     colors=colors)
    
    # ####################################################################################
    # confmat = joint_overlap_prob(loci_1, loci_2, w=0, symmetric=True)

    # p = sns.heatmap(
    #     confmat.astype(float), annot=True, fmt=".2f",
    #     linewidths=0.01,  cbar=True, annot_kws={"size": 8}, 
    #     vmin=0, vmax=1, cmap=custom_color_map)

    # sns.set(rc={'figure.figsize':(20,15)})
    # p.tick_params(axis='x', rotation=90, labelsize=10)
    # p.tick_params(axis='y', rotation=0, labelsize=10)

    # plt.title('IoU Overlap')
    # plt.xlabel('Replicate 2 Labels')
    # plt.ylabel("Replicate 1 Labels")
    # plt.tight_layout()
    # plt.show()

if __name__=="__main__":  
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
    
