from run import *
from indepth import *
from granul import *
from bio_valid import *
from overall import *

def load_data(posterior1_dir, posterior2_dir, subset=False, logit_transform=False):
    print("loading and intersecting")
    loci_1, loci_2 = intersect_parsed_posteriors(
        posterior1_dir, 
        posterior2_dir)

    if subset:
        loci_1 = loci_1.loc[loci_1["chr"]=="chr21"].reset_index(drop=True)
        loci_2 = loci_2.loc[loci_2["chr"]=="chr21"].reset_index(drop=True)

    print("the shapes of the input matrices are: {}, {}".format(str(loci_1.shape), str(loci_2.shape)))

    if logit_transform:
        loci_1.iloc[:,3:] = logit_array(np.array(loci_1.iloc[:,3:]))
        loci_2.iloc[:,3:] = logit_array(np.array(loci_2.iloc[:,3:]))

    return loci_1, loci_2

def process_data(loci_1, loci_2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False):
    print('generating confmat 1 ...')
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
            loci_1_mnemon = read_mnemonics("/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics.txt")
            loci_2_mnemon = read_mnemonics("/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics.txt")

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
                mnemon1_dict[str(i)] = str(i)
            if str(i) not in mnemon2_dict.keys():
                mnemon2_dict[str(i)] = str(i)
        
        if match:
            conf_mat = confusion_matrix(
                loci_1, loci_2, num_labels, 
                OE_transform=True, symmetric=False)

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
            conf_mat = confusion_matrix(
                loci_1, loci_2, num_labels, 
                OE_transform=True, symmetric=False)

            assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')
            loci_1, loci_2 = \
                connect_bipartite(loci_1, loci_2, assignment_pairs, mnemon=False)

            print('connected barpartite')

    return loci_1, loci_2

def ct_confus(loci_1, loci_2, savedir):
    """
    labels can be matched or not
    """
    num_labels = loci_1.shape[1]-3
        
    confmat = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=False)

    p = sns.heatmap(
        confmat.astype(int), annot=True, fmt="d",
        linewidths=0.01,  cbar=False)

    sns.set(rc={'figure.figsize':(15,20)})
    p.tick_params(axis='x', rotation=30, labelsize=7)
    p.tick_params(axis='y', rotation=30, labelsize=7)

    plt.title('Label Matching Heatmap (log(O/E) overlap)')
    plt.xlabel('Replicate 1 Labels')
    plt.ylabel("Replicate 2 Labels")
    plt.tight_layout()
    plt.savefig('{}/heatmap.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/heatmap.svg'.format(savedir), format='svg')
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    confmat.to_csv("{}/heatmap.csv".format(savedir))

def ct_granul(loci_1, loci_2, savedir):
    """
    for this function, the labels should not be matched
    """
    
    num_labels = loci_1.shape[1]-3
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[25, 16])
    label_being_plotted = 0

    for i in range(n_rows):
        for j in range(n_cols):
            c = list(loci_1.columns[3:])[label_being_plotted]
            cr, ar, ovr_rec = granularity_vs_agreement_nonsymmetric(loci_1.copy(), loci_2.copy(), k=c)

            perfect_agr = [0] + [1 for i in range(len(ar) - 1)]
            realAUC = metrics.auc(cr, ar)
            perfectAUC = metrics.auc(cr, perfect_agr)
            p_to_r_auc = float(realAUC/perfectAUC)

            axs[i,j].plot(cr, ar, c="yellowgreen")
            axs[i,j].set_title(
                c + str(" | Real/Perfect AUC = {:.2f}".format(p_to_r_auc)), 
                fontsize=11)

            axs[i,j].fill_between(cr, ar, color="yellowgreen", alpha=0.4)
            axs[i,j].fill_between(cr, perfect_agr, ar, color="palevioletred", alpha=0.4)
            axs[i,j].set_xticks(np.arange(0, 1.1, step=0.2))
            axs[i,j].set_yticks(np.arange(0, 1.1, step=0.2))
            label_being_plotted+=1
    
    # fig.text(0.5, 0.02, 'Coverage' , ha='center')
    # fig.text(0.02, 0.5, 'Agreement', va='center', rotation='vertical')
    plt.tight_layout()
    plt.savefig(savedir+"/granularity.pdf", format='pdf')
    plt.savefig(savedir+"/granularity.svg", format='svg')
    
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')

def ct_lable_calib(loci_1, loci_2, pltsavedir):
    """
    labels need to be matched
    """
    if os.path.exists(pltsavedir+"/calib") == False:
        os.mkdir(pltsavedir+"/calib")
    calb = posterior_calibration(
        loci_1, loci_2, log_transform=False, ignore_overconf=False, filter_nan=True, 
        oe_transform=True, savedir=pltsavedir+"/calib")
    calibrated_loci_1 = calb.perlabel_calibration_function(
        degree=5, num_bins=25, return_caliberated_matrix=True, scale_columnwise=True)
    
    plt.close("all")
    plt.style.use('default')

def ct_boundar(loci_1, loci_2, outdir, match_definition="BM", max_distance=50):
    """
    Here i'm gonna merge ECDF of length dist and boundary thing.
    """

    resolution = int(loci_1.iloc[0, 2] - loci_1.iloc[0, 1])

    max_distance = int(10000/resolution)
    num_labels = loci_1.shape[1]-3
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
    
    """
    we need to define what a "good match" is.
    by default, we will consider any label with log(OO/EO)>0 
    to be a match. this definition can be refined later.
    """

    confmat = enrichment_of_overlap_matrix(
            loci_1, loci_2, OE_transform=True)
    
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
        if MAP1[i] != MAP1[i+1]: 
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

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[16, 9])
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
            
            lendist = get_ecdf_for_label(loci_1, k, max=(max_distance+1)*resolution)
            if len(lendist)>0:
                sns.ecdfplot(lendist, ax=axs[i, j], color="red")

            axs[i, j].set_xticks(np.arange(0, (max_distance+1)*resolution, step=5*resolution))
            axs[i, j].tick_params(axis='both', labelsize=7)
            axs[i, j].tick_params(axis='x', rotation=90)
            axs[i, j].set_yticks(np.arange(0, 1.1, step=0.2))
            axs[i, j].set_xlabel('bp', fontsize=7)
            axs[i, j].set_ylabel('Proportion',fontsize=7)
            
            
            axs[i, j].set_title(k, fontsize=7)

            label_being_plotted += 1
            
    plt.tight_layout()
    plt.savefig(outdir+"/len_bound.pdf", format='pdf')
    plt.savefig(outdir+"/len_bound.svg", format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')

def label_granularity(loci_1, loci_2, savedir):
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
            fontsize=11)

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

def label_boundary(loci_1, loci_2, savedir, match_definition="BM", max_distance=50):
    savedir = savedir+"/len_bound/"
    if os.path.exists(savedir)==False:
        os.mkdir(savedir)

    """
    Here i'm gonna merge ECDF of length dist and boundary thing.
    """

    resolution = int(loci_1.iloc[0, 2] - loci_1.iloc[0, 1])

    max_distance = int(10000/resolution)
    num_labels = loci_1.shape[1]-3
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
    
    """
    we need to define what a "good match" is.
    by default, we will consider any label with log(OO/EO)>0 
    to be a match. this definition can be refined later.
    """

    confmat = enrichment_of_overlap_matrix(
            loci_1, loci_2, OE_transform=True)
    
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
        if MAP1[i] != MAP1[i+1]: 
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
        
        lendist = get_ecdf_for_label(loci_1, k, max=(max_distance+1)*resolution)
        if len(lendist)>0:
            sns.ecdfplot(lendist, color="red")

        plt.xticks(np.arange(0, (max_distance+1)*resolution, step=5*resolution))
        plt.tick_params(axis='both', labelsize=7)
        plt.tick_params(axis='x', rotation=90)
        plt.yticks(np.arange(0, 1.1, step=0.2))
        plt.xlabel('bp', fontsize=7)
        plt.ylabel('Proportion',fontsize=7)
        
        plt.title(k, fontsize=7)
            
        plt.tight_layout()
        plt.savefig(savedir+"/len_bound_{}.pdf".format(k), format='pdf')
        plt.savefig(savedir+"/len_bound_{}.svg".format(k), format='svg')

        sns.reset_orig
        plt.close("all")
        plt.style.use('default')

def label_merging_progression(loci_1, loci_2, savedir):
    savedir = savedir+"/prog/"
    if os.path.exists(savedir)==False:
        os.mkdir(savedir)

    for c in list(loci_1.columns[3:]):
        cr, ar, ovr_rec = granularity_vs_agreement_nonsymmetric(loci_1.copy(), loci_2.copy(), k=c)
        plot_progression(ar, cr, ovr_rec, c, savedir)

        sns.reset_orig
        plt.close("all")
        plt.style.use('default')

def get_all_labels(replicate_1_dir, replicate_2_dir):
    """
    -single granularity
    -single prog
    -single boundar
    """
    loci1, loci2 = load_data(
        replicate_1_dir+"/parsed_posterior.csv",
        replicate_2_dir+"/parsed_posterior.csv",
        subset=True)

    loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    label_granularity(loci1, loci2, replicate_1_dir)
    label_merging_progression(loci1, loci2, replicate_1_dir)
    label_boundary(loci1, loci2, replicate_1_dir, match_definition="BM", max_distance=50)

    label_granularity(loci2, loci1, replicate_2_dir)
    label_merging_progression(loci2, loci1, replicate_2_dir)
    label_boundary(loci2, loci1, replicate_2_dir, match_definition="BM", max_distance=50)

def get_all_ct(replicate_1_dir, replicate_2_dir):
    loci1, loci2 = load_data(
        replicate_1_dir+"/parsed_posterior.csv",
        replicate_2_dir+"/parsed_posterior.csv",
        subset=True)

    loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=True)
    ct_confus(loci1, loci2, replicate_1_dir)
    ct_confus(loci2, loci1, replicate_2_dir)

    ct_lable_calib(loci1, loci2, replicate_1_dir)
    ct_lable_calib(loci2, loci1, replicate_2_dir)

    loci1, loci2 = load_data(
        replicate_1_dir+"/parsed_posterior.csv",
        replicate_2_dir+"/parsed_posterior.csv",
        subset=True)
    loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    ct_granul(loci1, loci2, replicate_1_dir)
    ct_granul(loci2, loci1, replicate_2_dir)

    ct_boundar(loci1, loci2, replicate_1_dir, match_definition="BM", max_distance=50)
    ct_boundar(loci2, loci1, replicate_2_dir, match_definition="BM", max_distance=50)

def gather_labels(savedir):
    
    print("loading mnemonics...")
    if os.path.exists(
        "/".join(savedir.split("/")[:-1])+"/mnemonics_rep1.txt"):

        print("reading concat mnemonics")
        loci_1_mnemon = read_mnemonics("/".join(savedir.split("/")[:-1])+"/mnemonics_rep1.txt")
    else:
        print("reading mnemonics")
        loci_1_mnemon = read_mnemonics("/".join(savedir.split("/")[:-1])+"/mnemonics.txt")

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

def get_all_bioval(replicate_1_dir, replicate_2_dir, genecode_dir, rnaseq=None):
    loci1, loci2 = load_data(
        replicate_1_dir+"/parsed_posterior.csv",
        replicate_2_dir+"/parsed_posterior.csv",
        subset=True, logit_transform=False)

    loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    gene_coords = load_gene_coords(genecode_dir)
    if rnaseq != None:
        trans_data = load_transcription_data(rnaseq, gene_coords)

        trans_data = trans_data.drop(trans_data[trans_data.TPM==0].index).reset_index(drop=True)

        trans_data_exp = trans_data[np.log10(trans_data.TPM) > 2]
        trans_data_notexp = trans_data[np.log10(trans_data.TPM) < 0.5]

        trans_data_exp.TPM = np.log10(trans_data_exp.TPM)
        trans_data_notexp.TPM = np.log10(trans_data_notexp.TPM)
        trans_data.TPM = np.log10(trans_data.TPM)

        plot_general_transc(loci1, trans_data_exp, savedir=replicate_1_dir+"/general_transc_exp", exp=True)
        plot_general_transc(loci1, trans_data_notexp, savedir=replicate_1_dir+"/general_transc_notexp", exp=False)

        plot_general_transc(loci2, trans_data_exp, savedir=replicate_2_dir+"/general_transc_exp", exp=True)
        plot_general_transc(loci2, trans_data_notexp, savedir=replicate_2_dir+"/general_transc_notexp", exp=False)

        posterior_transcription_correlation(loci1, trans_data, savedir=replicate_1_dir+"/trans_post_correl")
        posterior_transcription_correlation(loci2, trans_data, savedir=replicate_2_dir+"/trans_post_correl")

    overal_TSS_enrichment(loci1, replicate_1_dir+"/tss_enr")
    overal_TSS_enrichment(loci2, replicate_2_dir+"/tss_enr")

def get_overalls(replicate_1_dir, replicate_2_dir):
    loci1, loci2 = load_data(
        replicate_1_dir+"/parsed_posterior.csv",
        replicate_2_dir+"/parsed_posterior.csv",
        subset=True, logit_transform=True)

    loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)
    
    bool_reprod_report = is_reproduced(loci1, loci2, enr_threshold=2, window_bp=1000)
    bool_reprod_report.to_csv(replicate_1_dir+"/boolean_reproducibility_report.csv")

    bool_reprod_report = is_reproduced(loci2, loci1, enr_threshold=2, window_bp=1000)
    bool_reprod_report.to_csv(replicate_2_dir+"/boolean_reproducibility_report.csv")

    reprep = reprod_score(loci1, loci2, window_bp=1000)
    pcorrel, scorrel = max_posterior_to_repr_score(loci1, reprep)
    print(pcorrel, scorrel)

    reprep = reprod_score(loci2, loci1, window_bp=1000)
    pcorrel, scorrel = max_posterior_to_repr_score(loci2, reprep)
    print(pcorrel, scorrel)

def get_contour(replicate_1_dir, replicate_2_dir):
    loci1, loci2 = load_data(
        replicate_1_dir+"/parsed_posterior.csv",
        replicate_2_dir+"/parsed_posterior.csv",
        subset=True, logit_transform=True)

    loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)
    contour_isrep(loci1, loci2, replicate_1_dir)
    contour_isrep(loci2, loci1, replicate_2_dir)

def GET_ALL(replicate_1_dir, replicate_2_dir, genecode_dir, rnaseq=None):
    get_all_ct(replicate_1_dir, replicate_2_dir)
    get_all_labels(replicate_1_dir, replicate_2_dir)

    gather_labels(replicate_1_dir)
    gather_labels(replicate_2_dir)

    get_all_bioval(
        replicate_1_dir, replicate_2_dir, 
        genecode_dir=genecode_dir, 
        rnaseq=rnaseq)
    
    get_overalls(replicate_1_dir, replicate_2_dir)

if __name__=="__main__":
    replicate_1_dir = "tests/cedar_runs/chmm/MCF7_R1/"
    replicate_2_dir = "tests/cedar_runs/chmm/MCF7_R2/"
    get_contour(replicate_1_dir, replicate_2_dir)
    exit()

    GET_ALL(
        replicate_1_dir="tests/cedar_runs/segway/GM12878_R1/", 
        replicate_2_dir="tests/cedar_runs/segway/GM12878_R2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv")

    # replicate_1_dir = "tests/cedar_runs/chmm/MCF7_R1/"
    # replicate_2_dir = "tests/cedar_runs/chmm/MCF7_R2/"
    # get_all_ct(replicate_1_dir, replicate_2_dir)
    # get_all_labels(replicate_1_dir, replicate_2_dir)
    # gather_labels(replicate_1_dir)
    # gather_labels(replicate_2_dir)

    # get_all_bioval(
    #     replicate_1_dir, replicate_2_dir, 
    #     genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
    #     rnaseq=None)

    # replicate_1_dir = "tests/cedar_runs/segway/MCF7_R1/"
    # replicate_2_dir = "tests/cedar_runs/segway/MCF7_R2/"
    # get_all_ct(replicate_1_dir, replicate_2_dir)
    # get_all_labels(replicate_1_dir, replicate_2_dir)
    # gather_labels(replicate_1_dir)
    # gather_labels(replicate_2_dir)

    # replicate_1_dir = "tests/cedar_runs/chmm/GM12878_R1/"
    # replicate_2_dir = "tests/cedar_runs/chmm/GM12878_R2/"
    # get_all_ct(replicate_1_dir, replicate_2_dir)
    # get_all_labels(replicate_1_dir, replicate_2_dir)
    # gather_labels(replicate_1_dir)
    # gather_labels(replicate_2_dir)

    # replicate_1_dir = "tests/cedar_runs/segway/GM12878_R1/"
    # replicate_2_dir = "tests/cedar_runs/segway/GM12878_R2/"
    # get_all_ct(replicate_1_dir, replicate_2_dir)
    # get_all_labels(replicate_1_dir, replicate_2_dir)
    # gather_labels(replicate_1_dir)
    # gather_labels(replicate_2_dir)