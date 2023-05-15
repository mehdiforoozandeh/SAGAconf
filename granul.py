import matplotlib.pyplot as plt
from _reproducibility import *
from _cluster_matching import *
from run import *
from sklearn import metrics
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from matplotlib.colors import LinearSegmentedColormap


"""
Part1:
    in this module, we try to skip the post clustering step 
    and we go straight to merging. 

    psdcode:

    def gran(loci_1, loci_2):
        get raw confusion matrix
        get enrichment confusion matrix
        for all labels in r1
            sort the enrichment values from high to low
            start merging them by the order 
            get genome coverage of the new label
            calculate r1vsr1 agreement in that particular label

            plot agreement as a function of genome coverage
"""

"""
Part2:
what if we skip matching too?

generate raw-unmatched logOE overlap vs r2
for  label k in r1:
    sort labels in r2 based on logOE overlap with k
    just merge labels in r2 
    at each step calculate agreement of k with merged labels
    <<in should be a non-decreasing trend>>
"""

def get_coverage(loci, query_label="all"):
    MAP = loci.iloc[:,3:].idxmax(axis=1)
    coverage = dict(zip(list(loci.columns[3:]), [0 for _ in range(len(loci.columns[3:]))]))

    if query_label == "all":
        for c in list(loci.columns[3:]):
            coverage[c] = len(MAP.loc[MAP == c]) / len(MAP)
        return coverage

    else:
        return len(MAP.loc[MAP == query_label]) / len(MAP)

def get_agr_symmetric(MAP1, MAP2, query_label):
    record = [0,0] # [agreement_count, disagreement_count]

    for i in range(len(MAP1)):
        if MAP1[i] == query_label:# or MAP2[i] == query_label:

            if MAP1[i] == MAP2[i]:
                record[0] += 1
            else:
                record[1] += 1

    agr = float(record[0] / (record[0] + record[1]))
    return agr

def get_agr_nonsymmetric(MAP1, MAP2, query_label1, query_label2):
    record = [0,0] # [agreement_count, disagreement_count]

    for i in range(len(MAP1)):
        if MAP1[i] == query_label1:

            if MAP2[i] == query_label2:
                record[0] += 1
            else:
                record[1] += 1

    try:
        agr = float(record[0] / (record[0] + record[1]))
    except:
        agr = 0
    return agr

def granularity_vs_agreement_symmetric(loci_1, loci_2, query_label, disregard_posterior=True):
    '''
    input locis must be matched (after hungarian alg) and optionally with updated mnemonics
    '''
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
    if disregard_posterior:
        # to avoid MAP changes during the merging process, disregard probability granularity and convert 
        # to hard zero vs one matrix
        

        for c in loci_1.columns[3:]:
            loci_1.loc[MAP1==c, c] = 1
            loci_1.loc[MAP1!=c, c] = 0
        
        for c in loci_2.columns[3:]:
            loci_2.loc[MAP2==c, c] = 1
            loci_2.loc[MAP2!=c, c] = 0
    
    k = query_label 

    confmat = overlap_matrix(loci_1, loci_2, type="IoU")
    sorted_k_vector = confmat.loc[k,:].sort_values(ascending=False)
    # print(sorted_k_vector)

    pre_merge_agr = get_agr_symmetric(MAP1, MAP2, query_label=k)
    pre_merge_coverage = get_coverage(loci_1, query_label=k)

    coverage_record = [pre_merge_coverage]
    agreement_record = [pre_merge_agr]

    for j in range(0, len(sorted_k_vector)):

        if str(sorted_k_vector.index[j]) not in query_label: #can't merge sth with itself :))

            loci_1[query_label + "+" + str(sorted_k_vector.index[j])] = \
                loci_1[query_label] + loci_1[str(sorted_k_vector.index[j])]

            loci_1 = loci_1.drop([query_label, str(sorted_k_vector.index[j])], axis=1)

            loci_2[query_label + "+" + str(sorted_k_vector.index[j])] = \
                loci_2[query_label] + loci_2[str(sorted_k_vector.index[j])]

            loci_2 = loci_2.drop([query_label, str(sorted_k_vector.index[j])], axis=1)

            query_label = query_label + "+" + str(sorted_k_vector.index[j])
        
            MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
            MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

            # print(loci_1.iloc[:,3:].max(axis=1).max(), '\n', loci_1.iloc[:,3:].sum(axis=1).sum())
            # print(loci_2.iloc[:,3:].max(axis=1).max(), '\n', loci_2.iloc[:,3:].sum(axis=1).sum())

            jth_order_agr = get_agr_symmetric(MAP1, MAP2, query_label=query_label)
            jth_order_cov = get_coverage(loci_1, query_label=query_label)

            coverage_record.append(jth_order_cov)
            agreement_record.append(jth_order_agr)

    return coverage_record, agreement_record
    
def granularity_vs_agreement_nonsymmetric(loci_1, loci_2, k, disregard_posterior=True):
    '''
    input locis must NOT be matched -- optionally with updated mnemonics
    '''
    num_labels = loci_1.shape[1]-3
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
    if disregard_posterior:
        # to avoid MAP changes during the merging process, disregard probability granularity and convert 
        # to hard zero vs one matrix
        
        for c in loci_1.columns[3:]:
            loci_1.loc[MAP1==c, c] = 1
            loci_1.loc[MAP1!=c, c] = 0
        
        for c in loci_2.columns[3:]:
            loci_2.loc[MAP2==c, c] = 1
            loci_2.loc[MAP2!=c, c] = 0

    confmat = IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False)
    
    sorted_k_vector = confmat.loc[k,:].sort_values(ascending=False)
    # print(k, sorted_k_vector)

    query_label = str(sorted_k_vector.index[0])
    coverage_record = [0]
    agreement_record = [0]

    for j in range(0, len(sorted_k_vector)):
        if j == 0:
            pass
        else:
            loci_2[query_label + "+" + str(sorted_k_vector.index[j])] = \
                loci_2[query_label] + loci_2[str(sorted_k_vector.index[j])]

            loci_2 = loci_2.drop([query_label, str(sorted_k_vector.index[j])], axis=1)

            query_label = query_label + "+" + str(sorted_k_vector.index[j])
    
        MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
        MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

        jth_order_agr = get_agr_nonsymmetric(MAP1, MAP2, query_label1=k, query_label2=query_label)
        jth_order_cov = get_coverage(loci_2, query_label=query_label)

        coverage_record.append(jth_order_cov)
        agreement_record.append(jth_order_agr)
    return coverage_record, agreement_record, sorted_k_vector.index

def plot_progression(ar, cr, ovr_rec, c, savedir):
    ff, aa = plt.subplots(figsize=(8,5))

    perfect_agr = [0] + [1 for i in range(len(ar) - 1)]
    realAUC = metrics.auc(cr, ar)
    perfectAUC = metrics.auc(cr, perfect_agr)
    p_to_r_auc = float(realAUC/perfectAUC)

    aa.plot(cr, ar, "--o" , c="grey")
    aa.set_title(
        c + str(" | Real/Perfect AUC = {:.2f}".format(p_to_r_auc)), 
        fontsize=12)
    
    aa.fill_between(cr, ar, perfect_agr, color="palevioletred", alpha=0.4)
    aa.fill_between(cr, ar, color="yellowgreen", alpha=0.4)
    aa.set_xticks(np.arange(0, 1.1, step=0.2))
    aa.set_yticks(np.arange(0, 1.1, step=0.2))

    ovr_rec = [str(i+1)+") "+ ovr_rec[i] for i in range(len(ovr_rec))]
    textstr = "\n".join(ovr_rec)
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)

    aa.text(1.02, 0.98, textstr, transform=aa.transAxes, fontsize=10, verticalalignment='top', bbox=props)

    plt.tight_layout()
    plt.savefig(savedir+"/prog_{}.pdf".format(c), format='pdf')
    plt.savefig(savedir+"/prog_{}.svg".format(c), format='svg')

    plt.close("all")
    plt.style.use('default')

def run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons, symmetric):
    print("loading and intersecting")
    loci_1, loci_2 = intersect_parsed_posteriors(
        replicate_1_dir+"/parsed_posterior.csv", 
        replicate_2_dir+"/parsed_posterior.csv")

    if run_on_subset:
        loci_1 = loci_1.loc[loci_1["chr"]=="chr1"].reset_index(drop=True)
        loci_2 = loci_2.loc[loci_2["chr"]=="chr1"].reset_index(drop=True)

    print("the shapes of the input matrices are: {}, {}".format(str(loci_1.shape), str(loci_2.shape)))

    num_labels = loci_1.shape[1]-3
    loci_1.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    loci_2.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    
    if symmetric:
        print('generating confmat 1')
        
        conf_mat = overlap_matrix(loci_1, loci_2, type="IoU")

        assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')

        new_columns = ["{}|{}".format(c[0], c[1]) for c in assignment_pairs]

        if mnemons:
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
                mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]

            mnemon2_dict = {}
            for i in loci_2_mnemon:
                mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]

            #handle missing mnemonics
            for i in range(num_labels):
                if str(i) not in mnemon1_dict.keys():
                    mnemon1_dict[str(i)] = str(i)
                if str(i) not in mnemon2_dict.keys():
                    mnemon2_dict[str(i)] = str(i)
                

            for i in range(len(assignment_pairs)):
                assignment_pairs[i] = (mnemon1_dict[str(assignment_pairs[i][0])], mnemon2_dict[str(assignment_pairs[i][1])])
            print(assignment_pairs)

            loci_1, loci_2 = \
                connect_bipartite(loci_1, loci_2, assignment_pairs)
        else:
            loci_1, loci_2 = \
                connect_bipartite(loci_1, loci_2, assignment_pairs, mnemon=False)
        
        print('connected barpartite')

        num_labels = loci_1.shape[1]-3
        n_cols = math.floor(math.sqrt(num_labels))
        n_rows = math.ceil(num_labels / n_cols)

        fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
        label_being_plotted = 0
        
        for i in range(n_rows):
            for j in range(n_cols):
                c = list(loci_1.columns[3:])[label_being_plotted]
                cr, ar = granularity_vs_agreement_symmetric(loci_2.copy(), loci_1.copy(), query_label=c)

                axs[i,j].plot(cr, ar, c="black")
                axs[i,j].set_title(c, fontsize=7)
                label_being_plotted+=1
        plt.ylabel("Agreement")
        plt.xlabel("Genome Coverage")
        plt.tight_layout()
        plt.show()

    else:

        if mnemons:
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
                mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]

            mnemon2_dict = {}
            for i in loci_2_mnemon:
                mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]

            #handle missing mnemonics
            for i in range(num_labels):
                if str(i) not in mnemon1_dict.keys():
                    mnemon1_dict[str(i)] = str(i)
                if str(i) not in mnemon2_dict.keys():
                    mnemon2_dict[str(i)] = str(i)

            loci_1.columns = list(loci_1.columns[:3]) + [mnemon1_dict[c.replace("posterior","")] for c in loci_1.columns[3:]]
            loci_2.columns = list(loci_2.columns[:3]) + [mnemon2_dict[c.replace("posterior","")] for c in loci_2.columns[3:]]


        # print(loci_1)
        # print(loci_2)

        # for c in list(loci_1.columns[3:]):
        #     cr, ar, ovr_rec = granularity_vs_agreement_nonsymmetric(loci_1.copy(), loci_2.copy(), k=c)
        #     plot_progression(ar, cr, ovr_rec, c)

        # exit()

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
        print(replicate_1_dir.split("/")[2]+"_"+replicate_2_dir.split("/")[2])

        plt.savefig(replicate_1_dir+"/granularity.pdf")

def row_similarity(iou_matrix):
    n = iou_matrix.shape[0]
    sim_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            sim_matrix[i, j] = np.sqrt(np.max(iou_matrix.iloc[i,:] * iou_matrix.iloc[j,:]))
    np.fill_diagonal(sim_matrix, 1)
    return sim_matrix

def merging_gain(joint_overlap):
    coverage1 = {k:sum(joint_overlap.loc[k,:]) for k in joint_overlap.index}
    coverage2 = {k:sum(joint_overlap.loc[:,k]) for k in joint_overlap.columns}

    # nmi0 = NMI_from_matrix(joint_overlap, return_MI=False)
    gain = pd.DataFrame(
        np.zeros([joint_overlap.shape[0], joint_overlap.shape[0]]), 
        index=joint_overlap.index,
        columns=joint_overlap.index)
    
    for i in joint_overlap.index:
        for j in joint_overlap.index:
            joint_overlap_prime = joint_overlap.copy()

            if i == j:
                gain.loc[i, j] = 0

            else:
                MI0 = NMI_from_matrix(joint_overlap_prime, return_MI=True)
                H0 = -1 * np.sum([
                    np.sum(joint_overlap_prime.loc[p, :]) * np.log(np.sum(joint_overlap_prime.loc[p, :])) for p in joint_overlap_prime.index
                ])

                joint_prime_ij = joint_overlap.loc[i, :] + joint_overlap.loc[j, :]
                coverage_ij = np.sum(joint_prime_ij)

                joint_overlap_prime.loc[str(i) + " + " + str(j)] = joint_prime_ij
                joint_overlap_prime = joint_overlap_prime.drop([i, j], axis=0)

                MI1 = NMI_from_matrix(joint_overlap_prime, return_MI=True)
                H1 = -1 * np.sum([
                    np.sum(joint_overlap_prime.loc[p, :]) * np.log(np.sum(joint_overlap_prime.loc[p, :])) for p in joint_overlap_prime.index
                ])
                
                # nmi1 = NMI_from_matrix(joint_overlap_prime, return_MI=False)

                # gain.loc[i, j] = nmi1 - nmi0
                gain.loc[i, j] = 1 - ((MI1 - MI0) / (H1 - H0))


    # gain = np.array((gain - gain.min().min()) / (gain.max().max() - gain.min().min()))
    # print(gain)
    # gain = np.array(gain)
    # np.fill_diagonal(gain, np.min(gain) - 1)
    p = sns.heatmap(
        gain.astype(float), annot=True, fmt=".2f",
        linewidths=0.01,  cbar=True, annot_kws={"size": 8})

    p.tick_params(axis='x', rotation=90, labelsize=10)
    p.tick_params(axis='y', rotation=0, labelsize=10)

    plt.tight_layout()
    plt.show()
    # exit()

    return gain

def merge_clusters(joint, loci_1, loci_2, r1=True):
    m = 0

    #########################################################################################################
    if r1:
        num_labels = loci_1.shape[1]-3
        r1_gain = pd.DataFrame(merging_gain(joint), columns=joint.index, index=joint.index)

        max_value = r1_gain.max().max()

        # Find the indices of the maximum value in the similarity matrix
        max_index = np.where(r1_gain == max_value)

        # Get the row and column labels of the maximum value
        row_label = r1_gain.index[max_index[0][0]]
        col_label = r1_gain.columns[max_index[1][0]]

        print(f"The most similar pair is ({row_label}, {col_label}) with a similarity of {max_value}")

        linkage_1 = np.array([[max_index[0][0], max_index[1][0]]])

        merged_label_ID_1 = {}
        labels = loci_1.iloc[:, 3:].columns
        for i in range(num_labels):
            merged_label_ID_1[i] = labels[i]
        
        to_be_merged_1 = [
            merged_label_ID_1[int(linkage_1[m, 0])],
            merged_label_ID_1[int(linkage_1[m, 1])],
        ]

        merged_label_ID_1[num_labels + m] = str(
            merged_label_ID_1[int(linkage_1[m, 0])] + "+" + merged_label_ID_1[int(linkage_1[m, 1])]
        )

        loci_1[merged_label_ID_1[num_labels + m]] = \
            loci_1[to_be_merged_1[0]] + loci_1[to_be_merged_1[1]]
        loci_1 = loci_1.drop(to_be_merged_1, axis=1)

    #########################################################################################################
    else:
        num_labels = loci_2.shape[1]-3
        r2_gain = pd.DataFrame(merging_gain(joint.T), columns=joint.columns, index=joint.columns)

        max_value = r2_gain.max().max()

        # Find the indices of the maximum value in the similarity matrix
        max_index = np.where(r2_gain == max_value)

        # Get the row and column labels of the maximum value
        row_label = r2_gain.index[max_index[0][0]]
        col_label = r2_gain.columns[max_index[1][0]]

        print(f"The most similar pair is ({row_label}, {col_label}) with a similarity of {max_value}")

        linkage_2 = np.array([[max_index[0][0], max_index[1][0]]])

        merged_label_ID_2 = {}
        labels = loci_2.iloc[:, 3:].columns
        for i in range(num_labels):
            merged_label_ID_2[i] = labels[i]

        to_be_merged_2 = [
            merged_label_ID_2[int(linkage_2[m, 0])],
            merged_label_ID_2[int(linkage_2[m, 1])],
        ]

        merged_label_ID_2[num_labels + m] = str(
            merged_label_ID_2[int(linkage_2[m, 0])] + "+" + merged_label_ID_2[int(linkage_2[m, 1])]
        )

        loci_2[merged_label_ID_2[num_labels + m]] = \
            loci_2[to_be_merged_2[0]] + loci_2[to_be_merged_2[1]]
        loci_2 = loci_2.drop(to_be_merged_2, axis=1)

    return loci_1, loci_2

if __name__=="__main__":
    run_on_subset = True
    mnemons = True
    symmetric = False

    replicate_1_dir = "tests/cedar_runs/chmm/MCF7_R1/"
    replicate_2_dir = "tests/cedar_runs/chmm/MCF7_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons, symmetric)
    run(replicate_2_dir, replicate_1_dir, run_on_subset, mnemons, symmetric)

    replicate_1_dir = "tests/cedar_runs/segway/MCF7_R1/"
    replicate_2_dir = "tests/cedar_runs/segway/MCF7_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons, symmetric)
    run(replicate_2_dir, replicate_1_dir, run_on_subset, mnemons, symmetric)

    replicate_1_dir = "tests/cedar_runs/chmm/GM12878_R1/"
    replicate_2_dir = "tests/cedar_runs/chmm/GM12878_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons, symmetric)
    run(replicate_2_dir, replicate_1_dir, run_on_subset, mnemons, symmetric)

    replicate_1_dir = "tests/cedar_runs/segway/GM12878_R1/"
    replicate_2_dir = "tests/cedar_runs/segway/GM12878_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons, symmetric)
    run(replicate_2_dir, replicate_1_dir, run_on_subset, mnemons, symmetric)

 