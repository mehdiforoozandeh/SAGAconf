from run import *
import seaborn as sns
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
from length_dist_analysis import *
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial import distance
from granul import *

def get_ecdf_for_label(loci, label, max=None):
    MAP = loci.iloc[:,3:].idxmax(axis=1)
    resolution = int(loci.iloc[0, 2] - loci.iloc[0, 1])

    label_segments = []

    midlabel = False

    for i in range(len(MAP)):
        if MAP[i] == label:
            if midlabel == False:
                label_segments.append(1)
                midlabel = True
            else:
                label_segments[-1] += 1
            
            if i+1 < len(MAP):
                if MAP[i+1] != label:
                    midlabel = False

    lendist = np.array(label_segments) * resolution
    if max==None:
        q = np.quantile(lendist, 0.99)
    else:
        q = max

    lendist = pd.Series(lendist, name="length_dist")
    lendist = lendist.loc[lendist < q]
    try:
        return lendist
    except:
        pass

def length_vs_boundary(loci_1, loci_2, outdir, match_definition="BM", max_distance=50):
    """
    Here i'm gonna merge ECDF of length dist and boundary thing.
    """

    resolution = int(loci_1.iloc[0, 2] - loci_1.iloc[0, 1])
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
    n_rows = num_labels
    n_cols = 1

    fig, axs = plt.subplots(n_rows, n_cols, sharex="col", sharey="col", figsize=[5, 30])
    label_being_plotted = 0
    
    for i in range(n_rows):
        # for j in range(n_cols):

        k = loci_1.columns[3:][label_being_plotted]

        matched_hist = np.histogram(boundary_distances[k], bins=max_distance, range=(0, max_distance)) #[0]is the bin size [1] is the bin edge

        cdf = []
        for jb in range(len(matched_hist[0])):
            if len(cdf) == 0:
                cdf.append(matched_hist[0][jb])
            else:
                cdf.append((cdf[-1] + matched_hist[0][jb]))

        cdf = np.array(cdf) / (np.sum(matched_hist[0]) + count_unmatched_boundaries[k])
        
        
        axs[i].plot(list(matched_hist[1][:-1]*resolution), list(cdf), color="black", label="distance from boundary")
        
        lendist = get_ecdf_for_label(loci_1, k, max=(max_distance+1)*resolution)
        if len(lendist)>0:
            sns.ecdfplot(lendist, ax=axs[i], color="red")

        axs[i].set_xticks(np.arange(0, (max_distance+1)*resolution, step=5*resolution))
        axs[i].tick_params(axis='both', labelsize=7)
        axs[i].tick_params(axis='x', rotation=90)
        axs[i].set_yticks(np.arange(0, 1.1, step=0.2))
        axs[i].set_xlabel('bp', fontsize=7)
        axs[i].set_ylabel('Proportion',fontsize=7)
        
        
        axs[i].set_title(k, fontsize=7)

        label_being_plotted += 1
            
    plt.tight_layout()
    plt.savefig(outdir)

def read_sigdist_file(file):
    return pd.read_csv(file, sep="\t")

def compare_signal_dist(loci_1, loci_2, sigdist1_file, sigdist2_file, k):
    """
    read signal_dist files
    for k in loci_1 labels:
        get granularity_vs_agreement_nonsymmetric(loci_1.copy(), loci_2.copy(), k=c)
        get sigdist1.k and sigdist2.k

        plot side by side

    """
    sigdist1 = read_sigdist_file(sigdist1_file)
    sigdist2 = read_sigdist_file(sigdist2_file)

    if np.min(np.array(sigdist1['label'].astype("int"))) == 1:
        for i in range(sigdist1.shape[0]):
            sigdist1["label"][i] = str(int(sigdist1["label"][i])-1)

    if np.min(np.array(sigdist2['label'].astype("int"))) == 1:
        for i in range(sigdist2.shape[0]):
            sigdist2["label"][i] = str(int(sigdist2["label"][i])-1)

    # update sigdist label names

    for i in range(sigdist1.shape[0]):
        sigdist1["label"][i] = [x for x in loci_1.columns[3:] if str(sigdist1["label"][i])== x.split("_")[0]][0]

    for i in range(sigdist2.shape[0]):
        sigdist2["label"][i] = [x for x in loci_2.columns[3:] if str(sigdist2["label"][i])== x.split("_")[0]][0]

    sigdist1 = sigdist1.drop(["sd"], axis=1)
    sigdist2 = sigdist2.drop(["sd"], axis=1)

    # print(sigdist1)
    # print(sigdist2)

    ## rescale each assay to 0-1
    list_of_assays = np.unique(np.array(sigdist1['trackname']))
    for a in list_of_assays:
        min_max_1 = (
            np.min(np.array(sigdist1.loc[sigdist1["trackname"]==a, ["mean"]])),
            np.max(np.array(sigdist1.loc[sigdist1["trackname"]==a, ["mean"]]))
            )
        
        min_max_2 = (
            np.min(np.array(sigdist2.loc[sigdist1["trackname"]==a, ["mean"]])),
            np.max(np.array(sigdist2.loc[sigdist1["trackname"]==a, ["mean"]]))
            )

        for i in range(len(sigdist1)):
            if sigdist1["trackname"][i] == a:
                sigdist1["mean"][i] = float(
                    (float(sigdist1["mean"][i]) - min_max_1[0]) / (min_max_1[1] - min_max_1[0])
                )

                sigdist2["mean"][i] = float(
                    (float(sigdist2["mean"][i]) - min_max_2[0]) / (min_max_2[1] - min_max_2[0])
                )

    ###
    # print(sigdist1)
    # print(sigdist2)
    confmat = enrichment_of_overlap_matrix(
            loci_1, loci_2, OE_transform=True)
    
    sorted_k_vector = confmat.loc[k,:].sort_values(ascending=False)
    print(sorted_k_vector)

    """
    label k in loci_1
    label sorted_k_vector[0] in loci_2alan 

    """
    v1 = sigdist1.loc[(sigdist1["label"] == k), ["trackname","mean"]]

    v1.index = v1["trackname"].values
    v1 = v1.drop(["trackname"], axis=1)

    v2 = sigdist2.loc[(sigdist2["label"] == sorted_k_vector.index[0]), ["trackname", "mean"]]

    v2.index = v2["trackname"].values
    v2 = v2.drop(["trackname"], axis=1)

    cosinesim = cosine_similarity(np.reshape(np.array(v1), (1, -1)), np.reshape(np.array(v2), (1, -1)))
    euc_distance = distance.euclidean(np.array(v1), np.array(v2))

    compar = pd.concat([v1, v2], axis=1)
    compar.columns = ["R1: "+k, "R2: "+sorted_k_vector.index[0]]
    compar = compar.T

    return compar, cosinesim, euc_distance

def plot_sigdist_compar(loci_1, loci_2, sigdistfile1, sigdistfile2, outdir):
    num_labels = len(loci_1.columns)-3

    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=False, sharey=False, figsize=[25, 16])
    label_being_plotted = 0
    
    for i in range(n_rows):
        for j in range(n_cols):

            k = loci_1.columns[3:][label_being_plotted]
            compar, cosine, eucld = compare_signal_dist(loci_1, loci_2, sigdistfile1, sigdistfile2, k)

            sns.heatmap(compar, ax=axs[i, j], cmap="crest", cbar=False)

            axs[i, j].tick_params(axis='both', labelsize=7)
            axs[i, j].tick_params(axis='x', rotation=90)
            axs[i, j].set_xlabel('Tracknames', fontsize=7)
            axs[i,j].set_title(k +" | cos: {:.2f} | euc: {:.2f}".format(float(cosine[0]), eucld), fontsize=7)
            label_being_plotted += 1

    plt.tight_layout()
    plt.savefig(outdir)

def sigdist_vs_overlap(sigdist1_file, sigdist2_file, loci_1, loci_2, outdir):
    sigdist1 = read_sigdist_file(sigdist1_file)
    sigdist2 = read_sigdist_file(sigdist2_file)

    if np.min(np.array(sigdist1['label'].astype("int"))) == 1:
        for i in range(sigdist1.shape[0]):
            sigdist1["label"][i] = str(int(sigdist1["label"][i])-1)

    if np.min(np.array(sigdist2['label'].astype("int"))) == 1:
        for i in range(sigdist2.shape[0]):
            sigdist2["label"][i] = str(int(sigdist2["label"][i])-1)


    #
    """
    add an option to match clusters based on hungarian algorithm
    """
    # confmat = enrichment_of_overlap_matrix(
    #     loci_1, loci_2, OE_transform=True)

    # assignmet_pairs = Hungarian_algorithm(confmat)

    # print(assignmet_pairs)
    # loci_1, loci_2 = connect_bipartite(loci_1, loci_2, assignmet_pairs)
    

    # update sigdist label names
    for i in range(sigdist1.shape[0]):
        sigdist1["label"][i] = [x for x in loci_1.columns[3:] if str(sigdist1["label"][i])== x.split("_")[0]][0]

    for i in range(sigdist2.shape[0]):
        sigdist2["label"][i] = [x for x in loci_2.columns[3:] if str(sigdist2["label"][i])== x.split("_")[0]][0]

    sigdist1 = sigdist1.drop(["sd"], axis=1)
    sigdist2 = sigdist2.drop(["sd"], axis=1)

    # print(sigdist1)
    # print(sigdist2)

    ## rescale each assay to 0-1
    list_of_assays = np.unique(np.array(sigdist1['trackname']))
    for a in list_of_assays:
        min_max_1 = (
            np.min(np.array(sigdist1.loc[sigdist1["trackname"]==a, ["mean"]])),
            np.max(np.array(sigdist1.loc[sigdist1["trackname"]==a, ["mean"]]))
            )
        
        min_max_2 = (
            np.min(np.array(sigdist2.loc[sigdist1["trackname"]==a, ["mean"]])),
            np.max(np.array(sigdist2.loc[sigdist1["trackname"]==a, ["mean"]]))
            )

        for i in range(len(sigdist1)):
            if sigdist1["trackname"][i] == a:
                sigdist1["mean"][i] = float(
                    (float(sigdist1["mean"][i]) - min_max_1[0]) / (min_max_1[1] - min_max_1[0])
                )

                sigdist2["mean"][i] = float(
                    (float(sigdist2["mean"][i]) - min_max_2[0]) / (min_max_2[1] - min_max_2[0])
                )

    ###

    confmat = enrichment_of_overlap_matrix(
        loci_1, loci_2, OE_transform=True)

    eucl = pd.DataFrame(
        np.zeros((len(loci_1.columns)-3, len(loci_2.columns)-3)),
        columns=loci_2.columns[3:], 
        index=loci_1.columns[3:])

    cossim = pd.DataFrame(
        np.zeros((len(loci_1.columns)-3, len(loci_2.columns)-3)),
        columns=loci_2.columns[3:], 
        index=loci_1.columns[3:])

    for i in range(len(loci_1.columns[3:])):
        for j in range(len(loci_2.columns[3:])):

            v1 = sigdist1.loc[(sigdist1["label"] == loci_1.columns[3:][i]), ["trackname", "mean"]]
            v1 = v1.sort_values(by=["trackname"])

            v2 = sigdist2.loc[(sigdist2["label"] == loci_2.columns[3:][j]), ["trackname", "mean"]]
            v2 = v2.sort_values(by=["trackname"])

            eucl.loc[loci_1.columns[3:][i], loci_2.columns[3:][j]] = distance.euclidean(
                np.array(v1["mean"]),
                np.array(v2["mean"]))

            cossim.loc[loci_1.columns[3:][i], loci_2.columns[3:][j]] = cosine_similarity(
                np.reshape(np.array(v1["mean"]) ,(1,-1)), 
                np.reshape(np.array(v2["mean"]) ,(1,-1))
            )

    # rescale euc
    # eucl = pd.DataFrame(
    #     (np.array(eucl) - np.min(np.array(eucl))) / (np.max(np.array(eucl)) - np.min(np.array(eucl))),
    #     columns=loci_2.columns[3:], 
    #     index=loci_1.columns[3:])

    fig, axs = plt.subplots(1, 3, sharex=False, sharey=False, figsize=[30, 10])
    
    sns.heatmap(confmat, annot=True, fmt=".2f", linewidths=0.01,  cbar=False, ax=axs[0])
    axs[0].set_title("Enrichment of Overlap")
    axs[0].set_xlabel("R1")
    axs[0].set_ylabel("R2")

    cmap_r = sns.cm.rocket_r
    sns.heatmap(eucl, annot=True, fmt=".2f", linewidths=0.01,  cbar=False, ax=axs[1], cmap = cmap_r)
    axs[1].set_title("Euclidean distance of signal distribution")
    axs[1].set_xlabel("R1")
    axs[1].set_ylabel("R2")

    sns.heatmap(cossim, annot=True, fmt=".2f", linewidths=0.01,  cbar=False, ax=axs[2])
    axs[2].set_title("Cosine similarity of signal distribution")
    axs[2].set_xlabel("R1")
    axs[2].set_ylabel("R2")

    plt.savefig(outdir)

def run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons):
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
    
    print('generating confmat 1')
    
    # conf_mat = confusion_matrix(
    #     loci_1, loci_2, num_labels, 
    #     OE_transform=True, symmetric=False)

    # assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')

    # new_columns = ["{}|{}".format(c[0], c[1]) for c in assignment_pairs]

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

        print(loci_1)
        print(loci_2)

        # R1 vs R2
        sigdist_vs_overlap(
            replicate_1_dir+"/sigdist/signal_distribution.tab",
            replicate_2_dir+"/sigdist/signal_distribution.tab", 
            loci_1, loci_2, 
            outdir=replicate_1_dir+"/overlap_eucd_cossim.pdf")
        
        plot_sigdist_compar(
                    loci_1, loci_2, 
                    replicate_1_dir+"/sigdist/signal_distribution.tab",
                    replicate_2_dir+"/sigdist/signal_distribution.tab", 
                    outdir=replicate_1_dir+"/sigdistcompar.pdf")
        
        length_vs_boundary(
            loci_1, loci_2, match_definition="BM", max_distance=100, 
            outdir=replicate_1_dir+"/len_bound.pdf")

        # R2 vs R1
        sigdist_vs_overlap(
            replicate_2_dir+"/sigdist/signal_distribution.tab",
            replicate_1_dir+"/sigdist/signal_distribution.tab", 
            loci_2, loci_1, 
            outdir=replicate_2_dir+"/overlap_eucd_cossim.pdf")
        
        plot_sigdist_compar(
                    loci_2, loci_1, 
                    replicate_2_dir+"/sigdist/signal_distribution.tab",
                    replicate_1_dir+"/sigdist/signal_distribution.tab", 
                    outdir=replicate_2_dir+"/sigdistcompar.pdf")
        
        length_vs_boundary(
            loci_2, loci_1, match_definition="BM", max_distance=50, 
            outdir=replicate_2_dir+"/len_bound.pdf")
        
        # exit()

def read_bedGraph_in_bp(bgfile, subset="chr21"):
    arr = pd.read_csv(bgfile, sep="\t", header=None)
    arr.columns = ["chr", "start", "end", "signal"]
    arr = arr.loc[arr["chr"]==subset, :]
    arr = np.array(arr)

    new_arr = []
    for i in range(len(arr)):
        _chr = arr[i][0]
        start = arr[i][1]
        end = arr[i][2]
        val = arr[i][3]
        
        while start < end:
            new_arr.append([str(_chr), int(start), int(start+1), float(val)])
            start +=1
    
    arr = np.array(new_arr)
    arr = pd.DataFrame(arr, columns=["chr", "start", "end", "signal"])
    arr['chr'] = arr["chr"].astype(str)
    arr['start'] = arr["start"].astype(int)
    arr['end'] = arr["end"].astype(int)
    arr['signal'] = arr["signal"].astype(float)
    return arr

def compare_tracks_heatmap(df1, df2):
    intersect = pd.merge(
        df1, 
        df2, 
        how='inner', 
        on=['chr', 'start', 'end'])

    df1 = [intersect.chr, intersect.start, intersect.end]
    df2 = [intersect.chr, intersect.start, intersect.end]

    for c in intersect.columns:
        if c[-1] == 'x':
            df1.append(intersect[c])
        elif c[-1] == 'y':
            df2.append(intersect[c])

    df1 = pd.concat(df1, axis=1)
    df1.columns = [c.replace("_x", "") for c in df1.columns]
    
    df2 = pd.concat(df2, axis=1)
    df2.columns = [c.replace("_y", "") for c in df2.columns]

    heatmap, xedges, yedges = np.histogram2d(df1.iloc[:,3], df2.iloc[:,3], bins=50)  
    heatmap = heatmap / len(df1)

    for i in range(heatmap.shape[0]):
        if np.sum(heatmap[i,:]) != 0:
            heatmap[i,:] = heatmap[i,:] / np.sum(heatmap[i,:])
    
    return heatmap, xedges, yedges

def compare_tracks_cdf(df1, df2):
    intersect = pd.merge(
        df1, 
        df2, 
        how='inner', 
        on=['chr', 'start', 'end'])

    df1 = [intersect.chr, intersect.start, intersect.end]
    df2 = [intersect.chr, intersect.start, intersect.end]

    for c in intersect.columns:
        if c[-1] == 'x':
            df1.append(intersect[c])
        elif c[-1] == 'y':
            df2.append(intersect[c])

    df1 = pd.concat(df1, axis=1)
    df1.columns = [c.replace("_x", "") for c in df1.columns]
    
    df2 = pd.concat(df2, axis=1)
    df2.columns = [c.replace("_y", "") for c in df2.columns]

    df1 = df1.iloc[:, 3]
    df2 = df2.iloc[:, 3]

    sorted_data_1 = np.sort(df1)
    cumulative_1 = np.cumsum(sorted_data_1)
    normalized_1 = cumulative_1 / cumulative_1[-1]

    sorted_data_2 = np.sort(df2)
    cumulative_2 = np.cumsum(sorted_data_2)
    normalized_2 = cumulative_2 / cumulative_2[-1]

    return sorted_data_1, normalized_1,  sorted_data_2, normalized_2

def compare_ct_alltracks(ctdir, savedir):
    list_of_tracks = [name for name in os.listdir(ctdir) if os.path.isdir(os.path.join(ctdir, name))]
    
    num_ct = len(list_of_tracks)
    n_cols = math.floor(math.sqrt(num_ct))
    n_rows = math.ceil(num_ct / n_cols)

    ##################################################################################################
    fig, axs = plt.subplots(n_rows, n_cols, sharex=False, sharey=False, figsize=[n_cols*4, n_rows*4])
    label_being_plotted = 0
    
    for i in range(n_rows):
        for j in range(n_cols):
            if label_being_plotted < num_ct:
                k = list_of_tracks[label_being_plotted]
                bgfiles = ["{}/{}/{}".format(ctdir, k, x) for x in os.listdir(ctdir+"/"+k) if ".bedGraph" in x]

                # Create a density heatmap with the pcolormesh() function
                df1 = read_bedGraph_in_bp(bgfiles[0], subset="chr21")
                df2 = read_bedGraph_in_bp(bgfiles[1], subset="chr21")
                heatmap, xedges, yedges = compare_tracks_heatmap(df1, df2)

                im = axs[i, j].pcolormesh(xedges, yedges, heatmap.T, cmap=plt.cm.get_cmap('inferno'))
                axs[i, j].set_title(k)

                # Add a colorbar to the heatmap
                fig.colorbar(im, ax=axs[i, j])
                
                label_being_plotted += 1
            
    plt.tight_layout()
    plt.savefig(savedir+"/2dhistogram.pdf", format='pdf')
    plt.savefig(savedir+"/2dhistogram.svg", format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')

    ##################################################################################################
    fig, axs = plt.subplots(n_rows, n_cols, sharex=False, sharey=False, figsize=[n_cols*3, n_rows*2])
    label_being_plotted = 0
    
    for i in range(n_rows):
        for j in range(n_cols):
            if label_being_plotted < num_ct:
                k = list_of_tracks[label_being_plotted]
                bgfiles = ["{}/{}/{}".format(ctdir, k, x) for x in os.listdir(ctdir+"/"+k) if ".bedGraph" in x]

                df1 = read_bedGraph_in_bp(bgfiles[0], subset="chr21")
                df2 = read_bedGraph_in_bp(bgfiles[1], subset="chr21")
                sorted_data_1, normalized_1,  sorted_data_2, normalized_2 = compare_tracks_cdf(df1, df2)

                axs[i, j].plot(sorted_data_1, normalized_1)
                axs[i, j].plot(sorted_data_2, normalized_2)
                axs[i, j].set_title(k)

                label_being_plotted += 1
            
    plt.tight_layout()
    plt.savefig(savedir+"/cdf.pdf", format='pdf')
    plt.savefig(savedir+"/cdf.svg", format='svg')

    sns.reset_orig
    plt.close("all")
    plt.style.use('default')

if __name__=="__main__":

    compare_ct_alltracks("files/CD14-positive_monocyte/", "files/CD14-positive_monocyte/")
    compare_ct_alltracks("files/GM12878/", "files/GM12878/")
    compare_ct_alltracks("files/HeLa-S3/", "files/HeLa-S3/")
    compare_ct_alltracks("files/K562/", "files/K562/")
    compare_ct_alltracks("files/MCF-7/", "files/MCF-7/")

    exit()
    run_on_subset = True
    mnemons = True

    replicate_1_dir = "tests/cedar_runs/chmm/MCF7_R1/"
    replicate_2_dir = "tests/cedar_runs/chmm/MCF7_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    replicate_1_dir = "tests/cedar_runs/segway/MCF7_R1/"
    replicate_2_dir = "tests/cedar_runs/segway/MCF7_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    replicate_1_dir = "tests/cedar_runs/chmm/GM12878_R1/"
    replicate_2_dir = "tests/cedar_runs/chmm/GM12878_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    replicate_1_dir = "tests/cedar_runs/segway/GM12878_R1/"
    replicate_2_dir = "tests/cedar_runs/segway/GM12878_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)