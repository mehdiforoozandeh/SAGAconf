from .run import *
from scipy.stats import norm

"""
This module contains function to investigate how "inexact" segment transition
affects reproducibility between two sets of annotations. 
"""

def BoundaryAnalysis_Hist(loci_1, loci_2, match_definition="BM"):
    '''
    input locis must NOT be matched -- optionally with updated mnemonics
    '''
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

    # now define a boundary
    boundary_distances = {}
    ratio_inexact_matched = {}
    for k in list(loci_1.columns[3:]):
        boundary_distances[k] = []
        ratio_inexact_matched[k] = [0,0] #matched within 100 segments, not matched within 100 segments

    for i in range(len(MAP1)-1):
        if MAP1[i] != MAP1[i+1]: 
            #this means transition
            """
            now check how far should we go from point i in the other replicate to see a good match.
            """
            
            if MAP2[i] not in per_label_matches[MAP1[i]]:
                matched = False
                dist = 0
                while matched == False and dist<100:
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
                    ratio_inexact_matched[MAP1[i]][0] += 1
                else:
                    ratio_inexact_matched[MAP1[i]][1] += 1
                    
    # plot histograms

    num_labels = len(boundary_distances.keys())
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
    label_being_plotted = 0

    for i in range(n_rows):
        for j in range(n_cols):
            label = loci_1.columns[3:][label_being_plotted]
            axs[i,j].hist(boundary_distances[label], bins=100)
            label_being_plotted += 1
            axs[i,j].set_title(label, fontsize=8)
            # axs[i,j].set_title(label+" | {:.2f}".format(
            #     float(ratio_inexact_matched[label][0] / (ratio_inexact_matched[label][0] + ratio_inexact_matched[label][1]))
            #     ), fontsize=8)

    fig.text(0.5, 0.005, 'Distance to first good match' , ha='center')
    fig.text(0.005, 0.5, 'Number of inexact match', va='center', rotation='vertical')
    plt.tight_layout()
    plt.show()

def BoundaryAnalysis_cumulative(loci_1, loci_2, match_definition="BM", max_distance=50):
    '''
    input locis must NOT be matched -- optionally with updated mnemonics
    '''
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
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
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

            print(list(matched_hist[1][:-1]))
            print(list(matched_hist[1][:-1]*resolution))
            print(cdf)
            axs[i,j].plot(list(matched_hist[1][:-1]*resolution), list(cdf))
            axs[i,j].set_xticks(np.arange(0, (max_distance+1)*resolution, step=5*resolution))
            axs[i,j].tick_params(axis='both', labelsize=7)
            axs[i,j].tick_params(axis='x', rotation=90)
            axs[i,j].set_yticks(np.arange(0, 1.1, step=0.2))
            
            axs[i,j].set_title(k, fontsize=8)

            label_being_plotted += 1
            
    fig.text(0.5, 0.005, 'Distance from boundary' , ha='center')
    fig.text(0.005, 0.5, 'Ratio matched boundaries', va='center', rotation='vertical')
    plt.tight_layout()
    plt.show()

def run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons):
    print("loading and intersecting")
    loci_1, loci_2 = intersect_parsed_posteriors(
        replicate_1_dir, 
        replicate_2_dir)

    if run_on_subset:
        # loci_1 = loci_1.iloc[[i for i in range(0, len(loci_1), 100)], :].reset_index(drop=True)
        # loci_2 = loci_2.iloc[[i for i in range(0, len(loci_2), 100)], :].reset_index(drop=True)

        loci_1 = loci_1.loc[loci_1["chr"]=="chr7"].reset_index(drop=True)
        loci_2 = loci_2.loc[loci_2["chr"]=="chr7"].reset_index(drop=True)
        print(loci_1)
        print(loci_2)
    
    print("the shapes of the input matrices are: {}, {}".format(str(loci_1.shape), str(loci_2.shape)))

    num_labels = loci_1.shape[1]-3
    loci_1.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    loci_2.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]

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

    # BoundaryAnalysis(loci_1, loci_2)
    BoundaryAnalysis_cumulative(loci_1, loci_2, match_definition="BM", max_distance=100)

if __name__=="__main__":
    run_on_subset = True
    mnemons = True

    # replicate_1_dir = "cedar_runs/chmm/gm12878_r1/parsed_posterior.csv"
    # replicate_2_dir = "cedar_runs/chmm/gm12878_r2/parsed_posterior.csv"
    # run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    # replicate_1_dir = "cedar_runs/chmm/cd14_r1/parsed_posterior.csv"
    # replicate_2_dir = "cedar_runs/chmm/cd14_r2/parsed_posterior.csv"
    # run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    # replicate_1_dir = "cedar_runs/chmm/hela_r1/parsed_posterior.csv"
    # replicate_2_dir = "cedar_runs/chmm/hela_r2/parsed_posterior.csv"
    # run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    # replicate_1_dir = "cedar_runs/chmm/k562_r1/parsed_posterior.csv"
    # replicate_2_dir = "cedar_runs/chmm/k562_r2/parsed_posterior.csv"
    # run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    # replicate_1_dir = "cedar_runs/chmm/mcf7_r1/parsed_posterior.csv"
    # replicate_2_dir = "cedar_runs/chmm/mcf7_r2/parsed_posterior.csv"
    # run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    # ###
    ###
    ###

    replicate_1_dir = "cedar_runs/segway/gm12878_r1/parsed_posterior.csv"
    replicate_2_dir = "cedar_runs/segway/gm12878_r2/parsed_posterior.csv"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    replicate_1_dir = "cedar_runs/segway/cd14_r1/parsed_posterior.csv"
    replicate_2_dir = "cedar_runs/segway/cd14_r2/parsed_posterior.csv"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    replicate_1_dir = "cedar_runs/segway/hela_r1/parsed_posterior.csv"
    replicate_2_dir = "cedar_runs/segway/hela_r2/parsed_posterior.csv"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    replicate_1_dir = "cedar_runs/segway/k562_r1/parsed_posterior.csv"
    replicate_2_dir = "cedar_runs/segway/k562_r2/parsed_posterior.csv"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)

    replicate_1_dir = "cedar_runs/segway/mcf7_r1/parsed_posterior.csv"
    replicate_2_dir = "cedar_runs/segway/mcf7_r2/parsed_posterior.csv"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons)