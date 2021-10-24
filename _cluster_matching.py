import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

def find_max_posteri(loci, num_labels=10):
    loci_posteri =  loci.loc[:, ['posterior'+str(i) for i in range(int(num_labels))]]
    for lpc in loci_posteri.columns:
        loci_posteri[lpc] = loci_posteri[lpc].astype('float')

    max_posteri = loci_posteri.idxmax(axis=1)
    return max_posteri

def confusion_matrix(loci_1, loci_2, num_labels):
    num_labels = int(num_labels)
    max_1posteri = find_max_posteri(loci_1, num_labels=num_labels)
    max_2posteri = find_max_posteri(loci_2, num_labels=num_labels)

    Confus_Mat = pd.DataFrame(
        np.zeros((int(num_labels), int(num_labels))),
        columns=['posterior'+str(i)for i in range(num_labels)],
        index=['posterior'+str(i)for i in range(num_labels)]
    )
    for i in range(len(loci_1)):
        Confus_Mat.loc[max_1posteri[i], max_2posteri[i]] += 1
    return Confus_Mat

def Hungarian_algorithm(confusion_matrix, verbose=True):
    confusion_matrix = np.array(confusion_matrix)
    best_assignments = linear_sum_assignment(confusion_matrix, maximize=True)

    if verbose:
        print(
            'Sum of optimal assignment sets / Sum of confusion matrix = {}/{}'.format(
            str(confusion_matrix[best_assignments[0],best_assignments[1]].sum()), 
            str(confusion_matrix.sum())))

    assignment_pairs = [(i, best_assignments[1][i]) for i in range(len(best_assignments[0]))]
    return assignment_pairs

def connect_bipartite(loci_1, loci_2, assignment_matching):
    corrected_loci_2 = []
    corrected_loci_2.append(loci_2.chr)
    corrected_loci_2.append(loci_2.window_start)
    corrected_loci_2.append(loci_2.window_end)

    for i in range(len(assignment_matching)):
        corrected_loci_2.append(loci_2['posterior'+str(assignment_matching[i][1])])

    corrected_loci_2 = pd.concat(corrected_loci_2, axis=1)

    corrected_loci_2.columns = loci_1.columns
    return loci_1, corrected_loci_2
    

def read_length_dist_files(len_file_1, len_file_2):
    length_dist_1 = pd.read_csv(len_file_1, sep='\t')
    length_dist_2 = pd.read_csv(len_file_2, sep='\t')

    return length_dist_1, length_dist_2

def gmtk_params_files(gmtk_file_1, gmtk_file_2):
    gmtk_params_1 = pd.read_csv(gmtk_file_1)
    gmtk_params_1 = gmtk_params_1.drop("Unnamed: 0", axis=1)
    gmtk_params_2 = pd.read_csv(gmtk_file_2)
    gmtk_params_2 = gmtk_params_2.drop("Unnamed: 0", axis=1)

    return gmtk_params_1, gmtk_params_2

def create_clustering_data(length_dist_1, length_dist_2, gmtk_params_1, gmtk_params_2):
    length_dist_1 = length_dist_1.drop(0, axis=0)
    length_dist_1 = length_dist_1.loc[:, ('mean.len', 'stdev.len')]

    length_dist_2 = length_dist_2.drop(0, axis=0)
    length_dist_2 = length_dist_2.loc[:, ('mean.len', 'stdev.len')]

    rep1_data = pd.concat([gmtk_params_1, length_dist_1], axis=1)
    rep2_data = pd.concat([gmtk_params_2, length_dist_2], axis=1)

    return rep1_data, rep2_data

def cluster(clustering_data):
    pass

def update_labels(loci_1, loci_2):
    """
    two ways:
    either connect using the confusion matrix with updated (clustered) labels
    OR pair clusters based on min euclidean distance
    """
    pass

