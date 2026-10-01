import math
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import cohen_kappa_score


def overall_overlap_ratio(loci_1, loci_2, w=0):
    joint = joint_overlap_prob(loci_1, loci_2, w=w, symmetric=False) * len(loci_1)
    sum_matched = np.sum([np.max(joint.loc[i, :]) for i in joint.index])
    return sum_matched / len(loci_1)


def overlap_matrix(loci_1, loci_2, type="IoU"):
    if type == "IoU":
        return IoU_overlap(loci_1, loci_2)
    elif type == "ck":
        return Cohens_Kappa_matrix(loci_1, loci_2)
    elif type == "enr":
        return enrichment_of_overlap_matrix(loci_1, loci_2, OE_transform=True)
    elif type == "conditional":
        return enrichment_of_overlap_matrix(loci_1, loci_2, OE_transform=False)
    elif type == "hard_joint":
        return joint_overlap_prob(loci_1, loci_2, w=0)
    elif type == "soft_joint":
        return soft_joint_prob(loci_1, loci_2)


def soft_joint_prob(loci_1, loci_2):
    """
    this kind of overlap, takes the posterior prob into account. 
    the posterior prob should not be in logit form.
        if 0<p<1:
            pass
        else:
            p = sigmoid(p)
    """

    num_labels = len(loci_1.columns) -3
    joint = pd.DataFrame(np.zeros((num_labels, num_labels)), columns=loci_2.columns[3:], index=loci_1.columns[3:])

    for k in loci_1.columns[3:]:
        for j in loci_2.columns[3:]:
            soft_overlap = sum(np.array(loci_1.loc[:, k]) * np.array(loci_2.loc[:, j])) / len(loci_1)
            joint.loc[k, j] = soft_overlap
    
    return joint


def soft_coverage(loci_1, loci_2):
    coverage1 = {k:sum(loci_1.loc[:,k])/len(loci_1) for k in loci_1.columns[3:]}
    coverage2 = {k:sum(loci_2.loc[:,k])/len(loci_2) for k in loci_2.columns[3:]}
    return coverage1, coverage2


def joint_overlap_prob(loci_1, loci_2, w=0, symmetric=True):
    num_labels1 = len(loci_1.columns) -3
    num_labels2 = len(loci_2.columns) -3
    resolution = loci_1["end"][0] - loci_1["start"][0] 
    w = math.ceil(w/resolution)

    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

    observed_overlap = {} 

    for k in loci_1.columns[3:]:
        for j in loci_2.columns[3:]:
            observed_overlap[str(k + "|" + j)] = 0

    oo_mat = pd.DataFrame(np.zeros((num_labels1, num_labels2)), columns=loci_2.columns[3:], index=loci_1.columns[3:])

    MAP1 = list(MAP1)
    MAP2 = list(MAP2)
                
    if w == 0:
        for i in range(len(MAP1)):
            k = MAP1[i]
            j = MAP2[i]
            observed_overlap[str(k + "|" + j)] += 1

    elif symmetric==False:
    ############################### NON-SYMMETRIC ##################################
        for i in range(len(MAP1)):
            k = MAP1[i]
            i_neighbors = MAP2[max(0, i-w) : min(i+w, len(MAP2)-1)]
            for j in set(i_neighbors):
                observed_overlap[str(k + "|" + j)] += 1
        
    else:
    ################################# SYMMETRIC ####################################
        for i in range(len(MAP1)):
            R1_window =  MAP1[max(0, i-w) : min(i+w, len(MAP1)-1)]
            R2_window =  MAP2[max(0, i-w) : min(i+w, len(MAP2)-1)]
            for k in R1_window:
                for j in R2_window:
                    observed_overlap[str(k + "|" + j)] += float(1/(len(R1_window)*len(R2_window)))

    for p in observed_overlap.keys():
        oo_mat.loc[p.split("|")[0], p.split("|")[1]] = observed_overlap[p]
    
    return oo_mat / len(loci_1)


def IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False, overlap_coeff=False):
    num_labels1 = len(loci_1.columns) - 3
    num_labels2 = len(loci_2.columns) - 3
    
    IoU = pd.DataFrame(np.zeros((num_labels1, num_labels2)), columns=loci_2.columns[3:], index=loci_1.columns[3:])

    if soft and w == 0:
        joint = soft_joint_prob(loci_1, loci_2)
        coverage1, coverage2 = soft_coverage(loci_1, loci_2)

    else:
        joint = joint_overlap_prob(loci_1, loci_2, w=w, symmetric=symmetric)

        MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
        MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

        coverage1 = {k:len(MAP1.loc[MAP1 == k])/len(loci_1) for k in loci_1.columns[3:]}
        coverage2 = {k:len(MAP2.loc[MAP2 == k])/len(loci_2) for k in loci_2.columns[3:]}

        # coverage1 = {k: joint.loc[k,:].sum() for k in joint.index}
        # coverage2 = {k: joint.loc[:,k].sum() for k in joint.columns}

    for A in loci_1.columns[3:]:
        for B in loci_2.columns[3:]:
            if overlap_coeff:
                IoU.loc[A, B] = (joint.loc[A, B])/np.min([coverage1[A], coverage2[B]])

            else:
                IoU.loc[A, B] = (joint.loc[A, B])/(coverage1[A] + coverage2[B] - ((joint.loc[A, B])))

    return IoU


def Cohens_Kappa_matrix(loci_1, loci_2):
    num_labels = len(loci_1.columns) -3
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

    coverage1 = {k:len(MAP1.loc[MAP1 == k])/len(loci_1) for k in loci_1.columns[3:]}
    coverage2 = {k:len(MAP2.loc[MAP2 == k])/len(loci_2) for k in loci_2.columns[3:]}

    ck_mat = pd.DataFrame(
        np.zeros((num_labels, num_labels)), 
        columns=loci_2.columns[3:], index=loci_1.columns[3:])

    for k in loci_1.columns[3:]:
        for l in loci_2.columns[3:]:
            o = 0
            for i in range(len(MAP1)):
                if MAP1[i] == k and MAP2[i] == l:
                    o+=1
            
            ck_mat.loc[k, l] = cohen_kappa_score( 
                list(MAP1 == k),
                list(MAP2 == l)
            )

    return ck_mat


def enrichment_of_overlap_matrix(loci_1, loci_2, OE_transform=True):
    num_labels = len(loci_1.columns) -3
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

    coverage1 = {k: len(MAP1.loc[MAP1 == k]) / len(MAP1) for k in loci_1.columns[3:]}
    coverage2 = {k: len(MAP2.loc[MAP2 == k]) / len(MAP2) for k in loci_2.columns[3:]}

    observed_overlap = {} # np.zeros((num_labels, num_labels))
    expected_overlap = {} # np.zeros((num_labels, num_labels))

    for k in loci_1.columns[3:]:
        for j in loci_2.columns[3:]:
            observed_overlap[str(k + "|" + j)] = 0
            expected_overlap[str(k + "|" + j)] = coverage1[k] * coverage2[j] * len(MAP1)

    for i in range(MAP1.shape[0]):
        k = MAP1[i]
        j = MAP2[i]

        observed_overlap[str(k + "|" + j)] += 1

    oo_mat = pd.DataFrame(np.zeros((num_labels, num_labels)), columns=loci_2.columns[3:], index=loci_1.columns[3:])
    eo_mat = pd.DataFrame(np.zeros((num_labels, num_labels)), columns=loci_2.columns[3:], index=loci_1.columns[3:])

    for p in observed_overlap.keys():
        oo_mat.loc[p.split("|")[0], p.split("|")[1]] = observed_overlap[p]
        eo_mat.loc[p.split("|")[0], p.split("|")[1]] = expected_overlap[p]
    
    if OE_transform:
        epsilon = 1e-3

        return np.log(
            (oo_mat + epsilon) / (eo_mat + epsilon)
        )
    
    else:
        for k in oo_mat.index:
            oo_mat.loc[k, :] = oo_mat.loc[k, :] / (coverage1[k]*len(loci_1))
        return oo_mat


def Hungarian_algorithm(matrix, conf_or_dis='conf'):

    if conf_or_dis == 'conf':
        confusion_matrix = np.array(matrix)
        best_assignments = linear_sum_assignment(confusion_matrix, maximize=True)

        print('Sum of optimal assignment sets / Sum of confusion matrix = {}/{}'.format(
            str(confusion_matrix[best_assignments[0],best_assignments[1]].sum()), 
            str(confusion_matrix.sum())))

        assignment_pairs = [(i, best_assignments[1][i]) for i in range(len(best_assignments[0]))]
        return assignment_pairs

    elif conf_or_dis == 'dist':
        distance_matrix = np.array(matrix)
        best_assignments = linear_sum_assignment(distance_matrix, maximize=False)

        assignment_pairs = [(i, best_assignments[1][i]) for i in range(len(best_assignments[0]))]
        return assignment_pairs


def connect_bipartite(loci_1, loci_2, assignment_matching, mnemon=True):
    corrected_loci_1 = loci_1.iloc[:, :3]
    corrected_loci_2 = loci_2.iloc[:, :3]

    for i in range(len(assignment_matching)):
        new_name = str(assignment_matching[i][0]) + "|" + str(assignment_matching[i][1])

        if mnemon:
            corrected_loci_1[new_name] = loci_1['posterior'+str(assignment_matching[i][0].split('_')[0])]
            corrected_loci_2[new_name] = loci_2['posterior'+str(assignment_matching[i][1].split('_')[0])]
        else:
            corrected_loci_1[new_name] = loci_1['posterior'+str(assignment_matching[i][0])]
            corrected_loci_2[new_name] = loci_2['posterior'+str(assignment_matching[i][1])]
    
    return corrected_loci_1, corrected_loci_2
