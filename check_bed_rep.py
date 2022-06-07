"""
README:

this script will compare two bedfiles (assumably from replicated or roughly similar runs)
and outputs a boolean label 1 for every genomic position if that annotation is replicated between the two or 0 otherwise.
Also, saves the results in bed format.
"""

import pandas as pd
import numpy as np
import multiprocessing as mp
from scipy.optimize import linear_sum_assignment
import sys

def inplace_binning(bedfile, resolution): 
    # posterior_ID = bedfile.split('/')[-1].replace('.bed', '')
    with open(bedfile,"r") as bed_file:
        lines = bed_file.readlines()
        annot_list = []
        for il in range(1, int(len(lines))):
            l = lines[il]
            ll = l.split('\t')
            ll = ll[:4]
            annot_list.append(ll)

    new_annot_list = [] 
    for i in range(len(annot_list)):
        if int(annot_list[i][2]) - int(annot_list[i][1]) > resolution:
            inner_list = []

            for j in range(int(annot_list[i][1]), int(annot_list[i][2]), resolution):
                if int(j + resolution) > int(annot_list[i][2]):
                    inner_list.append([annot_list[i][0], j, (annot_list[i][2]), int(annot_list[i][3])])
                else:
                    inner_list.append([annot_list[i][0], j, int(j + resolution), int(annot_list[i][3])])

            for k in inner_list:
                new_annot_list.append(k)
        else:
            new_annot_list.append(annot_list[i])
        
    new_annot_list = pd.DataFrame(new_annot_list, columns=['chr', 'start', 'end', "label"])
    return new_annot_list

def intersect(df1, df2):
    if "_1" in df1.iloc[0, 0] or "_2" in df1.iloc[0, 0]:
        chrdf1 = list(df1.chr)
        for i in range(len(chrdf1)):
            if "_1" in chrdf1[i]:
                chrdf1[i] = chrdf1[i].replace("_1", "")
            elif "_2" in chrdf1[i]:
                chrdf1[i] = chrdf1[i].replace("_2", "")
        df1.chr = np.array(chrdf1)

    if "_1" in df2.iloc[0, 0] or "_2" in df2.iloc[0, 0]:
        chrdf2 = list(df2.chr)
        for i in range(len(chrdf2)):
            if "_2" in chrdf2[i]:
                chrdf2[i] = chrdf2[i].replace("_2", "")
            elif "_1" in chrdf2[i]:
                chrdf2[i] = chrdf2[i].replace("_1", "")
        df2.chr = np.array(chrdf2)

    intersect = pd.merge(
        df1, 
        df2, 
        how='inner', on=['chr', 'start', 'end'])

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

    return df1, df2


def match_labels_and_get_reproducibility(rep1, rep2):
    
    print('generating confusion matrix')    

    num_labels = rep1.label.astype(int).max() - rep1.label.astype(int).min() + 1
    print("number of labels =   ", num_labels)

    labels_range = [rep1.label.astype(int).min(), rep1.label.astype(int).max()+1]
    print("matching labels")
    observed_overlap = np.zeros((num_labels, num_labels))
    expected_overlap = np.zeros((num_labels, num_labels))

    r1_label_bin_count = {}
    r2_label_bin_count = {}

    for i in range(rep1.label.astype(int).min(), rep1.label.astype(int).max()+1):
        r1_label_bin_count[i] = 0
        r2_label_bin_count[i] = 0
    
    rep1_labels = np.array(rep1.label)
    rep2_labels = np.array(rep2.label)
        
    for i in range(len(rep1)):
        try:
            r1_label_bin_count[int(rep1_labels[i])] += 1
            r2_label_bin_count[int(rep2_labels[i])] += 1

            if int(labels_range[0]) == 0:
                observed_overlap[int(rep1_labels[i]), int(rep2_labels[i])] += 1

            elif int(labels_range[0]) == 1:
                observed_overlap[int(rep1_labels[i])-1, int(rep2_labels[i])-1] += 1

        except:
            pass
    
    observed_overlap = pd.DataFrame(observed_overlap,
        columns=list(range(rep1.label.astype(int).min(), rep1.label.astype(int).max()+1)),
        index=list(range(rep1.label.astype(int).min(), rep1.label.astype(int).max()+1)))

    epsilon = 1e-3
    
    # calculate ratio of coverage for each label
    for k, v in r1_label_bin_count.items():
        r1_label_bin_count[k] = float(v) /  rep1.shape[0]

    for k, v in r2_label_bin_count.items():
        r2_label_bin_count[k] = float(v) /  rep2.shape[0]

    for i in range(rep1.label.astype(int).min(), rep1.label.astype(int).max()+1):
        for j in range(rep1.label.astype(int).min(), rep1.label.astype(int).max()+1):
            
            if labels_range[0] == 0:
                expected_overlap[i, j] = (r1_label_bin_count[i] * r2_label_bin_count[j]) * rep1.shape[0]

            elif labels_range[0] == 1:
                expected_overlap[i-1, j-1] = (r1_label_bin_count[i] * r2_label_bin_count[j]) * rep1.shape[0]
            
    expected_overlap = pd.DataFrame(
        expected_overlap, 
        columns=list(range(rep1.label.astype(int).min(), rep1.label.astype(int).max()+1)), 
        index=list(range(rep1.label.astype(int).min(), rep1.label.astype(int).max()+1)))

    oe_overlap = (observed_overlap + epsilon) / (expected_overlap + epsilon)
    oe_overlap = np.log(oe_overlap)

    
    confusion_matrix = np.array(oe_overlap)
    best_assignments = linear_sum_assignment(confusion_matrix, maximize=True)

    print('Sum of optimal assignment sets / Sum of confusion matrix = {}/{}'.format(
        str(confusion_matrix[best_assignments[0],best_assignments[1]].sum()), 
        str(confusion_matrix.sum())))

    assignment_pairs = {}
    for i in range(len(best_assignments[0])):
        if labels_range[0] == 0:
            assignment_pairs[int(i)] = int(best_assignments[1][i])

        elif labels_range[0] == 1:
            assignment_pairs[int(i)+1] = int(best_assignments[1][i])+1
    
    inv_assignment_pairs = {v: k for k, v in assignment_pairs.items()}
            
    print("label matching pairs:  ", inv_assignment_pairs)
    reprod_list_1 = []
    reprod_list_2 = []
    print("getting the reproducibility results")
    for i in range(len(rep1)):
        label_1_i = rep1_labels[i]
        label_2_i = rep2_labels[i]

        if int(label_1_i) == int(inv_assignment_pairs[int(label_2_i)]):
            reprod_list_1.append("1")
        else:
            reprod_list_1.append("0")
        
        if int(label_2_i) == int(assignment_pairs[int(label_1_i)]):
            reprod_list_2.append("1")
        else:
            reprod_list_2.append("0")
        
    rep1["reproduced"] = np.array(reprod_list_1)
    rep2["reproduced"] = np.array(reprod_list_2)

    return rep1, rep2

if __name__=="__main__":
    rep1_dir = sys.argv[1]
    rep2_dir = sys.argv[2]
    resolution = int(sys.argv[3])
    rep1_res = sys.argv[4]
    rep2_res = sys.argv[5]

    print("binning into resolution-sized bins")
    rep1 = inplace_binning(rep1_dir, resolution)
    rep2 = inplace_binning(rep2_dir, resolution)

    print("getting the intersection")
    rep1, rep2 = intersect(rep1, rep2)

    rep1, rep2 = match_labels_and_get_reproducibility(rep1, rep2)
    print("ratio of reproducibility for replicate 1=    ",len(rep1.loc[(rep1["reproduced"]=="1"), :])/len(rep1))
    print("ratio of reproducibility for replicate 2=    ",len(rep2.loc[(rep2["reproduced"]=="1"), :])/len(rep2))

    print("saving the results")
    rep1.to_csv(rep1_res, sep="\t")
    rep1.to_csv(rep2_res, sep="\t")