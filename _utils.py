from math import log
import pandas as pd
import numpy as np
from scipy.optimize import linear_sum_assignment
import os 
import multiprocessing as mp
import matplotlib.pyplot as plt

def read_bedgraph(bg_file):
    df = []
    with open(bg_file, 'r') as bgf:
        lines = bgf.readlines()
    for i in range(len(lines)): 
        if i!=0:
            tl = lines[i].split('\t')
            tl[-1] = float(tl[-1].replace('\n', ''))/100
            df.append(tl)
    return pd.DataFrame(df, columns=['chr', 'start', 'end', 'posterior'])

def define_positions(regions_file, windowsize):
    windowsize = int(windowsize)
    with open(regions_file, 'r') as rf:
        lines = rf.readlines()
        posmap = ([i.split('\t')[:-1] for i in lines])
    
    loci = []
    for idx in range(len(posmap)):
        for st in range(int(posmap[idx][1]), int(posmap[idx][2]), windowsize):

            if st + windowsize - 1 > int(posmap[idx][2]):
                loci.append([posmap[idx][0], st, int(posmap[idx][2])])
                
            else:
                loci.append([posmap[idx][0], st, st + windowsize - 1])
    
    loci = pd.DataFrame(loci, columns=['chr', 'window_start', 'window_end'])
    return loci

def align_scores(input_dict):
    # parse_input
    k, v, loci, windowsize, num_labels = input_dict['k'], input_dict['v'], input_dict['loci'], input_dict['windowsize'], input_dict['num_labels']
    
    print('aligning posteriorfile -> {}/{}. number of segments: {}'.format(str(int(k)+1), str(num_labels), str(len(v))))
    loci.insert(
        loc=loci.shape[1], column='posterior'+str(k), 
        value=[None for _ in range(loci.shape[0])]
        )

    # print(loci)

    num_filled = 0
    for ind in range(loci.shape[0]):
        
        if num_filled % 10000 == 0:
            print('{}th index, filled {}/{}'.format(ind, num_filled, loci.shape[0]))

        subset_df = v.loc[v['chr'] == loci['chr'][ind]]
        subset_df = subset_df.reset_index(drop=True)

        bin_filling_status = False
        for ssdf in range(len(subset_df)):
            statement1 = bool(int(subset_df['start'][ssdf]) <= int(loci['window_start'][ind]) <= int(subset_df['end'][ssdf]))
            statement2 = bool(int(loci['window_start'][ind]) <= int(subset_df['start'][ssdf]) <= int(loci['window_end'][ind]))

            if statement1 or statement2:
                range_ref_map = range(int(loci['window_start'][ind]), int(loci['window_end'][ind]))
                range_posteri = range(int(subset_df['start'][ssdf]), int(subset_df['end'][ssdf]))
                
                set_r1 = set(range_ref_map)
                overlap = set_r1.intersection(range_posteri)

                if len(overlap) >= windowsize/2:
                    loci.loc[ind, 'posterior'+str(k)] = subset_df['posterior'][ssdf]
                    num_filled += 1
                    bin_filling_status = True

        if bin_filling_status == False:
            print("could not fill {}th bin".format(ind))
        
    return loci['posterior'+str(k)]


def position_labels(posterior_df_list, position_map, windowsize, num_labels, n_subset=None): 
    windowsize = int(windowsize)
    '''
    chooses the predicted label for each position based on posterior 
    probability value.
    '''
    # define all positions (windowsize=resolution)

    loci = position_map

    if n_subset != None:
        loci = loci.iloc[list(range(0, loci.shape[0], int(loci.shape[0]/n_subset))), :]
        loci = loci.reset_index(drop=True)

    print('number of positions: {}'.format(loci.shape[0]))

    # align scores to windows
    mp_inputs = []
    
    # k as in posterior_k, and v is the k-th posterior df
    for k, v in posterior_df_list.items():
        mp_inputs.append({'k':k, 'v':v, 'loci':loci, 'windowsize':windowsize, 'num_labels':num_labels})

    with mp.Pool(num_labels) as pool_obj:
        mpresults = pool_obj.map(align_scores, mp_inputs)
    
    mpresults = pd.concat(mpresults, axis=1)
    loci = pd.concat([loci, mpresults], axis=1)
    return loci
    
def find_max_posteri(loci, num_labels=10):
    loci_posteri =  loci.loc[:, ['posterior'+str(i) for i in range(int(num_labels))]]
    for lpc in loci_posteri.columns:
        loci_posteri[lpc] = loci_posteri[lpc].astype('float')

    max_posteri = loci_posteri.idxmax(axis=1)
    return max_posteri

def confusion_matrix(loci_1, loci_2, num_labels):
    num_labels = int(num_labels)
    max_1posteri =  find_max_posteri(loci_1, num_labels=num_labels)
    max_2posteri =  find_max_posteri(loci_2, num_labels=num_labels)

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