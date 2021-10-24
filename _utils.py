from math import log
import pandas as pd
import numpy as np
import shutil
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

def read_posteriordir(posteriordir):
    ls_dir = os.listdir(posteriordir)
    ls = []

    for f in ls_dir:
        if ".bedGraph" in f:
            ls.append(f)

    posterior_df_list = {}

    for f in ls:
        posterior_id = f.replace('.bedGraph','')
        posterior_df_list[posterior_id.replace('posterior','')] = read_bedgraph(posteriordir+'/'+f)

    print(posterior_df_list.keys())

    return posterior_df_list #dictionary actually :))

"""
name_sig is the name_signiture of each run
which basically indicates the name of the 
final result directory containing posterior 
values, segway.bed, and other analysis.
"""

def gather_results(traindir, posteriordir, name_sig):
    os.system('mkdir {}'.format(name_sig))
    os.system('cp {}/*.gz {}'.format(posteriordir, name_sig))
    os.system('gzip -d {}/*.gz'.format(name_sig))
    os.system('cp {}/log/segway.sh {}/segway.sh'.format(traindir, name_sig))

def clean_up(traindir, posteriordir):
    shutil.rmtree(traindir)
    shutil.rmtree(posteriordir)

def run_segtools(name_sig):
    os.system("cd {} && segtools-length-distribution segway.bed".format(name_sig))
    os.system("cd {} && segtools-gmtk-parameters params.params".format(name_sig))

'''
TERMINOLOGY NOTE: posterior_df_list IS THE SAME AS bg_dfs
'''

def QC(posterior_df_list):
    '''
    returns 1 - (ratio of (bad values * length))
    higher values indicate better quality
    '''

    total_length_bp = 0
    bad_val_lengthcount = 0

    test_df = posterior_df_list['0']
    for j in range(len(test_df)):
        total_length_bp += int(test_df['end'][j]) - int(test_df['start'][j])


    for post_df in posterior_df_list.values():
        for i in range(len(post_df)):

            if post_df['posterior'][i] > 0.5:
                if post_df['posterior'][i] == 1:
                    bad_val_lengthcount += int(post_df['end'][i]) - int(post_df['start'][i])

    return 1 - float(bad_val_lengthcount/ total_length_bp)            
    
def diagnose_zero_one_ratio(bg_dfs):
    ones_count=0
    zeros_count=0
    posterior_count=0

    for i in bg_dfs.keys():
        post = np.array(bg_dfs[i].posterior)
        posterior_count += len(post)
        ones_count += len(post[post == 1])
        zeros_count += len(post[post == 0])
    
    return float(ones_count/posterior_count), float(zeros_count/posterior_count)

def diagnose_label_ratios(bg_dfs):
    label_ratios = {}
    for i in bg_dfs.keys():
        post = np.array(bg_dfs[i].posterior)
        label_ratios[i] = len(post[post >= 0.5])/len(post)
    
    return label_ratios

def read_include_file(file): 
    '''to parse the pilotregion file '''
    with open(file, 'r') as inc_file:
        lines = inc_file.readlines()
        
    include = []
    for l in lines:
        splitted = l.split('\t')[:-1]
        # splitted[-1] = splitted[-1].replace('\n', '')
        include.append(splitted)
    
    return pd.DataFrame(include, columns=['chr', 'start', 'end'])


def initialize_bins(coords, res): 
    '''initializes empty bins according to the genome positions specified in the pilotregions file'''

    supercontig_in_progress = []
    for i in range(len(coords)):
        region_in_progress = list(coords.iloc[i,:])
        for j in range(int(region_in_progress[1]), int(region_in_progress[2]), res):

            if int(j + res) > int(region_in_progress[2]):
                supercontig_in_progress.append(
                    [region_in_progress[0], int(j), int(region_in_progress[2]) + int(res - (int(region_in_progress[2]) % res))])

            else: 
                supercontig_in_progress.append([region_in_progress[0], int(j), int(j + res)])

    return pd.DataFrame(supercontig_in_progress, columns=['chr', 'start', 'end'])

def instance_binning(input_dict):
    """This function assigns the posterior value to each bin based on the corresponding 
    value in the posterior*.bedGraph file. I used a "search margin" to prevent the function from 
    performing in O(n^2). M denotes search margin and limits the inner loop to M indices 
    only, starting from where the last bin was aligned with the raw data. The final posterior 
    value of each bin is weighted by the length of alignment. """

    bins_df, posterior_df, M, posteriorID = input_dict['bins_df'], input_dict['posterior_df'], input_dict['M'], input_dict['posteriorID']

    bins_df.insert(3, posteriorID, np.zeros(len(bins_df)))
    notfilled = 0
    c = 0 # denotes the center of search space

    for i in range(len(bins_df)): # i denotes index of empty bin
        filledbool = False
        # define search space
        if i%100 == 0:
            print("filled {} bins. could not fill {} bins".format(i-notfilled, notfilled))
        
        if c < M:
            search_start_loc = 0
        else:
            search_start_loc = c - int(M/10)
    
        if c + M > len(posterior_df):
            search_end_loc = len(posterior_df)
        else:
            search_end_loc = c + M

        for j in range(search_start_loc, search_end_loc): # j denotes index of posterior segment

            if bins_df.iloc[i, 0] == posterior_df.iloc[j, 0]: # check chr match

                statement1 = bool(int(bins_df.iloc[i,1]) <= int(posterior_df.iloc[j, 1]) <= int(bins_df.iloc[i,2]) )
                statement2 = bool(int(posterior_df.iloc[j, 1]) <= int(bins_df.iloc[i,1]) <= int(posterior_df.iloc[j, 2]) )

                statement3 = bool(int(bins_df.iloc[i,2]) < int(posterior_df.iloc[j, 1])) # passed statement

                if statement1 or statement2:
                    bin_range = range(int(bins_df.iloc[i,1]), int(bins_df.iloc[i,2]))
                    signal_range = range(int(posterior_df.iloc[j, 1]), int(posterior_df.iloc[j, 2]))
                    
                    set_r1 = set(bin_range)
                    overlap = set_r1.intersection(signal_range)

                    if len(overlap) > 0:
                        filledbool = True
                        c = j # update the center of search space to the latest matched index
                        bin_len = int(bins_df.iloc[i, 2]) - int(bins_df.iloc[i, 1]) #almost always is equal to resolution
                        bins_df.iloc[i, 3] += float(posterior_df.iloc[j, 3] * (len(overlap) / bin_len))
                
                elif statement3:
                    break

        if filledbool == False:
            notfilled +=1

    return bins_df[posteriorID]

def mp_binning(posterior_df_list, empty_bins, M):
    print('number of positions: {}'.format(empty_bins.shape[0]))
    mp_inputs = []
    
    # k as in posterior_k, and v is the k-th posterior df
    for k, v in posterior_df_list.items():
        mp_inputs.append(
            {"bins_df": empty_bins, "posterior_df": v, "M": M, "posteriorID": str(empty_bins)})

    with mp.Pool(len(list(posterior_df_list.keys()))) as pool_obj:
        mpresults = pool_obj.map(instance_binning, mp_inputs)
    
    mpresults = pd.concat(mpresults, axis=1)
    parsed_bins = pd.concat([empty_bins, mpresults], axis=1)
    return parsed_bins
