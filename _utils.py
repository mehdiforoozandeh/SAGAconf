import pandas as pd
import numpy as np
import shutil, os, functools, itertools
import multiprocessing as mp

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

    # print(posterior_df_list.keys())

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
    os.system('cp {}/params/params.params {}/params.params'.format(traindir, name_sig))

def clean_up(traindir, posteriordir):
    try:
        shutil.rmtree(traindir)
    except:
        pass

    try:
        shutil.rmtree(posteriordir)
    except:
        pass

def run_segtools(name_sig):
    os.system("cd {} && segtools-length-distribution segway.bed".format(name_sig))
    os.system("cd {} && segtools-gmtk-parameters params.params".format(name_sig))

'''
TERMINOLOGY_NOTE: posterior_df_list IS THE SAME AS bg_dfs
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

def read_chrom_sizes_file(file):
    sizes = {}
    with open(file, 'r') as sizesfile:
        lines = sizesfile.readlines()
        for l in lines:
            splitted_l = l.split('\t')
            sizes[splitted_l[0]] = int(splitted_l[1][:-1])
    return sizes

def initialize_bins_includefile(coords, res): 
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

def initialize_bins_sizesfile(coords, res): 
    '''initializes empty bins according to the chromosome sizes specified in hg38.chrom.sizes file'''

    supercontig_in_progress = []
    for chr, size in coords.items():
        for i in range(0, int(size), int(res)):
            if i + int(res) > int(size):
                supercontig_in_progress.append([chr, int(i), int(size)])
            else:
                supercontig_in_progress.append([chr, int(i), int(i) + int(res)])

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
        if i%5000 == 0:
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
            {"bins_df": empty_bins, "posterior_df": v, "M": M, "posteriorID": 'posterior'+str(k)})
        
    with mp.Pool(len(list(posterior_df_list.keys()))) as pool_obj:
        mpresults = pool_obj.map(instance_binning, mp_inputs)
    
    mpresults = pd.concat(mpresults, axis=1)
    parsed_bins = pd.concat([empty_bins, mpresults], axis=1)
    return parsed_bins


def check_coord_match(df_1, df_2):
    matched_sofar=True

    if len(df_1) == len(df_2):
        for i in range(len(df_1)):
            match_statement = bool(
                df_1['chr'][i] == df_2['chr'][i] and
                df_1['start'][i] == df_2['start'][i] and
                df_1['end'][i] == df_2['end'][i] 
            )
            if match_statement == False:
                matched_sofar == False
                break
    else:
        matched_sofar == False
    
    return matched_sofar

def inplace_binning(posterior_file, resolution): 
    posterior_ID = posterior_file.split('/')[-1].replace('.bedGraph', '')
    with open(posterior_file,"r") as posterior_file:
        lines = posterior_file.readlines()
        posterior_list = []
        for il in range(1, len(lines)):
            l = lines[il]
            ll = l.split('\t')
            ll[-1] = float(int(ll[-1][:-1])/100)
            posterior_list.append(ll)

    new_posterior_list = []
    for i in range(len(posterior_list)):
        if int(posterior_list[i][2]) - int(posterior_list[i][1]) > resolution:
            inner_list = []

            for j in range(int(posterior_list[i][1]), int(posterior_list[i][2]), resolution):
                if int(j + resolution) > int(posterior_list[i][2]):
                    inner_list.append([posterior_list[i][0], j, (posterior_list[i][2]), int(posterior_list[i][3])])
                else:
                    inner_list.append([posterior_list[i][0], j, int(j + resolution), int(posterior_list[i][3])])

            for k in inner_list:
                new_posterior_list.append(k)
        else:
            new_posterior_list.append(posterior_list[i])

        if i % int(len(posterior_list)/5000) == 0:
            print(len(new_posterior_list))

    new_posterior_list = pd.DataFrame(new_posterior_list, columns=['chr', 'start', 'end', posterior_ID])
    return new_posterior_list

def mp_inplace_binning(posterior_dir, resolution, assert_coord_match=False, mp=False):
    ls_dir = os.listdir(posterior_dir)
    posterior_file_list = []

    for f in ls_dir:
        if ".bedGraph" in f:
            posterior_file_list.append(f)

    posterior_file_list = ["{}/posterior{}.bedGraph".format(posterior_dir, i) for i in range(len(posterior_file_list))]

    print("starting the parsing")
    if mp:
        p_obj = mp.Pool(len(posterior_file_list))
        parsed_list = p_obj.map(
            functools.partial(inplace_binning, resolution=resolution), posterior_file_list)
    else:
        parsed_list = []
        for pfile in posterior_file_list:
            parsed_list.append(inplace_binning(pfile, resolution))

    print("parsed all individual files")
    match_statement = True
    if assert_coord_match:
        print("asserting coordinate match between different posteriorIDs")
        for comb in itertools.combinations(range(len(parsed_list)), 2):
            if check_coord_match(parsed_list[comb[0]], parsed_list[comb[1]]) == False:
                match_statement = False

    if match_statement:
        print("concatenating results")
        parsed_posterior = [parsed_list[0]['chr'], parsed_list[0]['start'], parsed_list[0]['end']]
        for d in parsed_list:
            parsed_posterior.append(d.iloc[:, -1])
        
        parsed_posterior = pd.concat(parsed_posterior, axis=1)
    
        return parsed_posterior


def logit(p):

    if 1e-16 < p < (1- 1e-16):
        return np.log((p)/(1-p))

    elif p < 1e-16:
        return -36.84136

    elif (1- 1e-16) < p:
        return 36.84136

def logit_array(array2D):
    array2D = np.array(array2D)

    for i in range(array2D.shape[0]):
        for j in range(array2D.shape[1]):
            array2D[i,j] = logit(array2D[i,j])

    return array2D