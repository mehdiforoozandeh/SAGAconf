import functools
import itertools
import multiprocessing as mp
import numpy as np
import os
import pandas as pd


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

        # if i % int(len(posterior_list)/1000000) == 0:
        #     print(len(new_posterior_list))

    new_posterior_list = pd.DataFrame(new_posterior_list, columns=['chr', 'start', 'end', posterior_ID])
    return new_posterior_list


def mp_inplace_binning(posterior_dir, resolution, assert_coord_match=False, m_p=False):
    ls_dir = os.listdir(posterior_dir)
    posterior_file_list = []

    for f in ls_dir:
        if ".bedGraph" in f:
            posterior_file_list.append(f)

    posterior_file_list = ["{}/posterior{}.bedGraph".format(posterior_dir, i) for i in range(len(posterior_file_list))]

    print("starting the parsing")
    if m_p:
        p_obj = mp.Pool(len(posterior_file_list))
        parsed_list = p_obj.map(
            functools.partial(inplace_binning, resolution=resolution), posterior_file_list)
    else:
        parsed_list = []
        for pfile in posterior_file_list:
            print(pfile)
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
