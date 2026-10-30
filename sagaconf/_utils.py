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

    if array2D.ndim == 2 and array2D.dtype.kind == "f" and array2D.size > 0 and not np.isnan(array2D).any():
        # values that are all float16 numbers (the loaded posteriors) index a table by their bits
        with np.errstate(all="ignore"):
            half = array2D.astype(np.float16)
        if np.array_equal(half, array2D):
            bits = half.view(np.uint16)
            present = np.zeros(1 << 16, dtype=bool)
            present[bits] = True
            codes = np.flatnonzero(present)
            values = codes.astype(np.uint16).view(np.float16).astype(array2D.dtype)
            table = np.empty(1 << 16, dtype=array2D.dtype)
            for k in range(len(codes)):
                table[codes[k]] = logit(values[k])
            return table[bits]

        # logit depends only on the value: apply it once per distinct value
        uniq, inverse = np.unique(array2D, return_inverse=True)
        mapped = np.empty(len(uniq), dtype=array2D.dtype)
        for k in range(len(uniq)):
            mapped[k] = logit(uniq[k])
        return mapped[inverse].reshape(array2D.shape)

    for i in range(array2D.shape[0]):
        for j in range(array2D.shape[1]):
            array2D[i,j] = logit(array2D[i,j])

    return array2D


def to_tsv(df, path, header, na_rep="", chunk=1 << 20):
    """Same file as df.to_csv(path, sep='\t', header=header, index=False, na_rep=na_rep).

    Columns are formatted as to_csv formats them (floats and ints with astype(str), NaN as
    na_rep) and joined row by row, a chunk at a time; a frame with any other kind of column, or
    with a string to_csv would quote, is written by to_csv.
    """
    names = [str(c) for c in df.columns]
    arrays = [df.iloc[:, j].to_numpy() for j in range(df.shape[1])]
    plain = len(df) > 0 and len(names) > 1 and all(type(c) is str for c in df.columns)
    for v in arrays:
        if not plain:
            break
        if v.dtype == object:
            plain = all(type(x) is str for x in pd.unique(v))
            plain = plain and not any(_quoted(x) for x in pd.unique(v))
        else:
            plain = v.dtype.kind in "iuf"
    if not plain or any(_quoted(c) for c in names):
        return df.to_csv(path, sep='\t', header=header, index=False, na_rep=na_rep)

    with open(path, "w", encoding="utf-8", newline="") as fh:
        if header:
            fh.write("\t".join(names) + "\n")
        for start in range(0, len(df), chunk):
            cols = []
            for v in arrays:
                v = v[start:start + chunk]
                if v.dtype.kind == "f":
                    s = v.astype(str)
                    nan = np.isnan(v)
                    if nan.any():
                        s = s.astype(object)
                        s[nan] = na_rep
                    v = s
                elif v.dtype.kind in "iu":
                    v = v.astype(str)
                cols.append(v.tolist())
            fh.write("\n".join(map("\t".join, zip(*cols))) + "\n")


def _quoted(s):
    return s == "" or any(ch in s for ch in '\t"\n\r')
