import functools
import multiprocessing as mp
import numpy as np
import os
import pandas as pd


def inplace_binning(posterior_file, resolution): 
    posterior_ID = posterior_file.split('/')[-1].replace('.bedGraph', '')
    with open(posterior_file,"r") as posterior_file:
        lines = posterior_file.readlines()
        posterior_list = []
        for il in range(1, len(lines)):
            l = lines[il]
            ll = l.split('\t')
            ll[-1] = float(int(ll[-1].rstrip("\n"))/100)  # the last line may lack a newline
            posterior_list.append(ll)

    new_posterior_list = []
    for i in range(len(posterior_list)):
        if int(posterior_list[i][2]) - int(posterior_list[i][1]) > resolution:
            inner_list = []

            for j in range(int(posterior_list[i][1]), int(posterior_list[i][2]), resolution):
                if int(j + resolution) > int(posterior_list[i][2]):
                    inner_list.append([posterior_list[i][0], j, (posterior_list[i][2]), posterior_list[i][3]])
                else:
                    inner_list.append([posterior_list[i][0], j, int(j + resolution), posterior_list[i][3]])

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
    # The state columns are joined by row position below, so every file must have the same
    # bins; otherwise values land on the wrong coordinates. This check always runs, so
    # assert_coord_match is ignored (the old optional check could never fail).
    def coords(d):
        return (np.asarray(d['chr'].astype(str)), np.asarray(pd.to_numeric(d['start'])),
                np.asarray(pd.to_numeric(d['end'])))

    def describe(c, f, i):
        if i >= len(c[0]):
            return "no bin in %s" % f
        return "%s:%d-%d in %s" % (c[0][i], c[1][i], c[2][i], f)

    ref_file, ref = os.path.basename(posterior_file_list[0]), coords(parsed_list[0])
    for pfile, d in zip(posterior_file_list[1:], parsed_list[1:]):
        c = coords(d)
        n = min(len(ref[0]), len(c[0]))
        same = (ref[0][:n] == c[0][:n]) & (ref[1][:n] == c[1][:n]) & (ref[2][:n] == c[2][:n])
        if same.all() and len(ref[0]) == len(c[0]):
            continue
        i = int(np.argmin(same)) if not same.all() else n
        f = os.path.basename(pfile)
        raise ValueError(
            "Segway posterior files are not aligned: %s has %d bins and %s has %d bins after "
            "binning at %d bp. First difference at bin %d: %s, but %s. Every posterior file must "
            "cover the same regions with the same bin boundaries."
            % (ref_file, len(ref[0]), f, len(c[0]), resolution, i,
               describe(ref, ref_file, i), describe(c, f, i)))

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
