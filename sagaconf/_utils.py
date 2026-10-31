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


def to_tsv(df, path, header, na_rep=""):
    """Same file as df.to_csv(path, sep='\t', header=header, index=False, na_rep=na_rep).

    Columns are formatted as to_csv formats them (floats and ints with astype(str), NaN as
    na_rep) and joined row by row, a chunk at a time, ending lines with os.linesep as to_csv
    does. A frame with any other kind of column, or with a string to_csv would quote, a path
    to_csv would compress, or an na_rep longer than "nan" (to_csv may cut it), go to to_csv.
    """
    names = [str(c) for c in df.columns]
    arrays = [df.iloc[:, j].to_numpy() for j in range(df.shape[1])]
    plain = (isinstance(path, str) and len(df) > 0 and len(names) > 1 and len(na_rep) <= 3
             and all(type(c) is str for c in df.columns)
             and not path.lower().endswith((".gz", ".bz2", ".zip", ".xz", ".zst", ".tar")))
    for v in arrays:
        if not plain:
            break
        if v.dtype == object:
            try:
                u = pd.unique(v)
            except TypeError:
                u = [None]
            plain = all(type(x) is str for x in u) and not any(_quoted(x) for x in u)
        else:
            plain = v.dtype.kind in "iuf"
    if not plain or any(_quoted(c) for c in names):
        return df.to_csv(path, sep='\t', header=header, index=False, na_rep=na_rep)

    eol = os.linesep
    with open(path, "w", encoding="utf-8", newline="") as fh:
        if header:
            fh.write("\t".join(names) + eol)
        for start in range(0, len(df), 1 << 20):
            cols = []
            for v in arrays:
                v = v[start:start + (1 << 20)]
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
            fh.write(eol.join(map("\t".join, zip(*cols))) + eol)


def _quoted(s):
    # to_csv quotes a field with the separator, a quote or a line break, and an empty field
    return s == "" or any(ch in s for ch in '\t"\n\r')
