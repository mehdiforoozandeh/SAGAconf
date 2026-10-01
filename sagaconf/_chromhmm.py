import os
import pandas as pd


def chrhmm_initialize_bin(chrom, numbins, res):
    empty_bins = []
    next_start = 0
    for _ in range(numbins):
        empty_bins.append(
            [chrom, next_start, int(next_start+res)])
        next_start = int(next_start+res)
    empty_bins = pd.DataFrame(empty_bins, columns=['chr', 'start', 'end'])
    empty_bins['start'] = empty_bins['start'].astype("int32")
    empty_bins['end'] = empty_bins['end'].astype("int32")
    return empty_bins


def read_posterior_file(filepath):
    with open(filepath,'r') as posteriorfile:
        lines = posteriorfile.readlines()
    vals = []
    for il in range(len(lines)):
        ilth_vals = lines[il].split('\t')
        ilth_vals[-1] = ilth_vals[-1].replace("\n","")
        vals.append(ilth_vals)
    vals = pd.DataFrame(vals[2:], columns=["posterior{}".format(i.replace("E","")) for i in vals[1]])
    vals = vals.astype("float32")
    return vals


def ChrHMM_read_posteriordir(posteriordir, resolution=200):
    '''
    for each file in posteriordir
    Initialize emptybins based on chromsizes
    fill in the posterior values for each slot
    return DF
    '''
    ls = os.listdir(posteriordir)
    to_parse = []
    for f in ls:
        # if rep in f:
        to_parse.append(f)  

    parsed_posteriors = {}
    for f in to_parse:
        fileinfo = f.split("_")
        if "chr" in fileinfo[-2]:
            posteriors = read_posterior_file(posteriordir + '/' + f)
            bins = chrhmm_initialize_bin(fileinfo[-2], len(posteriors), resolution)
            posteriors = pd.concat([bins, posteriors], axis=1)
            parsed_posteriors[fileinfo[-2]] = posteriors 
        
    parsed_posteriors = pd.concat([parsed_posteriors[c] for c in sorted(list(parsed_posteriors.keys()))], axis=0)
    parsed_posteriors = parsed_posteriors.reset_index(drop=True)
    return parsed_posteriors
