from ast import parse
import os, sys
import pandas as pd
import numpy as np

def create_cellmarkfiletable(celltype_dir, out_filename, suffix_to_look_for="bedGraph"):
    '''tab delimited file each row contains the 
    cell type the all the associated marks'''
    walk_obj = os.walk(celltype_dir)
    dirs_and_subdirs = [x for x in walk_obj]
    with open(out_filename, 'w') as cmft_file:
        for i in range(1, len(dirs_and_subdirs)):
            for j in dirs_and_subdirs[i][2]:
                if suffix_to_look_for in j:
                    cmft_file.write('{}\t{}\t{}\n'.format(
                        celltype_dir.split('/')[-1], dirs_and_subdirs[i][0].split('/')[-1],
                        dirs_and_subdirs[i][0].split('/')[-1]+'/'+j
                    ))

def binarize_data(inputbeddir, cellmarkfiletable, outputdir, resolution=100, chromlength='CHROMSIZES/hg19.txt'):
    cmdline = "java -Xmx10g -jar ChromHMM.jar BinarizeBam -b {} {} {} {} {}".format(
        resolution, chromlength, inputbeddir, cellmarkfiletable, outputdir
    )
    os.system(cmdline)

def learnModel(binary_input_dir, output_dir, num_labels='16', assembly='hg19', n_threads='0', random_seed=73):
    learnmodel_cmdline = "java -Xmx10g -jar ChromHMM.jar LearnModel -init random -s {} -printposterior -p {} {} {} {} {}".format(
        random_seed, n_threads, binary_input_dir, output_dir, num_labels, assembly
    )
    os.system(learnmodel_cmdline)

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

def ChrHMM_read_posteriordir(posteriordir, rep, resolution=100):
    '''
    for each file in posteriordir
    Initialize emptybins based on chromsizes
    fill in the posterior values for each slot
    return DF
    '''
    ls = os.listdir(posteriordir)
    to_parse = []
    for f in ls:
        if rep in f:
            to_parse.append(f)
    parsed_posteriors = {}
    for f in to_parse:
        fileinfo = f.split("_")
        posteriors = read_posterior_file(posteriordir + '/' + f)
        bins = chrhmm_initialize_bin(fileinfo[3], len(posteriors), resolution)
        posteriors = pd.concat([bins, posteriors], axis=1)
        parsed_posteriors[fileinfo[3]] = posteriors 
        # print(fileinfo[3], posteriors.shape, end=" | ")
    parsed_posteriors = pd.concat([parsed_posteriors[c] for c in sorted(list(parsed_posteriors.keys()))], axis=0)
    # parsed_posteriors = parsed_posteriors.sort_values(by=["chr"])
    parsed_posteriors = parsed_posteriors.reset_index(drop=True)
    return parsed_posteriors

"""
1. relocate (mv/cp) .bam files of each cell type to a new directory called "chmmfiles"
2. create 3 cmft files at each chmmfiles dir :
    a. concat
    b. rep1
    c. rep2
3. run three instances of learnModel for each celltype
    a. concat
    b. rep1
    c. rep2
"""
def relocate_bamfiles(CellType_dir):
    pass

def ChromHMM_rep():
    pass

def ChromHMM_param_init():
    pass

def ChromHMM_concat():
    pass
