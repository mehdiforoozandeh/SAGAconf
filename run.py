from _utils import *
from _pipeline import *
from _reproducibility import *
from _interpret import *
from _visualize import *
from get_data import *
from matchlabels import *
import pandas as pd
import glob, os 
import multiprocessing as mp
from functools import partial

'''

In this file, the step-by-step execution of analysis-funcitons will be managed.
Basically, by importing and managing functions from different modules
and classes, we will execute every step of the methods section. 

1 - get data from encode : USE 
    search_encode(cell, download_dir, target_assembly="GRCh38")

SEGWAY:
    2 - create genomedata (separate and concatenated) files and trackname_assay.txt
    3 - run segway (separate and concatenated)

ChromHMM:
    2 - Binarize data
    3 - run chromhmm
    
4 - parse and bin results
4.5  - biological label interpretation
5 - match labels (for separate runs)
6 - cluster labels
7 - reproducibility analysis (at different levels of abstraction)

'''

def download_encode_files(celltype_list, download_dir, target_assembly):
    for celltype in celltype_list:
        search_encode(celltype, download_dir, target_assembly=target_assembly)

def check_if_data_exists(celltype_list, download_dir):
    exists_bool = []
    for ct in celltype_list:
        if ' ' in ct:
            ct = ct.replace(' ', '_')

        if os.path.isdir(download_dir+'/'+ct):
            exists_bool.append(True)

        else:
            exists_bool.append(False)

    return exists_bool

def read_list_of_assays(celltype_dir):
    ls = os.listdir(celltype_dir)
    list_of_assays = []
    for d in ls:
        if os.path.isdir(celltype_dir + '/' + d):
            list_of_assays.append(d)
    
    return list_of_assays

def Convert_all_BW2BG(celltype_dir):
    '''
    Creats BedGraph files from all BigWig files
    found within a celltype directory
    '''

    ls = os.listdir(celltype_dir)
    if " " in celltype_dir:
        celltype_dir = celltype_dir.replace(" ","\ ")
    for f in ls:
        if ".bigWig" in f:
            if f.replace(".bigWig",".bedGraph") not in ls:
                os.system(
                    "./bigWigToBedGraph {}.bigWig {}.bedGraph".format(
                        celltype_dir+'/'+f.replace(".bigWig",""), 
                        celltype_dir+'/'+f.replace(".bigWig","")))

def read_metadata(donwload_dir):
    metadata = {}
    for ct in CellType_list:
        tracks_list = [tr for tr in os.listdir(download_dir+ct) if os.path.isdir(download_dir+ct+'/'+tr)]
        for tr in tracks_list:
            tfmd = pd.read_csv(download_dir+ct+'/'+tr+'/track_files_metadata.csv')
            tfmd.index = list(tfmd['Unnamed: 0'])
            tfmd = tfmd.drop('Unnamed: 0', axis=1)
            if ct in list(metadata.keys()):
                metadata[ct].append(tfmd)
            else:
                metadata[ct] = [tfmd]

    return metadata

def gather_replicates(celltype_dir):
    '''
    for each celltype
    create two files (one for each replicate)
    each file contains the following

    segway_rep1_{biosampleID}.txt:
        {ASSAY 1}     {File 1}
        {ASSAY 2}     {File 2}
        ...

    segway_rep2_{biosampleID}.txt:
        {ASSAY 1}     {File 1}
        {ASSAY 2}     {File 2}
        ...
    
    SAME TWO FILES FOR CHROMHMM...
    '''

    tracks_list = [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)]
    with open(celltype_dir+'/replicate_number.txt', 'r') as repguide_file:
        lines = repguide_file.readlines()
        rep1_biosampleID = lines[0][5:-1]
        rep2_biosampleID = lines[1][5:]

    print("rep1_biosampleID:\t{}\nrep2_biosampleID:\t{}\n".format(rep1_biosampleID, rep2_biosampleID))
    nav_rep1_seg = {}
    nav_rep2_seg = {}

    nav_rep1_chmm = {}
    nav_rep2_chmm = {}

    for tr in tracks_list:
        tfmd = pd.read_csv(celltype_dir+'/'+tr+'/track_files_metadata.csv')
        tfmd.index = list(tfmd['Unnamed: 0'])
        tfmd = tfmd.drop('Unnamed: 0', axis=1)

        assert str(tfmd.loc['biosample', 'rep1_fcoc']) == rep1_biosampleID
        assert str(tfmd.loc['biosample', 'rep2_fcoc']) == rep2_biosampleID

        nav_rep1_seg[tfmd.loc['assay', 'rep1_fcoc']] = '{}/{}/{}.bedGraph'.format(
            celltype_dir, tr, tfmd.loc['accession', 'rep1_fcoc'])

        nav_rep2_seg[tfmd.loc['assay', 'rep2_fcoc']] = '{}/{}/{}.bedGraph'.format(
            celltype_dir, tr, tfmd.loc['accession', 'rep2_fcoc'])

        nav_rep1_chmm[tfmd.loc['assay', 'rep1_spv']] = '{}/{}/{}.bedGraph'.format(
            celltype_dir, tr, tfmd.loc['accession', 'rep1_spv'])

        nav_rep2_chmm[tfmd.loc['assay', 'rep2_spv']] = '{}/{}/{}.bedGraph'.format(
            celltype_dir, tr, tfmd.loc['accession', 'rep2_spv'])

    with open(celltype_dir+"/seg_rep1_{}.txt".format(rep1_biosampleID), "w") as segrep1_file:
        for k, v in nav_rep1_seg.items():
            segrep1_file.write("{}\t{}\n".format(k, v))

    with open(celltype_dir+"/seg_rep2_{}.txt".format(rep2_biosampleID), "w") as segrep2_file:
        for k, v in nav_rep2_seg.items():
            segrep2_file.write("{}\t{}\n".format(k, v))

    with open(celltype_dir+"/chmm_rep1_{}.txt".format(rep1_biosampleID), "w") as chmmrep1_file:
        for k, v in nav_rep1_chmm.items():
            chmmrep1_file.write("{}\t{}\n".format(k, v))

    with open(celltype_dir+"/chmm_rep2_{}.txt".format(rep2_biosampleID), "w") as chmmrep2_file:
        for k, v in nav_rep2_chmm.items():
            chmmrep2_file.write("{}\t{}\n".format(k, v))

def create_genomedata(celltype_dir, sequence_file):
    '''
    read segway replicate files.
    gather bedGraph file paths

    create genomedata file_rep1 in celltype_dir
    create genomedata file_rep2 in celltype_dir
    '''
    print('creating genomedata files for {}...'.format(celltype_dir))
    tracks_bgs_rep1 = {}
    with open(glob.glob(celltype_dir+'/seg_rep1_*')[0], 'r') as rep1_guidefile:
        lines_1 = rep1_guidefile.readlines()
        for l in lines_1:
            tracks_bgs_rep1[l.split('\t')[0]] = l.split('\t')[1][:-1]
    
    tracklist_rep1 = ''
    for k, v in tracks_bgs_rep1.items():
        tracklist_rep1 = tracklist_rep1 + '-t {}={} '.format(k, v)
    
    os.system(
        'genomedata-load -s {} --sizes {} --verbose {}.genomedata'.format(
            sequence_file, tracklist_rep1, celltype_dir+'/rep1'))

    tracks_bgs_rep2 = {}
    with open(glob.glob(celltype_dir+'/seg_rep2_*')[0], 'r') as rep2_guidefile:
        lines_2 = rep2_guidefile.readlines()
        for l in lines_2:
            tracks_bgs_rep2[l.split('\t')[0]] = l.split('\t')[1][:-1]
    
    tracklist_rep2 = ''
    for k, v in tracks_bgs_rep2.items():
        tracklist_rep2 = tracklist_rep2 + '-t {}={} '.format(k, v)
    
    os.system(
        'genomedata-load -s {} --sizes {} --verbose {}.genomedata'.format(
            sequence_file, tracklist_rep2, celltype_dir+'/rep2'))

    

def segway_parameters(celltype, replicate_number, random_seed=73, param_init_test=False):
    # replicate number should be in format "repN" -> i.e. rep1, rep2 
    '''
    For a cell type with M available datasets(tracks), we ask Segway to assign 10 + 2*sqrt(M) different states(labels).
    REF: https://link.springer.com/content/pdf/10.1186/s13059-019-1784-2.pdf
    '''

    num_tracks = None ###TEMP###

    if param_init_test:
        name_sig = celltype+'_'+replicate_number+'_'+"param_init_rs{}".format(random_seed)
    else:
        name_sig = celltype+'_'+replicate_number

    params_dict = {
        "random_seed":random_seed, "include":pilot_regions_file, "track_weight":0.01,
        "stws":1, "ruler_scale":100, "prior_strength":1, "resolution":100, "mini_batch_fraction":0.5,
        "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
        "name_sig":name_sig, "genomedata_file":celltype+'/'+replicate_number, 
        "traindir":name_sig+'_train', "posteriordir":name_sig+'_posterior'
    }

    return params_dict

def run_segway_parse_results(params_dict):
    run_segway_and_post_process(params_dict)

    parse_posterior_results(
        params_dict['name_sig'], params_dict['include'], params_dict['resolution'], M=100)

    pass

def segway_param_inits(celltype, replicate_number, random_seed_list):
    '''
    train segway with different parameter initializations (random seeds) 
    to compare how similar/reproducible results are with differet param inits
    '''

    for rs in random_seed_list:
        params = segway_parameters(celltype, replicate_number, random_seed=int(rs))
        run_segway_parse_results(params)

def concat_genomedata():
    pass

def segway_concat(rep1dir, rep2dir, concat_seg_res, sizesfile):
    '''
    navigate fcoc files in rep1
    navigate fcoc files in rep2

    rename_chroms()
    virtualize_chroms()

    virtualize_sizes_file()

    make genomedata_concat
    make genomedata_rep1_renamed
    make genomedata_rep2_renamed

    define params_dict
    run concat seg + post clean up

    parse posterior1
    parse posterior2
    '''
    pass
    
def get_biological_pnemonics(results_directory, segway=True):
    '''
    - run segtools signal-distribution in the segway's output directory and 
        move the results to label_interpretation/segwayOutput/*.
    - run segtools feature-aggregation in label_interpretation/segwayOutput/*
    - move the corresponding trackname_assays.txt
    - run apply_samples.py
    - get pnemonics and move bach to segway's output directory'''
    pass

def matching_clustering(seg_results_directory_1, seg_results_directory_2, cluster=True, concat=False):
    '''
    provides options to simply match labels (without clustering)

    OR

    order-based hierarchical cluster matching.

    if concatenated:
        no matching
        just clustering
    
    '''
    pass

def reproducibility_analysis():
    pass

"""when running the whole script from start to end to generate (and reproduce) results
remember to put label interpretation in try blocks (skippable) to prevent any kind of
dependency issue of segtools to cause issues in reproducibility of results"""


if __name__=="__main__":

    CellType_list = np.array(
        ['K562', 'MCF-7', 'GM12878', 'HeLa-S3', 'CD14-positive monocyte'])

    download_dir = 'files/'

    print('list of target celltypes', CellType_list)
    existing_data = np.array(check_if_data_exists(CellType_list, download_dir))
    CellType_list = [CellType_list[i] for i in range(len(CellType_list)) if existing_data[i]==False]

    if len(CellType_list) != 0:
        download_encode_files(CellType_list, download_dir, "GRCh38")
    else:
        print('No download required!')

    CellType_list = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+ct)]

    # clean up potential space characters in directory names to prevent later issues
    for ct in CellType_list:
        if " " in ct:
            os.system("mv {} {}".format(
                ct.replace(' ', '\ '), ct.replace(" ", "_")
            ))

    CellType_list = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+ct)]
    create_trackname_assay_file(download_dir)

    assays = {}
    for ct in CellType_list:
        assays[ct] = read_list_of_assays(download_dir+ct)

    print(assays)

    # convert all bigwigs to bedgraphs (for segway)
    for k, v in assays.items():
        for t in v:
            Convert_all_BW2BG(download_dir+k+'/'+t)

    # metadata = read_metadata(download_dir)
    for c in CellType_list:
        gather_replicates(celltype_dir=download_dir+c)

    
    # download chromosome sizes file for hg38
    if os.path.exists(download_dir+"hg38.chrom.sizes") == False:
        sizes_url = 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes'
        sizes_file_dl_response = requests.get(sizes_url, allow_redirects=True)
        open(download_dir+"hg38.chrom.sizes", 'wb').write(sizes_file_dl_response.content)
        print('downloaded the hg38.chrom.sizes file')

    # check for existence of genomedata files
    gd_exists = []
    for ct in CellType_list:
        if os.path.exists(download_dir+ct+'/rep1.genomedata') == False or \
            os.path.exists(download_dir+ct+'/rep2.genomedata') == False:
            gd_exists.append(False)
        
        else:
            gd_exists.append(True)

    gd_to_create = [CellType_list[i] for i in range(len(CellType_list)) if gd_exists[i]==False]

    if len(gd_to_create) != 0:
        p_obj = mp.Pool(len(gd_to_create))
        p_obj.map(partial(
            create_genomedata, sequence_file=download_dir+"hg38.chrom.sizes"), 
            [download_dir + ct for ct in gd_to_create])