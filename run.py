from operator import truediv
from os import lseek
from types import CellType
from unicodedata import name
from _utils import *
from _pipeline import *
from _reproducibility import *
from _interpret import *
from _visualize import *
from get_data import *
from matchlabels import *
import pandas as pd

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
    for f in ls:
        if ".bigWig" in f:
            if f.replace(".bigWig",".bedGraph") not in ls:
                if " " in f:
                    f = f.replace(" ","\ ")
                print(
                    "./bigWigToBedGraph {}.bigWig {}.bedGraph".format(
                        celltype_dir+'/'+f.replace(".bigWig",""), 
                        celltype_dir+'/'+f.replace(".bigWig","")))
                os.system(
                    "./bigWigToBedGraph {}.bigWig {}.bedGraph".format(
                        celltype_dir+'/'+f.replace(".bigWig",""), 
                        celltype_dir+'/'+f.replace(".bigWig","")))

def gather_segway_replicates(celltype_dir):
    pass

def create_genomedata(celltype_dir, sequence_file):
    create_genomedata_rep(celltype_dir + '/rep1', sequence_file)
    create_genomedata_rep(celltype_dir + '/rep2', sequence_file)
    pass

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
        ['K562', 'MCF-7', 'GM12878', 'HeLa-S3', 'HepG2', 'CD14-positive monocyte']
        )

    download_dir = 'files/'

    print('list of target celltypes', CellType_list)
    existing_data = np.array(check_if_data_exists(CellType_list, download_dir))
    CellType_list = np.delete(CellType_list, np.where(CellType_list[existing_data==True]))

    if len(CellType_list) != 0:
        download_encode_files(CellType_list, download_dir, "GRCh38")
    else:
        print('No download required!')

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



    