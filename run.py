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
    cell_track_nav = {}

    for celltype in celltype_list:
        track_nav = search_encode(celltype, download_dir, target_assembly=target_assembly)
        cell_track_nav[str(celltype)] = track_nav
    
    return cell_track_nav
        
def gather_segway_replicates(celltype_dir, create_trackname_assay=True):
    if create_trackname_assay:
        pass
    pass

def create_genomedata(celltype_dir, sequence_file):
    create_genomedata_rep(celltype_dir + '/rep1', sequence_file)
    create_genomedata_rep(celltype_dir + '/rep2', sequence_file)
    pass

def segway_parameters(celltype, replicate_number, random_seed=73):
    # replicate number should be in format "repN" -> i.e. rep1, rep2 
    '''
    For a cell type with M available datasets(tracks), we ask Segway to assign 10 + 2*sqrt(M) different states(labels).
    REF: https://link.springer.com/content/pdf/10.1186/s13059-019-1784-2.pdf
    '''

    num_tracks = None ###TEMP###
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
    # parse_posterior_results 1
    # parse_posterior_results 2
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


