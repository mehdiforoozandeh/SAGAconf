from _utils import *
from _pipeline import *
import os, random,sys
from sklearn.model_selection import ParameterGrid
import multiprocessing as mp    
from datetime import datetime


# ======== param_dict sample ======== #

genome_data_files = ['batch_runs/rep1.genomedata', 'batch_runs/rep2.genomedata']

paramdict_1 = {"name_sig_main":"batch_runs/concat_11", "name_sig_aux":"batch_runs/concat_12",
    "random_seed":73, "include":'batch_runs/encodePilotRegions.hg19.bed', "track_weight":float(1e-2/12),
    "stws":1e-2, "ruler_scale":1000, "prior_strength":1, "resolution":1000,
    "num_labels":16, "traindir":'batch_runs/train_dir_1', "mini_batch_fraction":0.1,
    "genomedata_file_main":genome_data_files[0],  "genomedata_file_aux":genome_data_files[1],
    "posteriordir_main":"batch_runs/posterior_concat_11",  "posteriordir_aux":"batch_runs/posterior_concat_12"}

paramdict_2 = {"name_sig_main":"batch_runs/concat_21", "name_sig_aux":"batch_runs/concat_22",
    "random_seed":73, "include":'batch_runs/encodePilotRegions.hg19.bed', "track_weight":float(1e-2/12),
    "stws":1e-2, "ruler_scale":1000, "prior_strength":1, "resolution":1000,
    "num_labels":16, "traindir":'batch_runs/train_dir_2', "mini_batch_fraction":0.1,
    "genomedata_file_main":genome_data_files[1],  "genomedata_file_aux":genome_data_files[0],
    "posteriordir_main":"batch_runs/posterior_concat_21",  "posteriordir_aux":"batch_runs/posterior_concat_22"}

segway_concatenated_and_postprocess(paramdict_1)
segway_concatenated_and_postprocess(paramdict_2)

parse_posterior_results("batch_runs/concat_11", 'batch_runs/encodePilotRegions.hg19.bed', 1000, M=50)
parse_posterior_results("batch_runs/concat_12", 'batch_runs/encodePilotRegions.hg19.bed', 1000, M=50)
parse_posterior_results("batch_runs/concat_21", 'batch_runs/encodePilotRegions.hg19.bed', 1000, M=50)
parse_posterior_results("batch_runs/concat_22", 'batch_runs/encodePilotRegions.hg19.bed', 1000, M=50)