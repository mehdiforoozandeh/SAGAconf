from _utils import *
from _pipeline import *
import os, random,sys
from sklearn.model_selection import ParameterGrid
import multiprocessing as mp    
from datetime import datetime


def grid_search(grid_space_dict):
    return list(ParameterGrid(grid_space_dict))

def batch_run_mp(grid_list, static_params, n_threads=1):
    folder_prefix = datetime.now().strftime("%Y%m%d_%H%M")
    run_instances = []

    counter = 0
    for instance in grid_list:
        counter += 1 
        run_instances.append({
            "name_sig": "batch_runs/"+folder_prefix+'/'+str(counter), 
            "traindir": "batch_runs/"+folder_prefix+"/train_"+str(counter), 
            "posteriordir": "batch_runs/"+folder_prefix+"/posterior_"+str(counter),
            **instance, **static_params})
        
    random.shuffle(run_instances)
    pool = mp.Pool(n_threads)
    pool.map(run_segway_and_post_process, run_instances)

    list_of_name_sigs = [i["name_sig"] for i in run_instances]
    return list_of_name_sigs

def parse_posterior_batch(list_of_name_sigs, include, resolution, M=50):
    for i in list_of_name_sigs:
        parse_posterior_results(i, include, resolution, M)

if __name__ == "__main__":
    # ======== param_dict sample ======== #

    # {"name_sig":None, "random_seed":None, "include":None, "track_weight":None,
    # "stws":None, "ruler_scale":None, "prior_strength":None, "resolution":None,
    # "mini_batch_fraction":None,
    # "num_labels":None, "genomedata_file":None, "traindir":None, "posteriordir":None}'''

    genomedatafile = sys.argv[1]
    include = sys.argv[2]
    # ======== param_grid_space ======== #
    grid_space = {
        "track_weight":[1, 1e-1, 1e-2, 1e-3], 
        "stws":[1e-2, 1e-1, 1, 10, 100], "ruler_scale":[100], 
        "prior_strength":[0, 1], "num_labels":[16]
        }

    # ======== static parameters ======== #
    static_params = {
        "random_seed":73, "include":include, 
        "genomedata_file":genomedatafile, "resolution":100,
        "mini_batch_fraction":0.5,
    }

    grid_list = grid_search(grid_space)
    list_of_name_sigs = batch_run_mp(grid_list, static_params, n_threads=5)
    parse_posterior_batch(list_of_name_sigs, include, static_params['resolution'], M=50)