from .._pipeline import run_segway_and_post_process
import multiprocessing as mp
import sys

def main(n_threads, rep):
    param_dicts = []
    
    static_params = {
        "random_seed":73, "include":'encodePilotRegions.hg19.bed', "track_weight":1/(100*12),
        "stws":1, "ruler_scale":100, "prior_strength":1, 
        "resolution":100, "genomedata_file":'{}.genomedata'.format(rep)
    }

    non_static_params ={
        'name_sig':None, 'num_labels':None,
        "traindir":None, "posteriordir":None
    }
    
    for i in range(8,17): #want to test num_labels in range [8,16]
        non_static_params['name_sig'] = 'rep2_{}labels'.format(i)
        non_static_params['traindir'] = 'train_rep2_{}labels'.format(i)
        non_static_params['posteriordir'] = 'posterior_rep2_{}labels'.format(i)
        non_static_params['num_labels'] = int(i)

        param_dicts.append({**static_params, **non_static_params})
        
    pool = mp.Pool(n_threads)
    pool.map(run_segway_and_post_process, param_dicts)

    
if __name__=="__main__":
    if sys.argv[1] == 'rep1':
        main(n_threads=3, rep='rep1')
        
    elif sys.argv[1] == 'rep2':
        main(n_threads=3, rep='rep2')