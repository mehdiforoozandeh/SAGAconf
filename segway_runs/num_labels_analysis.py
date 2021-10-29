import multiprocessing as mp
import sys
sys.path.insert(0, '..')
import _pipeline

def main(n_threads, rep, res, n_labels_list):
    param_dicts = []

    static_params = {
        "random_seed":73, "include":'encodePilotRegions.hg19.bed', "track_weight":1/((res/10)*12),
        "stws":0.01, "ruler_scale":res, "prior_strength":1,
        "resolution":res, "genomedata_file":'{}.genomedata'.format(rep)
    }

    non_static_params ={
        'name_sig':None, 'num_labels':None,
        "traindir":None, "posteriordir":None
    }

    for i in n_labels_list:
        non_static_params['name_sig'] = '{}_{}labels_res{}'.format(rep, i, res)
        non_static_params['traindir'] = 'train_{}_{}labels_res{}'.format(rep, i, res)
        non_static_params['posteriordir'] = 'posterior_{}_{}labels_res{}'.format(rep, i, res)
        non_static_params['num_labels'] = int(i)

        param_dicts.append({**static_params, **non_static_params})

    pool = mp.Pool(n_threads)
    pool.map(_pipeline.run_segway_and_post_process, param_dicts)


if __name__=="__main__":
    if sys.argv[1] == 'rep1':
        main(n_threads=3, rep='rep1', res=1000, n_labels_list=[13, 14])
        main(n_threads=3, rep='rep1', res=100, n_labels_list=list(range(8,17)))

    elif sys.argv[1] == 'rep2':
        main(n_threads=3, rep='rep2', res=1000, n_labels_list=[10, 13, 14, 15, 16])
        main(n_threads=3, rep='rep1', res=100, n_labels_list=list(range(8,17)))