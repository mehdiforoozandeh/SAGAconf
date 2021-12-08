from _cluster_matching import *
from _pipeline import *
from _visualize import visualize
from matchlabels import *

def concat_label_matching(replicate_1_dir, replicate_2_dir, num_labels, num_clusters=None):

    if num_clusters != None:
        if int(num_clusters) <= num_labels: 
            print('post clst => conf_mat // number of clusters = {}'.format(num_clusters))
            cluster_matching(\
                replicate_1_dir, replicate_2_dir, n_clust=num_clusters ,matching_strategy='conf_matrix')

    else:
        print('plain seg => conf_mat // number of labels = {}'.format(num_labels))
        plain_seg_matching(\
            replicate_1_dir, replicate_2_dir, matching_strategy='conf_matrix')

def concat_reproducibility():
    pass

if __name__=="__main__":
    concat_label_matching("batch_runs/concat_11", "batch_runs/concat_12", 16, num_clusters=None)
    concat_label_matching("batch_runs/concat_21", "batch_runs/concat_22", 16, num_clusters=None)

    concat_label_matching("batch_runs/concat_11", "batch_runs/concat_12", 16, num_clusters=9)
    concat_label_matching("batch_runs/concat_21", "batch_runs/concat_22", 16, num_clusters=9)