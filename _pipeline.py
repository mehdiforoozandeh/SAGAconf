from _utils import *
import os

"""
create genomefiles using the track data 
and reference genome position files
"""

def prepare_track(infile, outfile):
    '''prepare Wig or bedGraph track files for creating genomedata file'''
    pass
    

def create_genomefile(seq_file, track_file, out_dir):
    '''
    seq_file and track_file both should be list
    '''

    seq_cmd = ' '.join("-s " + sf for sf in list(seq_file))
    track_cmd = ' '.join("-t " + tf for tf in list(track_file))
    
    os.system(
        'genomedata-load {} {} {} --verbose'.format(
            seq_cmd, track_cmd, out_dir
            ))

'''
run segway train and segway posterior
'''

def run_segway_train_and_posterior(slurm_job_filename,
    genomedata_file, traindir, posteriordir, include_coords_file, 
    num_labels, resolution, track_weight, segtransition_weight, 
    cluster=False,
    base_job_file=None, memory="30G", n_cpus='5', time="0-10:00", 
    num_instances=10, submit=True):

    '''
    initiates .sh job file to be submitted using slurm's "sbatch"
    - if cluster == True AND base_job_file == None use default base for Cedar cluster
    
    '''

    if base_job_file == None:
        default_base = '''#!/bin/bash\n#SBATCH --account=def-maxwl\n'''
        
        base_args = '''#SBATCH --cpus-per-task={}\n#SBATCH --mem={}\n#SBATCH --time={}\n'''.format(
            n_cpus, memory, time
        )

        base = default_base + base_args

    else:
        with open(base_job_file, 'r') as fobj:
            base = fobj.read()


    if os.path.exists(traindir):
        os.system('rm -r {}'.format(traindir))
        os.system('mkdir {}'.format(traindir))
    
    train_cmd = '''segway train --include-coords={} --num-instances={}\
         --track-weight={} --segtransition-weight-scale={} --resolution={}\
              --num-labels={} {} {}'''.format(
                  include_coords_file, num_instances, track_weight, segtransition_weight, 
                  resolution, num_labels, genomedata_file, traindir
              )

    if os.path.exists(posteriordir):
        os.system('rm -r {}'.format(posteriordir))
        os.system('mkdir {}'.format(posteriordir))

    posterior_cmd = '''segway posterior --include-coords={} {} {} {}'''.format(
        include_coords_file, genomedata_file, traindir, posteriordir
    )
        
    final_cmd = base + '\n' + traindir + '\n' + posteriordir

    with open(slurm_job_filename, 'w') as job_fobj:
        job_fobj.write(final_cmd)

    if submit:
        os.system('sbatch {}'.format(slurm_job_filename))
    
'''
parse results into resolution sized bins

[[ NEEDS OPTIMIZATION! CURRENTLY NOT FAST ENOUGH! ]]
'''

def posterior_results_into_bins(posterior_df_list, windowsize, regions_file, num_labels, n_subset):
    position_map = define_positions(regions_file, windowsize)
    # loci is the position map, updated with posterior values (completed position map)

    loci = position_labels(posterior_df_list, position_map, windowsize, num_labels=int(num_labels), n_subset=None)
    return loci

'''
match labels using hungarian algorithms & return corrected labels
'''

def match_labels(loci_1, loci_2, num_labels):
    conf_mat = confusion_matrix(loci_1, loci_2, num_labels)
    assignment_pairs = Hungarian_algorithm(conf_mat, verbose=True)
    print("label_mapping:\t", assignment_pairs)

    loci_1, corrected_loci_2 = connect_bipartite(loci_1, loci_2, assignment_pairs)

    return loci_1, corrected_loci_2