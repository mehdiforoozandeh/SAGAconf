from ._utils import *
from ._cluster_matching import *

import os


"""
create genomefiles using the track data 
and reference genome position files
"""

def bigwig_to_bedgraph_rep(repath):
    '''
    provide a directroy containing sub-directories named after the track-name.
    these subdirectories must contain bigWig files of the track.
    the bigWig file will be converted to bedGraph
    '''
    lsr = os.listdir(repath)

    for i in lsr:    

        track_ls = os.listdir(repath + '/' + i)
        if len(track_ls) > 1:
            for j in range(len(track_ls)):
                if 'bigWig' in track_ls[j]:
                    bw_path = repath + '/' + i + '/' + track_ls[j]
        else: 
            bw_path = repath + '/' + i + '/' + track_ls[0]

        accname = bw_path[bw_path.find('ENC'):bw_path.find('.')]
        os.system('../bigWigToBedGraph {} {}'.format(bw_path, bw_path.replace("bigWig", "bedGraph")))
        print("made BedGraph for {}".format(i))

def create_genomedata_rep(repath, seqfile):
    '''
    provide a directroy containing sub-directories named after the track-name.
    these subdirectories must contain bedGraph files of the track.
    this function integrate all of these data and loads them on a genomedata file
    '''
    lsr = os.listdir(repath)
    tracks_bgs = {}

    for i in lsr:
        track_ls = os.listdir(repath + '/' + i)
        if len(track_ls) > 1:
            for j in range(len(track_ls)):
                if 'bedGraph' in track_ls[j]:
                    tracks_bgs[i] = repath + '/' + i + '/' + track_ls[j]
    
    tracklist = ''
    for k, v in tracks_bgs.items():
        tracklist = tracklist + '-t {}={} '.format(k, v)
    
    os.system('genomedata-load -s {} --sizes {} --verbose {}.genomedata'.format(seqfile, tracklist, repath))


def run_segway_and_post_process(param_dict):
    '''
    run segway train and segway posterior

    ======== param_dict sample ======== #

    {"name_sig":None, "random_seed":None, "track_weight":None,
    "stws":None, "ruler_scale":None, "prior_strength":None, "resolution":None,
    "mini_batch_fraction":None, "num_labels":None, "genomedata_file":None, 
    "traindir":None, "posteriordir":None}'''

    for k, v in param_dict.items():
        param_dict[k] = str(v) 

    os.system('mkdir {}'.format(param_dict['traindir']))

    os.system('SEGWAY_RAND_SEED={} segway train --max-train-rounds=100 --num-instances=10\
          --track-weight={} --segtransition-weight-scale={}\
             --ruler-scale={} --prior-strength={} --resolution={}\
                  --minibatch-fraction={} --num-labels={} {} {}'.format(
                      param_dict["random_seed"], param_dict["track_weight"],
                      param_dict["stws"], param_dict["ruler_scale"], param_dict["prior_strength"], 
                      param_dict["resolution"], param_dict["mini_batch_fraction"],
                      param_dict["num_labels"], param_dict["genomedata_file"], 
                      param_dict["traindir"]
                  ))

    os.system('mkdir {}'.format(param_dict['posteriordir']))
    if os.path.isfile(param_dict['traindir']+'/params/params.params'):
        os.system('segway posterior {} {} {}'.format(
            param_dict["genomedata_file"], param_dict["traindir"], 
            param_dict["posteriordir"]
        ))

    if os.path.isfile(param_dict['posteriordir']+'/segway.bed.gz'):
        gather_results(param_dict["traindir"], param_dict['posteriordir'], param_dict["name_sig"])
        run_segtools(param_dict["name_sig"])
        clean_up(param_dict["traindir"], param_dict['posteriordir'])

        with open(param_dict["name_sig"]+"/run_parameters", 'w') as write_params:
            write_params.write(str(param_dict))

        posterior_df_list = read_posteriordir(param_dict["name_sig"])
        print('{}: loaded posterior bedgraphs'.format(param_dict['name_sig']))

        write_qc = open('{}/QC.txt'.format(param_dict['name_sig']), 'w') 

        try:
            posterior_quality = QC(posterior_df_list)
            print('{}: calculated QC'.format(param_dict['name_sig']))
            write_qc.write('posterior_QC(length dependent): {}\n'.format(posterior_quality))
        except:
            write_qc.write('FAILED AT CALCULATING QC score')

        try:
            number_of_segments = sum(1 for line in open(param_dict["name_sig"]+'/segway.bed', 'r'))
            print('{}: counted # of segments'.format(param_dict['name_sig']))
            write_qc.write('number_of_segments: {}\n'.format(number_of_segments))
        except:
            write_qc.write('FAILED AT COUNTING NUMBER OF SEGMENTS')
            
        write_qc.close()

        print('wrote the posterior granularity QC results into {}/QC.txt'.format(param_dict['name_sig']))

    else:
        print("FAILED AT RUNNING SEGWAY!")
        clean_up(param_dict["traindir"], param_dict['posteriordir'])


def parse_posterior_results_old(posterior_dir, chrom_sizes_file, resolution, M):
    '''DEPRICATED'''
    posterior_df_list = read_posteriordir(posterior_dir)
    print(posterior_df_list.keys())

    coords = read_chrom_sizes_file(chrom_sizes_file)
    empty_bins = initialize_bins_sizesfile(coords, resolution)
    parsed_df = mp_binning(posterior_df_list, empty_bins, M)

    parsed_df.to_csv(posterior_dir+'/parsed_posterior.csv')
    return parsed_df

def parse_posterior_results(posterior_dir, resolution, mp=True):
    '''
    parse results into resolution sized bins
    '''
    print('Parsing posterior results for {}'.format(posterior_dir))
    parsed_df = mp_inplace_binning(
        posterior_dir, resolution, assert_coord_match=True, mp=mp)

    print('saving results to parsed_posterior.csv file')
    parsed_df.to_csv(posterior_dir+'/parsed_posterior.csv')
    return parsed_df

'''
match labels using hungarian algorithms & return corrected labels
'''

# def Hungarian_label_matching(loci_1, loci_2, num_labels):
#     conf_mat = confusion_matrix(loci_1, loci_2, num_labels)
#     assignment_pairs = Hungarian_algorithm(conf_mat, verbose=True)
#     print("label_mapping:\t", assignment_pairs)
#     loci_1, corrected_loci_2 = connect_bipartite(loci_1, loci_2, assignment_pairs)

#     return loci_1, corrected_loci_2