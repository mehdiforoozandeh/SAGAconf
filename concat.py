import os
from _utils import *

def navigate_bgs(repath): 
    '''
    finds original bgs inside replicate path/directory'''
    lsr = os.listdir(repath)
    tracks_bgs = {}

    for i in lsr:
        track_ls = os.listdir(repath + '/' + i)
        if len(track_ls) > 1:
            for j in range(len(track_ls)):
                if 'concat' in track_ls[j]:
                    tracks_bgs[i] = repath + '/' + i + '/' + track_ls[j]

                elif 'bedGraph' in track_ls[j]:
                    tracks_bgs[i] = repath + '/' + i + '/' + track_ls[j]
    
    tracklist = {}
    for k, v in tracks_bgs.items():
        tracklist[k] = v

    return tracklist


def rename_chroms(in_bg_file, out_bg_file, chr_suffix): 
    # chr_suffix should be either 1 or 2
    with open(in_bg_file, 'r') as inf:
        with open(out_bg_file, 'w') as outf:
            for l in inf:
                splitted_l = l.split('\t')
                splitted_l[0] = splitted_l[0]+ "_"+ str(chr_suffix)
                outf.write("\t".join(splitted_l))

def virtualize_chroms(bg_file1, bg_file2, out_bg):
    with open(out_bg, 'w') as outf:

        with open(bg_file1, 'r') as inf1:
            for l in inf1:
                outf.write(l)

        with open(bg_file2, 'r') as inf2:
            for l in inf2:
                outf.write(l)

def virtualize_include_file(in_incl_file, out_incl_file):
    with open(in_incl_file, 'r') as inf:
        lines = inf.readlines()

        with open(out_incl_file, 'w') as outf:
            for l in lines:
                splitted_l = l.split('\t')
                splitted_l[0] = splitted_l[0]+ "_1"
                outf.write("\t".join(splitted_l))
            
            for l in lines:
                splitted_l = l.split('\t')
                splitted_l[0] = splitted_l[0]+ "_2"
                outf.write("\t".join(splitted_l))

def make_concat_genomedata(bedgraph_file_dict, seqfile, out_gd):
    tracklist = ''
    for k, v in bedgraph_file_dict.items():
        tracklist = tracklist + '-t {}={} '.format(k, v)  

    os.system('genomedata-load -s {} --sizes {} --verbose {}.genomedata'.format(seqfile, tracklist, out_gd))

def segway_concatenated_and_postprocess(concat_param_dict):
    '''
    use concatenation of both replicates to train the segway model with 
    and then use that single model to annotate/get-posterior of both replicates 
    separately.

    ============== concat_param_dict EXAMPLE ==============

    {"name_sig_concat":None, "name_sig_rep1":None, "name_sig_rep2":None,
    "random_seed":None, "include":None, "track_weight":None,
    "stws":None, "ruler_scale":None, "prior_strength":None, "resolution":None,
    "num_labels":None, "traindir":None, "mini_batch_fraction":None,
    "genomedata_file_concat":None,  "genomedata_file_rep1":None, "genomedata_file_rep2":None,
    "posteriordir_rep1":None, "posteriordir_rep2":None}
    =======================================================
    '''
    for k, v in concat_param_dict.items():
        concat_param_dict[k] = str(v) 

    os.system('mkdir {}'.format(concat_param_dict['traindir']))

    os.system('SEGWAY_RAND_SEED={} segway train --include-coords={}\
         --num-instances=10 --track-weight={} --segtransition-weight-scale={}\
             --ruler-scale={} --prior-strength={} --resolution={}\
                  --minibatch-fraction={} --num-labels={} {} {}'.format(
                      concat_param_dict["random_seed"], concat_param_dict["include"], concat_param_dict["track_weight"],
                      concat_param_dict["stws"], concat_param_dict["ruler_scale"], concat_param_dict["prior_strength"], 
                      concat_param_dict["resolution"], concat_param_dict["mini_batch_fraction"],
                      concat_param_dict["num_labels"], concat_param_dict["genomedata_file_concat"], 
                      concat_param_dict["traindir"]
                  ))

    # run posterior rep1
    os.system('mkdir {}'.format(concat_param_dict['posteriordir_rep1']))
    if os.path.isfile(concat_param_dict['traindir']+'/params/params.params'):
        os.system('segway posterior --include-coords={} {} {} {}'.format(
            concat_param_dict["include"], concat_param_dict["genomedata_file_rep1"], 
            concat_param_dict["traindir"], concat_param_dict["posteriordir_rep1"]
        ))
    
    # run posterior rep2
    os.system('mkdir {}'.format(concat_param_dict['posteriordir_rep2']))
    if os.path.isfile(concat_param_dict['traindir']+'/params/params.params'):
        os.system('segway posterior --include-coords={} {} {} {}'.format(
            concat_param_dict["include"], concat_param_dict["genomedata_file_rep2"], 
            concat_param_dict["traindir"], concat_param_dict["posteriordir_rep2"]
        ))

    # clean-up
    if os.path.isfile(concat_param_dict['posteriordir_rep1']+'/segway.bed.gz') and os.path.isfile(concat_param_dict['posteriordir_rep2']+'/segway.bed.gz'):

        gather_results(concat_param_dict["traindir"], concat_param_dict['posteriordir_rep1'], concat_param_dict["name_sig_rep1"])
        run_segtools(concat_param_dict["name_sig_rep1"])

        gather_results(concat_param_dict["traindir"], concat_param_dict['posteriordir_rep2'], concat_param_dict["name_sig_rep2"])
        run_segtools(concat_param_dict["name_sig_rep2"])

        clean_up(concat_param_dict["traindir"], concat_param_dict['posteriordir_rep1'])
        clean_up(concat_param_dict["traindir"], concat_param_dict['posteriordir_rep2'])

        with open(concat_param_dict["name_sig_rep1"]+"/run_parameters", 'w') as write_params:
            write_params.write(str(concat_param_dict))
        
        with open(concat_param_dict["name_sig_rep2"]+"/run_parameters", 'w') as write_params:
            write_params.write(str(concat_param_dict))

        # QC rep1
        posterior_df_list = read_posteriordir(concat_param_dict["name_sig_rep1"])
        print('{}: loaded posterior bedgraphs'.format(concat_param_dict['name_sig_rep1']))
        write_qc = open('{}/QC.txt'.format(concat_param_dict['name_sig_rep1']), 'w') 

        try:
            posterior_quality = QC(posterior_df_list)
            print('{}: calculated QC'.format(concat_param_dict['name_sig_rep1']))
            write_qc.write('posterior_QC(length dependent): {}\n'.format(posterior_quality))
        except:
            write_qc.write('FAILED AT CALCULATING QC score')

        try:
            number_of_segments = sum(1 for line in open(concat_param_dict["name_sig_rep1"]+'/segway.bed', 'r'))
            print('{}: counted # of segments'.format(concat_param_dict['name_sig_rep1']))
            write_qc.write('number_of_segments: {}\n'.format(number_of_segments))
        except:
            write_qc.write('FAILED AT COUNTING NUMBER OF SEGMENTS')
            
        write_qc.close()
        print('wrote the results into {}/QC.txt'.format(concat_param_dict['name_sig_rep1']))

        # QC rep2
        posterior_df_list = read_posteriordir(concat_param_dict["name_sig_rep2"])
        print('{}: loaded posterior bedgraphs'.format(concat_param_dict['name_sig_rep2']))
        write_qc = open('{}/QC.txt'.format(concat_param_dict['name_sig_rep2']), 'w') 

        try:
            posterior_quality = QC(posterior_df_list)
            print('{}: calculated QC'.format(concat_param_dict['name_sig_rep2']))
            write_qc.write('posterior_QC(length dependent): {}\n'.format(posterior_quality))
        except:
            write_qc.write('FAILED AT CALCULATING QC score')

        try:
            number_of_segments = sum(1 for line in open(concat_param_dict["name_sig_rep2"]+'/segway.bed', 'r'))
            print('{}: counted # of segments'.format(concat_param_dict['name_sig_rep2']))
            write_qc.write('number_of_segments: {}\n'.format(number_of_segments))
        except:
            write_qc.write('FAILED AT COUNTING NUMBER OF SEGMENTS')
            
        write_qc.close()
        print('wrote the results into {}/QC.txt'.format(concat_param_dict['name_sig_rep2']))

    else:
        print("FAILED AT RUNNING SEGWAY!")
        clean_up(concat_param_dict["traindir"], concat_param_dict['posteriordir_rep1'])
        clean_up(concat_param_dict["traindir"], concat_param_dict['posteriordir_rep2'])


def main(rep1dir, rep2dir, concat_dir, include_file, seqfile):
    # search both replicate folders for bedgraph files to edit (make 2 dictionaries)
    rep1_tracks = navigate_bgs(rep1dir)
    rep2_tracks = navigate_bgs(rep2dir)

    # update all intended bedgraphs using rename_chroms()
    for k, v in rep1_tracks.items():
        if 'concat' not in v:
            rename_chroms(v, v.replace('.bedGraph', '_concat.bedGraph'), "1")
            rep1_tracks[k] = v.replace('.bedGraph', '_concat.bedGraph')

    for k, v in rep2_tracks.items():
        if 'concat' not in v:
            rename_chroms(v, v.replace('.bedGraph', '_concat.bedGraph'), "2")
            rep2_tracks[k] = v.replace('.bedGraph', '_concat.bedGraph')

    # virtualize bedgraphs into concatenated tracks    
    # (this means that we will have 1 bg for each track instead of initial 2)

    assert list(rep1_tracks.keys()) == list(rep2_tracks.keys())

    concat_tracks = {}
    for k in rep1_tracks.keys():
        #create concat_dir and subdirs
        if not os.path.exists(concat_dir):
            os.mkdir(concat_dir)
        
        if not os.path.exists(concat_dir+'/'+str(k)):
            os.mkdir(concat_dir+'/'+str(k))
            
        virtualize_chroms(
            rep1_tracks[k], rep2_tracks[k], concat_dir+'/'+str(k)+'/'+str(k)+'_concatenated.bedGraph')

        concat_tracks[k] = concat_dir+'/'+str(k)+'/'+str(k)+'_concatenated.bedGraph'
        

    # update/virtualize include file
        virtualize_include_file(
            include_file, include_file.replace('.', '_concat.'))

    # make 3 genomedata files for 1)concat 2)rep1 3)rep2 (all using updated bgs)
    make_concat_genomedata(concat_tracks, seqfile, concat_dir+'/concat.gd')
    make_concat_genomedata(rep1_tracks, seqfile, concat_dir+'/rep1.gd')
    make_concat_genomedata(rep2_tracks, seqfile, concat_dir+'/rep2.gd')
    
    # create concat_param_dict
    concat_param_dict = {
        "name_sig_concat":"seg_concat", "name_sig_rep1":"seg_rep1", "name_sig_rep2":"seg_rep2",
        "random_seed":73, "include":include_file, "track_weight":0.01,
        "stws":10, "ruler_scale":100, "prior_strength":0, "resolution":100,
        "num_labels":16, "traindir":"seg_concat_train", "mini_batch_fraction":0.2,
        "genomedata_file_concat":concat_dir+'/concat.gd',  "genomedata_file_rep1":concat_dir+'/rep1.gd', 
        "genomedata_file_rep2":concat_dir+'/rep2.gd',
        "posteriordir_rep1":'seg_rep1_posterior', "posteriordir_rep2":'seg_rep2_posterior'
        }

    segway_concatenated_and_postprocess(concat_param_dict)

if __name__=="__main__":
    main(
        'rep1', 'rep2', 'concatenated_files', 
        'encodePilotRegions.hg19.bed', 'hg38.chrom.sizes')