import os
from ._utils import *

def concat_rename_chroms(in_bg_file, out_bg_file, chr_suffix): 
    # chr_suffix should be either 1 or 2
    with open(in_bg_file, 'r') as inf:
        with open(out_bg_file, 'w') as outf:
            for l in inf:
                splitted_l = l.split('\t')
                splitted_l[0] = splitted_l[0]+ "_"+ str(chr_suffix)
                outf.write("\t".join(splitted_l))

def concat_virtualize_chroms(bg_file1, bg_file2, out_bg):
    with open(out_bg, 'w') as outf:

        with open(bg_file1, 'r') as inf1:
            for l in inf1:
                outf.write(l)

        with open(bg_file2, 'r') as inf2:
            for l in inf2:
                outf.write(l)

def concat_virtualize_include_file(in_incl_file, out_incl_file, write_separate=True):
    #can be either include-coords file or chromsizes file
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

        if write_separate: 
            #also writes separate include files for later annotation or posterior calculation
            with open(out_incl_file.replace("_concatenated.sizes", "_VirtRep1_concatenated.sizes"), 'w') as outf:
                for l in lines:
                    splitted_l = l.split('\t')
                    splitted_l[0] = splitted_l[0]+ "_1"
                    outf.write("\t".join(splitted_l))

            with open(out_incl_file.replace("_concatenated.sizes", "_VirtRep2_concatenated.sizes"), 'w') as outf:
                for l in lines:
                    splitted_l = l.split('\t')
                    splitted_l[0] = splitted_l[0]+ "_2"
                    outf.write("\t".join(splitted_l))

def concat_segwayrun_and_postprocess(concat_param_dict):
    '''
    use concatenation of both replicates to train the segway model with 
    and then use that single model to annotate/get-posterior of both replicates 
    separately.

    ============== concat_param_dict EXAMPLE ==============

    {"name_sig_concat":None, "name_sig_rep1":None, "name_sig_rep2":None,
    "random_seed":None, "track_weight":None,"stws":None, "ruler_scale":None, 
    "prior_strength":None, "resolution":None,
    "num_labels":None, "traindir":None, "mini_batch_fraction":None,
    "genomedata_file_concat":None,  "genomedata_file_rep1":None, "genomedata_file_rep2":None,
    "posteriordir_rep1":None, "posteriordir_rep2":None}
    =======================================================
    '''

    for k, v in concat_param_dict.items():
        concat_param_dict[k] = str(v) 

    os.system('mkdir {}'.format(concat_param_dict['traindir']))

    os.system('SEGWAY_RAND_SEED={} segway train --max-train-rounds=100 --num-instances=10\
         --track-weight={} --segtransition-weight-scale={}\
             --ruler-scale={} --prior-strength={} --resolution={}\
                  --minibatch-fraction={} --num-labels={} {} {}'.format(
                      concat_param_dict["random_seed"], concat_param_dict["track_weight"],
                      concat_param_dict["stws"], concat_param_dict["ruler_scale"], concat_param_dict["prior_strength"], 
                      concat_param_dict["resolution"], concat_param_dict["mini_batch_fraction"],
                      concat_param_dict["num_labels"], concat_param_dict["genomedata_file_concat"], 
                      concat_param_dict["traindir"]
                  ))

    # run posterior rep1
    os.system('mkdir {}'.format(concat_param_dict['posteriordir_rep1']))
    if os.path.isfile(concat_param_dict['traindir']+'/params/params.params'):
        os.system('segway posterior {} {} {}'.format(
            concat_param_dict["genomedata_file_rep1"], 
            concat_param_dict["traindir"], concat_param_dict["posteriordir_rep1"]
        ))
    
    # run posterior rep2
    os.system('mkdir {}'.format(concat_param_dict['posteriordir_rep2']))
    if os.path.isfile(concat_param_dict['traindir']+'/params/params.params'):
        os.system('segway posterior {} {} {}'.format(
            concat_param_dict["genomedata_file_rep2"], 
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