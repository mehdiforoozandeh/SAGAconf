import os
from run import *

def run_clean_parse(param_dict):
    if os.path.exists(param_dict["traindir"]) == False:
        os.system("segway train  --num-instances=10  --resolution={} --minibatch-fraction=0.03 --num-labels=16 {} {}".format(param_dict["resolution"], param_dict["genomedata_file"], param_dict["traindir"]))
    
    if os.path.exists(param_dict["posteriordir"]) == False:
        os.system('segway posterior {} {} {}'.format(param_dict["genomedata_file"], param_dict["traindir"], param_dict["posteriordir"]))

    gather_results(param_dict["traindir"], param_dict['posteriordir'], param_dict["name_sig"])
    run_segtools(param_dict["name_sig"])

    if os.path.exists(param_dict["traindir"]):
        shutil.rmtree(param_dict["traindir"])
    
    if os.path.exists(param_dict["posteriordir"]):
        shutil.rmtree(param_dict["posteriordir"])

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

    parse_posterior_results(param_dict["name_sig"], param_dict["resolution"], mp=False)

if __name__=="__main__":
    param_dict = {"name_sig":"testsegway/res100_GM12878_1", "resolution":100, "genomedata_file":"_protect_files_/GM12878/rep1.genomedata", "traindir":"testsegway/res100_GM12878_1_train", "posteriordir":"testsegway/res100_GM12878_1_posterior"}
    run_clean_parse(param_dict)

    param_dict = {"name_sig":"testsegway/res100_GM12878_2", "resolution":100, "genomedata_file":"_protect_files_/GM12878/rep2.genomedata", "traindir":"testsegway/res100_GM12878_2_train", "posteriordir":"testsegway/res100_GM12878_2_posterior"}
    run_clean_parse(param_dict)

    param_dict = {"name_sig":"testsegway/res200_GM12878_1", "resolution":200, "genomedata_file":"_protect_files_/GM12878/rep1.genomedata", "traindir":"testsegway/res200_GM12878_1_train", "posteriordir":"testsegway/res200_GM12878_1_posterior"}
    run_clean_parse(param_dict)

    param_dict = {"name_sig":"testsegway/res200_GM12878_2", "resolution":200, "genomedata_file":"_protect_files_/GM12878/rep2.genomedata", "traindir":"testsegway/res200_GM12878_2_train", "posteriordir":"testsegway/res200_GM12878_2_posterior"}
    run_clean_parse(param_dict)