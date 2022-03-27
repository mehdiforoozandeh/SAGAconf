import os
import pandas as pd


def binarize_data(inputbeddir, cellmarkfiletable, outputdir, resolution=100, chromlength='ChromHMM/CHROMSIZES/hg19.txt'):
    cmdline = "java -Xmx10g -jar ChromHMM/ChromHMM.jar BinarizeBam -b {} {} {} {} {}".format(
        resolution, chromlength, inputbeddir, cellmarkfiletable, outputdir
    )
    os.system(cmdline)

def learnModel(binary_input_dir, output_dir, num_labels='16', assembly='hg19', n_threads='0', random_seed=None):
    if random_seed != None:
        learnmodel_cmdline = "java -Xmx10g -jar ChromHMM/ChromHMM.jar LearnModel -init random -s {} -printposterior -p {} {} {} {} {}".format(
            random_seed, n_threads, binary_input_dir, output_dir, num_labels, assembly)
    else:
        learnmodel_cmdline = "java -Xmx10g -jar ChromHMM/ChromHMM.jar LearnModel -init information -printposterior -p {} {} {} {} {}".format(
            n_threads, binary_input_dir, output_dir, num_labels, assembly)
    os.system(learnmodel_cmdline)

def chrhmm_initialize_bin(chrom, numbins, res):
    empty_bins = []
    next_start = 0
    for _ in range(numbins):
        empty_bins.append(
            [chrom, next_start, int(next_start+res)])
        next_start = int(next_start+res)
    empty_bins = pd.DataFrame(empty_bins, columns=['chr', 'start', 'end'])
    empty_bins['start'] = empty_bins['start'].astype("int32")
    empty_bins['end'] = empty_bins['end'].astype("int32")
    return empty_bins

def read_posterior_file(filepath):
    with open(filepath,'r') as posteriorfile:
        lines = posteriorfile.readlines()
    vals = []
    for il in range(len(lines)):
        ilth_vals = lines[il].split('\t')
        ilth_vals[-1] = ilth_vals[-1].replace("\n","")
        vals.append(ilth_vals)
    vals = pd.DataFrame(vals[2:], columns=["posterior{}".format(i.replace("E","")) for i in vals[1]])
    vals = vals.astype("float32")
    return vals

def ChrHMM_read_posteriordir(posteriordir, rep, resolution=100):
    '''
    for each file in posteriordir
    Initialize emptybins based on chromsizes
    fill in the posterior values for each slot
    return DF
    '''
    ls = os.listdir(posteriordir)
    to_parse = []
    for f in ls:
        if rep in f:
            to_parse.append(f)
    parsed_posteriors = {}
    for f in to_parse:
        fileinfo = f.split("_")
        posteriors = read_posterior_file(posteriordir + '/' + f)
        bins = chrhmm_initialize_bin(fileinfo[-2], len(posteriors), resolution)
        posteriors = pd.concat([bins, posteriors], axis=1)
        parsed_posteriors[fileinfo[-2]] = posteriors 
        
    parsed_posteriors = pd.concat([parsed_posteriors[c] for c in sorted(list(parsed_posteriors.keys()))], axis=0)
    parsed_posteriors = parsed_posteriors.reset_index(drop=True)
    return parsed_posteriors

"""
1. relocate (mv/cp) .bam files of each cell type to a new directory called "chmmfiles"
2. create 3 cmft files at each chmmfiles dir :
    a. concat
    b. rep1
    c. rep2
3. run three instances of learnModel for each celltype
    a. concat
    b. rep1
    c. rep2
"""

def prepare_chmm_inputdata(CellType_dir, assertion=False):
    celltype_name = CellType_dir.split("/")[-1]

    if "chmmfiles" not in os.listdir():
        os.mkdir("chmmfiles/")

    if os.path.exists("chmmfiles/{}".format(celltype_name)) == False:
        os.mkdir("chmmfiles/{}".format(celltype_name))

    with open(CellType_dir+'/replicate_number.txt', 'r') as repguide_file:
        lines = repguide_file.readlines()
        rep1_biosampleID = lines[0][5:-1]
        rep2_biosampleID = lines[1][5:]
        if "\n" in rep2_biosampleID:
            rep2_biosampleID.replace("\n","")

    assaylist = [tr for tr in os.listdir(CellType_dir) if os.path.isdir(CellType_dir+'/'+tr)]
    navigate = []
    for tr in assaylist:
        tfmd = pd.read_csv(CellType_dir+'/'+tr+'/track_files_metadata.csv')
        tfmd.index = list(tfmd['Unnamed: 0'])
        tfmd = tfmd.drop('Unnamed: 0', axis=1) 

        if assertion:
            assert str(tfmd.loc['biosample', 'rep1_alig']) == rep1_biosampleID
            assert str(tfmd.loc['biosample', 'rep2_alig']) == rep2_biosampleID
            assert tfmd.loc['assay', 'rep1_alig'] == tr
            assert tfmd.loc['assay', 'rep2_alig'] == tr

        navigate.append(
            [tfmd.loc['assay', 'rep1_alig'], "rep1", str(tfmd.loc['accession', 'rep1_alig'])+".bam"])  
        navigate.append(
            [tfmd.loc['assay', 'rep2_alig'], "rep2", str(tfmd.loc['accession', 'rep2_alig'])+".bam"])

    
    cmft_concat = open("chmmfiles/{}/cmft_concat.txt".format(celltype_name), "w")
    cmft_rep1 = open("chmmfiles/{}/cmft_rep1.txt".format(celltype_name), "w")
    cmft_rep2 = open("chmmfiles/{}/cmft_rep2.txt".format(celltype_name), "w")
    
    for ins in navigate:
        os.system("cp {} {}".format(
            CellType_dir + "/" + ins[0] + "/" + ins[2], "chmmfiles/{}/".format(celltype_name)))    
        cmft_concat.write("{}_{}\t{}\t{}\n".format(
            celltype_name, ins[1], ins[0], ins[2]))

        if ins[1] == "rep1":
            cmft_rep1.write("{}_{}\t{}\t{}\n".format(
                celltype_name, ins[1], ins[0], ins[2]))
        elif ins[1] == "rep2":
            cmft_rep2.write("{}_{}\t{}\t{}\n".format(
                celltype_name, ins[1], ins[0], ins[2]))

def ChromHMM_replicate_runs(chmm_celltype_dir, chmm_output_dir, n_thread='0'):
    namesig = chmm_celltype_dir.split("/")[-1]

    if os.path.exists(chmm_celltype_dir+"/binarized_rep1") == False:
        binarize_data(
            chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep1.txt", chmm_celltype_dir+"/binarized_rep1", 
            resolution=100)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1") == False:
        learnModel(
            chmm_celltype_dir+"/binarized_rep1", chmm_output_dir+"/"+namesig+"_rep1", 
            num_labels=16, assembly='hg19', n_threads=n_thread, random_seed=None)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1/parsed_posterior.csv") == False:
        parsed_posterior = ChrHMM_read_posteriordir(
            chmm_output_dir+"/"+namesig+"_rep1/POSTERIOR", "rep1", resolution=100)
        
        parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep1/parsed_posterior.csv")

    ####=================================================================================================####
    if os.path.exists(chmm_celltype_dir+"/binarized_rep2") == False:
        binarize_data(
            chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep2.txt", chmm_celltype_dir+"/binarized_rep2", 
            resolution=100)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2") == False:
        learnModel(
            chmm_celltype_dir+"/binarized_rep2", chmm_output_dir+"/"+namesig+"_rep2", 
            num_labels=16, assembly='hg19', n_threads=n_thread, random_seed=None)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2/parsed_posterior.csv") == False:
        parsed_posterior = ChrHMM_read_posteriordir(
            chmm_output_dir+"/"+namesig+"_rep2/POSTERIOR", "rep2", resolution=100)
        
        parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep2/parsed_posterior.csv")


def ChromHMM_paraminit_runs(chmm_celltype_dir, chmm_output_dir, random_seeds,  n_thread='0'):
    namesig = chmm_celltype_dir.split("/")[-1]
    if "/" in namesig:
        namesig.replace("/", "")

    for rs in random_seeds:
        if os.path.exists(chmm_celltype_dir+"/binarized_rep1") == False:
            binarize_data(
                chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep1.txt", chmm_celltype_dir+"/binarized_rep1", 
                resolution=100)
        
        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1_rs{}".format(rs)) == False:
            learnModel(
                chmm_celltype_dir+"/binarized_rep1", chmm_output_dir+"/"+namesig+"_rep1_rs{}".format(rs), 
                num_labels=16, assembly='hg19', n_threads=n_thread, random_seed=int(rs))

        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1_rs{}/parsed_posterior.csv".format(rs)) == False:
            parsed_posterior = ChrHMM_read_posteriordir(
                chmm_output_dir+"/"+namesig+"_rep1_rs{}/POSTERIOR".format(rs), "rep1", resolution=100)
            parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep1_rs{}/parsed_posterior.csv".format(rs))
        
        ####=================================================================================================####

        if os.path.exists(chmm_celltype_dir+"/binarized_rep2") == False:
            binarize_data(
                chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep2.txt", chmm_celltype_dir+"/binarized_rep2", 
                resolution=100)
        
        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2_rs{}".format(rs)) == False:
            learnModel(
                chmm_celltype_dir+"/binarized_rep2", chmm_output_dir+"/"+namesig+"_rep2_rs{}".format(rs), 
                num_labels=16, assembly='hg19', n_threads=n_thread, random_seed=int(rs))

        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2_rs{}/parsed_posterior.csv".format(rs)) == False:
            parsed_posterior = ChrHMM_read_posteriordir(
                chmm_output_dir+"/"+namesig+"_rep2_rs{}/POSTERIOR".format(rs), "rep2", resolution=100)
            parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep2_rs{}/parsed_posterior.csv".format(rs))

def ChromHMM_concat_runs(chmm_celltype_dir, chmm_output_dir, n_thread='0'):
    namesig = chmm_celltype_dir.split("/")[-1]
    if "/" in namesig:
        namesig.replace("/", "")
    
    if os.path.exists(chmm_celltype_dir+"/binarized_concat") == False:
        binarize_data(
            chmm_celltype_dir, chmm_celltype_dir+"/cmft_concat.txt", chmm_celltype_dir+"/binarized_concat", 
            resolution=100)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_concat") == False:
        learnModel(
            chmm_celltype_dir+"/binarized_concat", chmm_output_dir+"/"+namesig+"_concat", 
            num_labels=16, assembly='hg19', n_threads=n_thread, random_seed=None)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_concat/parsed_posterior_rep1.csv") == False:
        parsed_posterior = ChrHMM_read_posteriordir(
            chmm_output_dir+"/"+namesig+"_concat/POSTERIOR", "rep1", resolution=100)
        
        parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_concat/parsed_posterior_rep1.csv")
    
    if os.path.exists(chmm_output_dir+"/"+namesig+"_concat/parsed_posterior_rep2.csv") == False:
        parsed_posterior = ChrHMM_read_posteriordir(
            chmm_output_dir+"/"+namesig+"_concat/POSTERIOR", "rep2", resolution=100)

        parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_concat/parsed_posterior_rep2.csv")

