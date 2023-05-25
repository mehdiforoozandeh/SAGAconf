import os
import pandas as pd
import numpy as np

def binarize_data(inputbeddir, cellmarkfiletable, outputdir, resolution=200, chromlength='ChromHMM/CHROMSIZES/hg38.txt'):
    cmdline = "java -Xmx20g -jar ChromHMM/ChromHMM.jar BinarizeBed -b {} -t {} {} {} {} {}".format(
        resolution, outputdir+"/signals", chromlength, inputbeddir, cellmarkfiletable, outputdir
    )
    os.system(cmdline)

def learnModel(binary_input_dir, output_dir, num_labels='16', assembly='hg38', n_threads='0', random_seed=None):
    if random_seed != None:
        learnmodel_cmdline = "java -Xmx20g -jar ChromHMM/ChromHMM.jar LearnModel -printposterior -r 400 -init random -s {} -p {} {} {} {} {}".format(
            random_seed, n_threads, binary_input_dir, output_dir, num_labels, assembly)
    else:
        learnmodel_cmdline = "java -Xmx20g -jar ChromHMM/ChromHMM.jar LearnModel -printposterior -p {} {} {} {} {}".format(
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

def ChrHMM_read_posteriordir(posteriordir, resolution=200):
    '''
    for each file in posteriordir
    Initialize emptybins based on chromsizes
    fill in the posterior values for each slot
    return DF
    '''
    ls = os.listdir(posteriordir)
    to_parse = []
    for f in ls:
        # if rep in f:
        to_parse.append(f)  

    parsed_posteriors = {}
    for f in to_parse:
        fileinfo = f.split("_")
        if "chr" in fileinfo[-2]:
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
    # essential_tracks = ['H3K9me3', 'H3K27me3', 'H3K36me3', 'H3K4me3', 'H3K4me1', 'H3K27ac']
    # essential_tracks = ["H3K4me3"]
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

    assaylist = [tr for tr in os.listdir(CellType_dir) if os.path.isdir(CellType_dir+'/'+tr)]# and tr in essential_tracks]
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
            [tfmd.loc['assay', 'rep1_alig'], "rep1", str(tfmd.loc['accession', 'rep1_alig'])+".bed"])  
        navigate.append(
            [tfmd.loc['assay', 'rep2_alig'], "rep2", str(tfmd.loc['accession', 'rep2_alig'])+".bed"])

        navigate.append(
            [tfmd.loc['assay', 'rep1_alig'], "rep1_psd1", str(tfmd.loc['accession', 'rep1_alig'])+"_psdrep1.bed"])
        navigate.append(
            [tfmd.loc['assay', 'rep1_alig'], "rep1_psd2", str(tfmd.loc['accession', 'rep1_alig'])+"_psdrep2.bed"])
        navigate.append(
            [tfmd.loc['assay', 'rep2_alig'], "rep2_psd1", str(tfmd.loc['accession', 'rep2_alig'])+"_psdrep1.bed"])
        navigate.append(
            [tfmd.loc['assay', 'rep2_alig'], "rep2_psd2", str(tfmd.loc['accession', 'rep2_alig'])+"_psdrep2.bed"])

    
    cmft_concat = open("chmmfiles/{}/cmft_concat.txt".format(celltype_name), "w")
    cmft_rep1 = open("chmmfiles/{}/cmft_rep1.txt".format(celltype_name), "w")
    cmft_rep2 = open("chmmfiles/{}/cmft_rep2.txt".format(celltype_name), "w")

    cmft_rep1psd1 = open("chmmfiles/{}/cmft_rep1psd1.txt".format(celltype_name), "w")
    cmft_rep1psd2 = open("chmmfiles/{}/cmft_rep1psd2.txt".format(celltype_name), "w")
    cmft_rep2psd1 = open("chmmfiles/{}/cmft_rep2psd1.txt".format(celltype_name), "w")
    cmft_rep2psd2 = open("chmmfiles/{}/cmft_rep2psd2.txt".format(celltype_name), "w")

    cmft_rep1psd_concat = open("chmmfiles/{}/cmft_rep1psd_concat.txt".format(celltype_name), "w")
    cmft_rep2psd_concat = open("chmmfiles/{}/cmft_rep2psd_concat.txt".format(celltype_name), "w")
    
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

        elif ins[1] == "rep1_psd1":
            cmft_rep1psd1.write("{}_{}\t{}\t{}\n".format(
                celltype_name, ins[1], ins[0], ins[2]))
            cmft_rep1psd_concat.write("{}_{}\t{}\t{}\n".format(
                celltype_name, ins[1], ins[0], ins[2]))
        elif ins[1] == "rep1_psd2":
            cmft_rep1psd2.write("{}_{}\t{}\t{}\n".format(
                celltype_name, ins[1], ins[0], ins[2]))
            cmft_rep1psd_concat.write("{}_{}\t{}\t{}\n".format(
                celltype_name, ins[1], ins[0], ins[2]))
        elif ins[1] == "rep2_psd1":
            cmft_rep2psd1.write("{}_{}\t{}\t{}\n".format(
                celltype_name, ins[1], ins[0], ins[2]))
            cmft_rep2psd_concat.write("{}_{}\t{}\t{}\n".format(
                celltype_name, ins[1], ins[0], ins[2]))
        elif ins[1] == "rep2_psd2":
            cmft_rep2psd2.write("{}_{}\t{}\t{}\n".format(
                celltype_name, ins[1], ins[0], ins[2]))
            cmft_rep2psd_concat.write("{}_{}\t{}\t{}\n".format(
                celltype_name, ins[1], ins[0], ins[2]))
    
    cmft_concat.close()
    cmft_rep1.close()
    cmft_rep2.close()

    cmft_rep1psd1.close()
    cmft_rep1psd2.close()
    cmft_rep2psd1.close()
    cmft_rep2psd2.close()

    cmft_rep1psd_concat.close()
    cmft_rep2psd_concat.close()

def ChromHMM_replicate_runs(chmm_celltype_dir, chmm_output_dir, n_thread='0', num_labels=16):
    namesig = chmm_celltype_dir.split("/")[-1]

    if os.path.exists(chmm_celltype_dir+"/binarized_rep1") == False:
        binarize_data(
            chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep1.txt", chmm_celltype_dir+"/binarized_rep1", 
            resolution=200)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1") == False:
        learnModel(
            chmm_celltype_dir+"/binarized_rep1", chmm_output_dir+"/"+namesig+"_rep1", 
            num_labels=num_labels, assembly='hg38', n_threads=n_thread, random_seed=None)

    try:
        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1/parsed_posterior.csv") == False:
            parsed_posterior = ChrHMM_read_posteriordir(
                chmm_output_dir+"/"+namesig+"_rep1/POSTERIOR", "rep1", resolution=200)
            
            parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep1/parsed_posterior.csv")
    except:
        pass

    ####=================================================================================================####
    if os.path.exists(chmm_celltype_dir+"/binarized_rep2") == False:
        binarize_data(
            chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep2.txt", chmm_celltype_dir+"/binarized_rep2", 
            resolution=200)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2") == False:
        learnModel(
            chmm_celltype_dir+"/binarized_rep2", chmm_output_dir+"/"+namesig+"_rep2", 
            num_labels=num_labels, assembly='hg38', n_threads=n_thread, random_seed=None)

    try:
        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2/parsed_posterior.csv") == False:
            parsed_posterior = ChrHMM_read_posteriordir(
                chmm_output_dir+"/"+namesig+"_rep2/POSTERIOR", "rep2", resolution=200)
            
            parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep2/parsed_posterior.csv")
    except:
        pass

def ChromHMM_pseudoreplicate_runs(chmm_celltype_dir, chmm_output_dir, n_thread='0', num_labels=16):
    namesig = chmm_celltype_dir.split("/")[-1]

    if os.path.exists(chmm_celltype_dir+"/binarized_rep1psd_concat") == False:
        binarize_data(
            chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep1psd_concat.txt", chmm_celltype_dir+"/binarized_rep1psd_concat", 
            resolution=200)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1psd_concat") == False:
        learnModel(
            chmm_celltype_dir+"/binarized_rep1psd_concat", chmm_output_dir+"/"+namesig+"_rep1psd_concat", 
            num_labels=num_labels, assembly='hg38', n_threads=n_thread, random_seed=None)

    try:
        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1psd_concat/parsed_posterior_rep1_psd1.csv") == False:
            parsed_posterior = ChrHMM_read_posteriordir(
                chmm_output_dir+"/"+namesig+"_rep1psd_concat/POSTERIOR", "rep1_psd1", resolution=200)
            
            parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep1psd_concat/parsed_posterior_rep1_psd1.csv")
    except:
        pass

    try:
        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1psd_concat/parsed_posterior_rep1_psd2.csv") == False:
            parsed_posterior = ChrHMM_read_posteriordir(
                chmm_output_dir+"/"+namesig+"_rep1psd_concat/POSTERIOR", "rep1_psd2", resolution=200)
            
            parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep1psd_concat/parsed_posterior_rep1_psd2.csv")
    except:
        pass

    ####=================================================================================================####
    # namesig = chmm_celltype_dir.split("/")[-1]

    # if os.path.exists(chmm_celltype_dir+"/binarized_rep2psd_concat") == False:
    #     binarize_data(
    #         chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep2psd_concat.txt", chmm_celltype_dir+"/binarized_rep2psd_concat", 
    #         resolution=200)

    # if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2psd_concat") == False:
    #     learnModel(
    #         chmm_celltype_dir+"/binarized_rep2psd_concat", chmm_output_dir+"/"+namesig+"_rep2psd_concat", 
    #         num_labels=num_labels, assembly='hg38', n_threads=n_thread, random_seed=None)

    # try:
    #     if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2psd_concat/parsed_posterior_rep2_psd1.csv") == False:
    #         parsed_posterior = ChrHMM_read_posteriordir(
    #             chmm_output_dir+"/"+namesig+"_rep2psd_concat/POSTERIOR", "rep2_psd1", resolution=200)
            
    #         parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep2psd_concat/parsed_posterior_rep2_psd1.csv")
    # except:
    #     pass

    # try:
    #     if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2psd_concat/parsed_posterior_rep2_psd2.csv") == False:
    #         parsed_posterior = ChrHMM_read_posteriordir(
    #             chmm_output_dir+"/"+namesig+"_rep2psd_concat/POSTERIOR", "rep2_psd2", resolution=200)
            
    #         parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep2psd_concat/parsed_posterior_rep2_psd2.csv")
    # except:
    #     pass

def ChromHMM_paraminit_multi_subruns(chmm_celltype_dir, chmm_output_dir, random_seeds,  num_subruns=10, n_thread='0', num_labels=16):
    """
    initialize s_n=10 subrun_seeds with parent-seed {rs}
    run all 10 instances

    look at all of their -log-likelihood
    pick the best
    remove the others
    name it namesig+_rep*_{rs}
    parse posterior for namesig+_rep*_{rs} only
    """

    namesig = chmm_celltype_dir.split("/")[-1]
    if "/" in namesig:
        namesig.replace("/", "")

    if os.path.exists(chmm_celltype_dir+"/binarized_rep1") == False:
        binarize_data(
            chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep1.txt", 
            chmm_celltype_dir+"/binarized_rep1", resolution=200)

    for rs in random_seeds: 
        np.random.seed(rs)
        sub_randseeds = np.random.randint(1000, size=int(num_subruns))
        list_of_subrun_dirs = []
        for srs in sub_randseeds:
            learnModel(
                chmm_celltype_dir+"/binarized_rep1", chmm_output_dir+"/"+namesig+"_rep1_rs{}_srs{}".format(rs, srs), 
                num_labels=num_labels, assembly='hg38', n_threads=n_thread, random_seed=int(srs))

            list_of_subrun_dirs.append(chmm_output_dir+"/"+namesig+"_rep1_rs{}_srs{}".format(rs, srs))
    
        # now browse subruns
        nlls = {}
        for subrun in list_of_subrun_dirs:
            lfiles = os.listdir(subrun)
            for ll in lfiles:
                if "model" in ll:
                    nlls[subrun] = float(open(subrun+"/"+ll, 'r').readlines()[0].split("\t")[3])

        # pick the best run
        best_run = max(nlls, key=nlls.get)
        os.system("mv {} {}".format(
            best_run,
            chmm_output_dir+"/"+namesig+"_rep1_rs{}".format(rs)))

        for k in nlls.keys():
            if os.path.isdir(k):
                os.system("rm -r {}".format(k))
        
        # now parse posteriors
        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1_rs{}/parsed_posterior.csv".format(rs)) == False:
                parsed_posterior = ChrHMM_read_posteriordir(
                    chmm_output_dir+"/"+namesig+"_rep1_rs{}/POSTERIOR".format(rs), "rep1", resolution=200)
                parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep1_rs{}/parsed_posterior.csv".format(rs))
    
    ####=================================================================================================####
    # if os.path.exists(chmm_celltype_dir+"/binarized_rep2") == False:
    #     binarize_data(
    #         chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep2.txt", 
    #         chmm_celltype_dir+"/binarized_rep2", resolution=200)
            
    # for rs in random_seeds: 
        # list_of_subrun_dirs = []
        # np.random.seed(rs)
        # sub_randseeds = np.random.randint(1000, size=int(num_subruns))
    #     for srs in sub_randseeds:
    #         learnModel(
    #             chmm_celltype_dir+"/binarized_rep2", chmm_output_dir+"/"+namesig+"_rep2_rs{}_srs{}".format(rs, srs), 
    #             num_labels=num_labels, assembly='hg38', n_threads=n_thread, random_seed=int(srs))

    #         list_of_subrun_dirs.append(chmm_output_dir+"/"+namesig+"_rep2_rs{}_srs{}".format(rs, srs))
    
        # now browse subruns
        # nlls = {}
        # for subrun in list_of_subrun_dirs:
        #     lfiles = os.listdir(subrun)
        #     for ll in lfiles:
        #         if "model" in ll:
        #             nlls[subrun] = float(open(subrun+"/"+ll, 'r').readlines()[0].split("\t")[3])

        # # pick the best run
        # best_run = max(nlls, key=nlls.get)
        # os.system("mv {} {}".format(
        #     chmm_output_dir+"/"+namesig+"_rep2_rs{}_srs{}".format(rs, srs),
        #     chmm_output_dir+"/"+namesig+"_rep2_rs{}".format(rs)))

        # for k in nlls.keys():
        #     if os.path.isdir(k):
        #         os.system("rm -r {}".format(k))
        
        # # now parse posteriors
        # if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2_rs{}/parsed_posterior.csv".format(rs)) == False:
        #         parsed_posterior = ChrHMM_read_posteriordir(
        #             chmm_output_dir+"/"+namesig+"_rep2_rs{}/POSTERIOR".format(rs), "rep2", resolution=200)
        #         parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep2_rs{}/parsed_posterior.csv".format(rs))

def ChromHMM_paraminit_runs(chmm_celltype_dir, chmm_output_dir, random_seeds,  n_thread='0', num_labels=16):
    namesig = chmm_celltype_dir.split("/")[-1]
    if "/" in namesig:
        namesig.replace("/", "")

    for rs in random_seeds:
        if os.path.exists(chmm_celltype_dir+"/binarized_rep1") == False:
            binarize_data(
                chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep1.txt", chmm_celltype_dir+"/binarized_rep1", 
                resolution=200)
        
        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1_rs{}".format(rs)) == False:
            learnModel(
                chmm_celltype_dir+"/binarized_rep1", chmm_output_dir+"/"+namesig+"_rep1_rs{}".format(rs), 
                num_labels=num_labels, assembly='hg38', n_threads=n_thread, random_seed=int(rs))

        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep1_rs{}/parsed_posterior.csv".format(rs)) == False:
            parsed_posterior = ChrHMM_read_posteriordir(
                chmm_output_dir+"/"+namesig+"_rep1_rs{}/POSTERIOR".format(rs), "rep1", resolution=200)
            parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep1_rs{}/parsed_posterior.csv".format(rs))
        
        ####=================================================================================================####

        if os.path.exists(chmm_celltype_dir+"/binarized_rep2") == False:
            binarize_data(
                chmm_celltype_dir, chmm_celltype_dir+"/cmft_rep2.txt", chmm_celltype_dir+"/binarized_rep2", 
                resolution=200)
        
        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2_rs{}".format(rs)) == False:
            learnModel(
                chmm_celltype_dir+"/binarized_rep2", chmm_output_dir+"/"+namesig+"_rep2_rs{}".format(rs), 
                num_labels=num_labels, assembly='hg38', n_threads=n_thread, random_seed=int(rs))

        if os.path.exists(chmm_output_dir+"/"+namesig+"_rep2_rs{}/parsed_posterior.csv".format(rs)) == False:
            parsed_posterior = ChrHMM_read_posteriordir(
                chmm_output_dir+"/"+namesig+"_rep2_rs{}/POSTERIOR".format(rs), "rep2", resolution=200)
            parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_rep2_rs{}/parsed_posterior.csv".format(rs))

def ChromHMM_concat_runs(chmm_celltype_dir, chmm_output_dir, n_thread='0', num_labels=16):
    namesig = chmm_celltype_dir.split("/")[-1]
    if "/" in namesig:
        namesig.replace("/", "")
    
    if os.path.exists(chmm_celltype_dir+"/binarized_concat") == False:
        binarize_data(
            chmm_celltype_dir, chmm_celltype_dir+"/cmft_concat.txt", chmm_celltype_dir+"/binarized_concat", 
            resolution=200)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_concat") == False:
        learnModel(
            chmm_celltype_dir+"/binarized_concat", chmm_output_dir+"/"+namesig+"_concat", 
            num_labels=num_labels, assembly='hg38', n_threads=n_thread, random_seed=None)

    if os.path.exists(chmm_output_dir+"/"+namesig+"_concat/parsed_posterior_rep1.csv") == False:
        parsed_posterior = ChrHMM_read_posteriordir(
            chmm_output_dir+"/"+namesig+"_concat/POSTERIOR", "rep1", resolution=200)
        
        parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_concat/parsed_posterior_rep1.csv")
    
    if os.path.exists(chmm_output_dir+"/"+namesig+"_concat/parsed_posterior_rep2.csv") == False:
        parsed_posterior = ChrHMM_read_posteriordir(
            chmm_output_dir+"/"+namesig+"_concat/POSTERIOR", "rep2", resolution=200)

        parsed_posterior.to_csv(chmm_output_dir+"/"+namesig+"_concat/parsed_posterior_rep2.csv")
