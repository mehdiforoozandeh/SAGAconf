from signal import SIG_DFL
from _utils import *
from _pipeline import *
from _reproducibility import *
from _interpret import *
from _visualize import *
from _concat import *
from _chromhmm import *
from get_data import *
from matchlabels import *
import pandas as pd
import glob, os 
import multiprocessing as mp
from functools import partial
import matplotlib
'''

In this file, the step-by-step execution of analysis-funcitons will be managed.
Basically, by importing and managing functions from different modules
and classes, we will execute every step of the methods section. 

1 - get data from encode : USE 
    search_encode(cell, download_dir, target_assembly="GRCh38")

SEGWAY:
    2 - create genomedata (separate and concatenated) files and trackname_assay.txt
    3 - run segway (separate and concatenated)

ChromHMM:
    2 - Binarize data
    3 - run chromhmm
    
4 - parse and bin results
4.5  - biological label interpretation
5 - match labels (for separate runs)
6 - cluster labels
7 - reproducibility analysis (at different levels of abstraction)

'''

def download_encode_files(celltype_list, download_dir, target_assembly):
    for celltype in celltype_list:
        search_encode(celltype, download_dir, target_assembly=target_assembly)

def check_if_data_exists(celltype_list, download_dir):
    exists_bool = []
    for ct in celltype_list:
        if ' ' in ct:
            ct = ct.replace(' ', '_')

        if os.path.isdir(download_dir+'/'+ct):
            exists_bool.append(True)

        else:
            exists_bool.append(False)

    return exists_bool

def read_list_of_assays(celltype_dir):
    ls = os.listdir(celltype_dir)
    list_of_assays = []
    for d in ls:
        if os.path.isdir(celltype_dir + '/' + d):
            list_of_assays.append(d)
    
    return list_of_assays

def Convert_all_BW2BG(celltype_dir):
    '''
    Creats BedGraph files from all BigWig files
    found within a celltype directory
    '''
    if os.path.exists("bigWigToBedGraph") ==False:
        url = 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bigWigToBedGraph'
        file_dl_response = requests.get(url, allow_redirects=True)
        open(download_dir+"bigWigToBedGraph", 'wb').write(file_dl_response.content)

    ls = os.listdir(celltype_dir)
    if " " in celltype_dir:
        celltype_dir = celltype_dir.replace(" ","\ ")
    for f in ls:
        if ".bigWig" in f:
            if f.replace(".bigWig",".bedGraph") not in ls:
                os.system(
                    "./bigWigToBedGraph {}.bigWig {}.bedGraph".format(
                        celltype_dir+'/'+f.replace(".bigWig",""), 
                        celltype_dir+'/'+f.replace(".bigWig","")))

def read_metadata(donwload_dir):
    metadata = {}
    for ct in CellType_list:
        tracks_list = [tr for tr in os.listdir(download_dir+ct) if os.path.isdir(download_dir+ct+'/'+tr)]
        for tr in tracks_list:
            tfmd = pd.read_csv(download_dir+ct+'/'+tr+'/track_files_metadata.csv')
            tfmd.index = list(tfmd['Unnamed: 0'])
            tfmd = tfmd.drop('Unnamed: 0', axis=1)
            if ct in list(metadata.keys()):
                metadata[ct].append(tfmd)
            else:
                metadata[ct] = [tfmd]

    return metadata

def gather_replicates(celltype_dir):
    '''
    for each celltype
    create two files (one for each replicate)
    each file contains the following

    segway_rep1_{biosampleID}.txt:
        {ASSAY 1}     {File 1}
        {ASSAY 2}     {File 2}
        ...

    segway_rep2_{biosampleID}.txt:
        {ASSAY 1}     {File 1}
        {ASSAY 2}     {File 2}
        ...
    
    SAME TWO FILES FOR CHROMHMM...
    '''

    tracks_list = [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)]
    with open(celltype_dir+'/replicate_number.txt', 'r') as repguide_file:
        lines = repguide_file.readlines()
        rep1_biosampleID = lines[0][5:-1]
        rep2_biosampleID = lines[1][5:]

    # print("rep1_biosampleID:\t{}\nrep2_biosampleID:\t{}\n".format(rep1_biosampleID, rep2_biosampleID))
    nav_rep1_seg = {}
    nav_rep2_seg = {}

    nav_rep1_chmm = {}
    nav_rep2_chmm = {}

    for tr in tracks_list:
        tfmd = pd.read_csv(celltype_dir+'/'+tr+'/track_files_metadata.csv')
        tfmd.index = list(tfmd['Unnamed: 0'])
        tfmd = tfmd.drop('Unnamed: 0', axis=1)

        assert str(tfmd.loc['biosample', 'rep1_fcoc']) == rep1_biosampleID
        assert str(tfmd.loc['biosample', 'rep2_fcoc']) == rep2_biosampleID

        nav_rep1_seg[tfmd.loc['assay', 'rep1_fcoc']] = '{}/{}/{}.bedGraph'.format(
            celltype_dir, tr, tfmd.loc['accession', 'rep1_fcoc'])

        nav_rep2_seg[tfmd.loc['assay', 'rep2_fcoc']] = '{}/{}/{}.bedGraph'.format(
            celltype_dir, tr, tfmd.loc['accession', 'rep2_fcoc'])

        nav_rep1_chmm[tfmd.loc['assay', 'rep1_alig']] = '{}/{}/{}.bam'.format(
            celltype_dir, tr, tfmd.loc['accession', 'rep1_alig'])

        nav_rep2_chmm[tfmd.loc['assay', 'rep2_alig']] = '{}/{}/{}.bam'.format(
            celltype_dir, tr, tfmd.loc['accession', 'rep2_alig'])

    with open(celltype_dir+"/seg_rep1_{}.txt".format(rep1_biosampleID), "w") as segrep1_file:
        for k, v in nav_rep1_seg.items():
            segrep1_file.write("{}\t{}\n".format(k, v))

    with open(celltype_dir+"/seg_rep2_{}.txt".format(rep2_biosampleID), "w") as segrep2_file:
        for k, v in nav_rep2_seg.items():
            segrep2_file.write("{}\t{}\n".format(k, v))

    with open(celltype_dir+"/chmm_rep1_{}.txt".format(rep1_biosampleID), "w") as chmmrep1_file:
        for k, v in nav_rep1_chmm.items():
            chmmrep1_file.write("{}\t{}\n".format(k, v))

    with open(celltype_dir+"/chmm_rep2_{}.txt".format(rep2_biosampleID), "w") as chmmrep2_file:
        for k, v in nav_rep2_chmm.items():
            chmmrep2_file.write("{}\t{}\n".format(k, v))

def bamtobed(download_dir):
    if os.path.exists("bedtools")==False:
        bt_url = "https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary"
        bt_file_dl_response = requests.get(bt_url, allow_redirects=True)
        open("bedtools.static.binary", 'wb').write(bt_file_dl_response.content)
        print('downloaded the bedtools')
        os.system("mv bedtools.static.binary bedtools")
        os.system("chmod a+x bedtools")
    CellType_list = [download_dir+ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+"/"+ct)]
    for ct in CellType_list:
        assay_list = [assay for assay in os.listdir(ct) if os.path.isdir(ct+"/"+assay)]
        for ass in assay_list:
            bam_list = [bm for bm in os.listdir(ct+"/"+ass) if ".bam" in ct+"/"+ass+'/'+bm]
            for bm in bam_list:
                if os.path.exists(ct+"/"+ass+'/'+bm.replace(".bam", '.bed')) == False:
                    os.system(
                        "bedtools bamtobed -i {} > {}".format(
                            ct+"/"+ass+'/'+bm, ct+"/"+ass+'/'+bm.replace(".bam", '.bed')
                        )
                    )
                    print("converted ", ct+"/"+ass+'/'+bm, "->", ct+"/"+ass+'/'+bm.replace(".bam", '.bed'))

def create_genomedata(celltype_dir, sequence_file):
    '''
    read segway replicate files.
    gather bedGraph file paths

    create genomedata file_rep1 in celltype_dir
    create genomedata file_rep2 in celltype_dir
    '''
    print('creating genomedata files for {}...'.format(celltype_dir))
    tracks_bgs_rep1 = {}
    with open(glob.glob(celltype_dir+'/seg_rep1_*')[0], 'r') as rep1_guidefile:
        lines_1 = rep1_guidefile.readlines()
        for l in lines_1:
            tracks_bgs_rep1[l.split('\t')[0]] = l.split('\t')[1][:-1]
    
    tracklist_rep1 = ''
    for k, v in tracks_bgs_rep1.items():
        tracklist_rep1 = tracklist_rep1 + '-t {}={} '.format(k, v)
    
    os.system(
        'genomedata-load -s {} --sizes {} --verbose {}.genomedata'.format(
            sequence_file, tracklist_rep1, celltype_dir+'/rep1'))

    tracks_bgs_rep2 = {}
    with open(glob.glob(celltype_dir+'/seg_rep2_*')[0], 'r') as rep2_guidefile:
        lines_2 = rep2_guidefile.readlines()
        for l in lines_2:
            tracks_bgs_rep2[l.split('\t')[0]] = l.split('\t')[1][:-1]
    
    tracklist_rep2 = ''
    for k, v in tracks_bgs_rep2.items():
        tracklist_rep2 = tracklist_rep2 + '-t {}={} '.format(k, v)
    
    os.system(
        'genomedata-load -s {} --sizes {} --verbose {}.genomedata'.format(
            sequence_file, tracklist_rep2, celltype_dir+'/rep2'))

def concat_create_genomedata(celltype_dir, sequence_file):
    tracks_bgs_rep1 = {}
    with open(glob.glob(celltype_dir+'/seg_rep1_*')[0], 'r') as rep1_guidefile:
        lines_1 = rep1_guidefile.readlines()
        for l in lines_1:
            tracks_bgs_rep1[l.split('\t')[0]] = l.split('\t')[1][:-1]

    tracks_bgs_rep2 = {}
    with open(glob.glob(celltype_dir+'/seg_rep2_*')[0], 'r') as rep2_guidefile:
        lines_2 = rep2_guidefile.readlines()
        for l in lines_2:
            tracks_bgs_rep2[l.split('\t')[0]] = l.split('\t')[1][:-1]
    
    # rename chromosomes for bgs in dictionaries
    for k,v in tracks_bgs_rep1.items():
        if os.path.exists(v.replace("bedGraph","_concat.bedGraph")) == False:
            concat_rename_chroms(v, v.replace("bedGraph","_concat.bedGraph"), 1)

    for k,v in tracks_bgs_rep2.items():
        if os.path.exists(v.replace("bedGraph","_concat.bedGraph")) == False:
            concat_rename_chroms(v, v.replace("bedGraph","_concat.bedGraph"), 2)

    virtualized_tracks = {}
    # virtualize bedgraphs into concatenated tracks    
    for k in tracks_bgs_rep1.keys():
        if os.path.exists(celltype_dir+ "/" + k + "/{}_concatenated.bedGraph".format(k)) == False:
            concat_virtualize_chroms(
                tracks_bgs_rep1[k].replace("bedGraph","_concat.bedGraph"),
                tracks_bgs_rep2[k].replace("bedGraph","_concat.bedGraph"),
                celltype_dir+ "/" + k + "/{}_concatenated.bedGraph".format(k))

        virtualized_tracks[k] = celltype_dir+ "/" + k + "/{}_concatenated.bedGraph".format(k)

    if os.path.exists(sequence_file.replace(".sizes", "_concatenated.sizes")) == False:
        concat_virtualize_include_file(
            sequence_file, sequence_file.replace(".sizes", "_concatenated.sizes"), 
            write_separate=True)

    if os.path.exists(celltype_dir + "/" + "concatenated.genomedata") == False:
        tracklist = ''
        for k, v in virtualized_tracks.items():
            tracklist = tracklist + '-t {}={} '.format(k, v)  

        os.system('genomedata-load -s {} --sizes {} --verbose {}'.format(
            sequence_file.replace(".sizes", "_concatenated.sizes"), 
            tracklist, celltype_dir + "/" + "concatenated.genomedata"))

    if os.path.exists(celltype_dir + "/" + "concat_rep1.genomedata") == False:
        tracklist = ''
        for k, v in tracks_bgs_rep1.items():
            tracklist = tracklist + '-t {}={} '.format(k, v.replace("bedGraph","_concat.bedGraph"))  

        os.system('genomedata-load -s {} --sizes {} --verbose {}'.format(
            sequence_file.replace(".sizes", "_VirtRep1_concatenated.sizes"), 
            tracklist, celltype_dir + "/" + "concat_rep1.genomedata"))

    if os.path.exists(celltype_dir + "/" + "concat_rep2.genomedata") == False:
        tracklist = ''
        for k, v in tracks_bgs_rep2.items():
            tracklist = tracklist + '-t {}={} '.format(k, v.replace("bedGraph","_concat.bedGraph"))  

        os.system('genomedata-load -s {} --sizes {} --verbose {}'.format(
            sequence_file.replace(".sizes", "_VirtRep2_concatenated.sizes"), 
            tracklist, celltype_dir + "/" + "concat_rep2.genomedata"))


def concat_RunParse_segway(celltype_dir, output_dir, random_seed=73):
    '''
    define params_dict
    run concat seg + post clean up
    '''
    name_sig = celltype_dir.split('/')[-1]
    num_tracks = len(
        [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)])

    concat_param_dict = {
        "random_seed":random_seed, "track_weight":0.01,"stws":1, "ruler_scale":100, 
        "prior_strength":1, "resolution":100, "mini_batch_fraction":0.03,
        "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
        "name_sig_concat":output_dir+name_sig+'_concat', 
        "name_sig_rep1":output_dir+name_sig+'_concat_rep1', 
        "name_sig_rep2":output_dir+name_sig+'_concat_rep2',
        "genomedata_file_concat":celltype_dir+'/concatenated.genomedata',  
        "genomedata_file_rep1":celltype_dir+'/concat_rep1.genomedata', 
        "genomedata_file_rep2":celltype_dir+'/concat_rep2.genomedata',
        "traindir":output_dir+name_sig+'_concat_train', 
        "posteriordir_rep1":output_dir+name_sig+'_concat_posterior_1',
        "posteriordir_rep2":output_dir+name_sig+'_concat_posterior_2'}

    if os.path.exists(concat_param_dict['name_sig_rep1']) == False or \
         os.path.exists(concat_param_dict['name_sig_rep2']) == False:
        print('Running Segway concatenated celltype {}...'.format(celltype_dir))
        concat_segwayrun_and_postprocess(concat_param_dict)

def RunParse_segway_replicates(celltype_dir, output_dir, random_seed=73):
    name_sig = celltype_dir.split('/')[-1]
    num_tracks = len(
        [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)])

    params_dict_1 = {
        "random_seed":random_seed, "track_weight":0.01,
        "stws":1, "ruler_scale":100, "prior_strength":1, "resolution":100, 
        "mini_batch_fraction":0.03, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
        "name_sig":output_dir+name_sig+'_rep1', 
        "genomedata_file":celltype_dir+'/rep1.genomedata', 
        "traindir":output_dir+name_sig+'rep1'+'_train', 
        "posteriordir":output_dir+name_sig+'rep1'+'_posterior'
    }
    
    if os.path.exists(params_dict_1['name_sig']) == False:
        print('Running Segway celltype {} Rep1'.format(celltype_dir))
        run_segway_and_post_process(params_dict_1)

    else:
        print(params_dict_1['name_sig'], "already exists")

    params_dict_2 = {
        "random_seed":random_seed, "track_weight":0.01,
        "stws":1, "ruler_scale":100, "prior_strength":1, "resolution":100, 
        "mini_batch_fraction":0.03, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
        "name_sig":output_dir+name_sig+'_rep2', 
        "genomedata_file":celltype_dir+'/rep2.genomedata', 
        "traindir":output_dir+name_sig+'rep2'+'_train', 
        "posteriordir":output_dir+name_sig+'rep2'+'_posterior'
    }
    
    if os.path.exists(params_dict_2['name_sig']) == False:
        print('Running Segway celltype {} Rep2'.format(celltype_dir))
        run_segway_and_post_process(params_dict_2)

    else:
        print(params_dict_2['name_sig'], "already exists")

def RunParse_segway_param_init(celltype_dir, replicate_number, random_seeds, output_dir):
    # replicate number should be in format "repN" -> i.e. rep1, rep2
    name_sig = celltype_dir.split('/')[-1]
    num_tracks = len(
        [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)])

    params_dict_1 = {
        "random_seed":random_seeds[0], "track_weight":0.01,
        "stws":1, "ruler_scale":100, "prior_strength":1, "resolution":100, 
        "mini_batch_fraction":0.03, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
        "name_sig":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0]), 
        "genomedata_file":celltype_dir+'/{}.genomedata'.format(replicate_number), 
        "traindir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0])+'_train', 
        "posteriordir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0])+'_posterior'
    }

    params_dict_2 = {
        "random_seed":random_seeds[1], "track_weight":0.01,
        "stws":1, "ruler_scale":100, "prior_strength":1, "resolution":100, 
        "mini_batch_fraction":0.03, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
        "name_sig":output_dir+name_sig+'_{}_rs{}'.format(replicate_number,random_seeds[1]),
        "genomedata_file":celltype_dir+'/{}.genomedata'.format(replicate_number), 
        "traindir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0])+'_train', 
        "posteriordir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0])+'_posterior'
    }

    if os.path.exists(params_dict_1['name_sig']) == False:
        print('Running Segway parameter initialization test on celltype {}, {}, with random seed {}'.format(
            celltype_dir, replicate_number, random_seeds[0]))
        run_segway_and_post_process(params_dict_1)

    else:
        print(params_dict_1['name_sig'], "already exists")

    if os.path.exists(params_dict_2['name_sig']) == False:
        print('Running Segway parameter initialization test on celltype {}, {}, with random seed {}'.format(
            celltype_dir, replicate_number, random_seeds[1]))
        run_segway_and_post_process(params_dict_2)

    else:
        print(params_dict_2['name_sig'], "already exists")

def intersect_parsed_posteriors(parsed_df_dir_1, parsed_df_dir_2):
    intersect = pd.merge(
        pd.read_csv(parsed_df_dir_1).drop("Unnamed: 0", axis=1), 
        pd.read_csv(parsed_df_dir_2).drop("Unnamed: 0", axis=1), 
        how='inner', on=['chr', 'start', 'end'])

    df1 = [intersect.chr, intersect.start, intersect.end]
    df2 = [intersect.chr, intersect.start, intersect.end]

    for c in intersect.columns:
        if c[-1] == 'x':
            df1.append(intersect[c])
        elif c[-1] == 'y':
            df2.append(intersect[c])

    df1 = pd.concat(df1, axis=1)
    df1.columns = [c.replace("_x", "") for c in df1.columns]
    
    df2 = pd.concat(df2, axis=1)
    df2.columns = [c.replace("_y", "") for c in df2.columns]

    return df1, df2

def get_biological_mnemonics(results_directory, segway=True):
    '''
    - run segtools signal-distribution in the segway's output directory and 
        move the results to label_interpretation/segwayOutput/*.
    - run segtools feature-aggregation in label_interpretation/segwayOutput/*
    - move the corresponding trackname_assays.txt
    - run apply_samples.py
    - get mnemonics and move back to segway's output directory'''
    pass

def read_mnemonics(mnemon_file):
    df = pd.read_csv(mnemon_file, sep="\t")
    mnemon = []
    for i in range(len(df)):
        mnemon.append(str(df["old"][i])+"_"+df["new"][i])

    return mnemon

def _report_reproducibility(loci_1, loci_2, pltsavedir, general=True, num_bins=20, merge_cc_curves=False, full_report=True):
    '''
    TODO:
    1- instead of showing the plot, savefig (vector)
    2- save the visualized data for each plot in text,csv format to enable 
    easy reproduction of better visualize with other visualization tools
    '''
    if os.path.exists(pltsavedir)==False:
        os.mkdir(pltsavedir)

    if os.path.exists(pltsavedir+"/agreement")==False:
        os.mkdir(pltsavedir+"/agreement")

    if os.path.exists(pltsavedir+"/sankey")==False:
        os.mkdir(pltsavedir+"/sankey")

    agr = Agreement(loci_1, loci_2, pltsavedir+"/agreement")
    vis = sankey(loci_1, loci_2, pltsavedir+"/sankey")

    to_report = [
        agr.general_agreement(), agr.general_OE_ratio(log_transform=False), agr.general_cohens_kappa()]

    print('general agreement:    ', to_report[0])
    print('general o/e ratio:    ', to_report[1])
    print('general CK:    ', to_report[2])

    agr.plot_agreement()
    agr.plot_CK()
    agr.plot_OE()
    vis.sankey_diag()

    if full_report:
        if os.path.exists(pltsavedir+"/cc")==False:
            os.mkdir(pltsavedir+"/cc")

        if os.path.exists(pltsavedir+"/calib")==False:
            os.mkdir(pltsavedir+"/calib")

        cc = correspondence_curve(loci_1, loci_2, pltsavedir+"/cc")
        # repr = Reprodroducibility_vs_posterior(loci_1,loci_2, pltsavedir+"/calib",log_transform=False)

        cc.plot_curve(plot_general=general, merge_plots=merge_cc_curves)
        repr.per_label_count_independent(num_bins=num_bins)

        if general:
            repr.general_count_independent(num_bins=num_bins)
        
        return to_report


def _post_clustering(loci_1, loci_2, pltsavedir, OE_transform=True):
    num_labels = loci_1.shape[1]-3

    print('generating confmat 1')
    
    conf_mat = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=False)

    assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')

    print('connecting barpartite')

    corrected_loci_1, corrected_loci_2 = \
        connect_bipartite(loci_1, loci_2, assignment_pairs)
    
    print('connected barpartite')

    del loci_1, loci_2

    print('generating confmat 2')
    conf_mat = confusion_matrix(
        corrected_loci_1, corrected_loci_2, num_labels, 
        OE_transform=OE_transform, symmetric=True)  

    mat_max = conf_mat.max(axis=1).max(axis=0)
    mat_min = conf_mat.min(axis=1).min(axis=0)

    conf_mat = (conf_mat - mat_min) / (mat_max - mat_min)
    distance_matrix = 1 - conf_mat
    linkage = hc.linkage(distance_matrix, method='average')

    if os.path.exists(pltsavedir)==False:
        os.mkdir(pltsavedir)

    # initial reproducibility without cluster merging
    if os.path.exists(pltsavedir+"/{}_labels".format(num_labels)) == False:
            os.mkdir(pltsavedir+"/{}_labels".format(num_labels))

    stepwise_report = {}
    
    stepwise_report[num_labels] = report_reproducibility(
            corrected_loci_1, corrected_loci_2, pltsavedir+"/{}_labels".format(num_labels), 
            general=False, num_bins=20, merge_cc_curves=False)

    # sns.clustermap(
    #     distance_matrix, row_linkage=linkage, 
    #     col_linkage=linkage, annot=True)

    # if os.path.exists(pltsavedir+'/post_clustering/')==False:
    #     os.mkdir(pltsavedir+'/post_clustering/')

    # plt.savefig('{}/post_clustering/clustermap.pdf'.format(pltsavedir), format='pdf')
    # plt.savefig('{}/post_clustering/clustermap.svg'.format(pltsavedir), format='svg')
    # plt.clf()

    for m in range(len(linkage)):
        # merging clusters one at a time
        to_be_merged = [
            "posterior{}".format(int(linkage[m, 0])), 
            "posterior{}".format(int(linkage[m, 1]))]

        print("to be merged", to_be_merged)
        corrected_loci_1['posterior{}'.format(num_labels + m)] = \
            corrected_loci_1[to_be_merged[0]] + corrected_loci_1[to_be_merged[1]]

        corrected_loci_1 = corrected_loci_1.drop(to_be_merged, axis=1)

        corrected_loci_2['posterior{}'.format(num_labels + m)] = \
            corrected_loci_2[to_be_merged[0]] + corrected_loci_2[to_be_merged[1]]

        corrected_loci_2 = corrected_loci_2.drop(to_be_merged, axis=1)

        # print("merged loci_1", corrected_loci_1, "\nmerged loci_2", corrected_loci_2)

        if os.path.exists(pltsavedir+"/"+str(num_labels-(m+1))+'_labels') == False:
            os.mkdir(pltsavedir+"/"+str(num_labels-(m+1))+'_labels')

        stepwise_report[int(num_labels-(m+1))] = report_reproducibility(
            corrected_loci_1, corrected_loci_2, pltsavedir+"/"+str(num_labels-(m+1))+'_labels', 
            general=False, num_bins=20, merge_cc_curves=False, full_report=True)

    xaxis, yaxis = [], []
    for k, v in stepwise_report.items():
        xaxis.append(k)
        yaxis.append(v[0])

    plt.bar(xaxis, yaxis)
    plt.xlabel('Number of Labels')
    plt.ylabel('General Agreement')
    plt.savefig('{}/postclustering_agreement.svg'.format(pltsavedir), format='svg')
    plt.savefig('{}/postclustering_agreement.pdf'.format(pltsavedir), format='pdf')
    plt.clf()
    
    xaxis, yaxis = [], []
    for k, v in stepwise_report.items():
        xaxis.append(k)
        yaxis.append(v[1])

    plt.bar(xaxis, yaxis)
    plt.xlabel('Number of Labels')
    plt.ylabel('O/E Agreement')
    plt.savefig('{}/postclustering_oe_agreement.svg'.format(pltsavedir), format='svg')
    plt.savefig('{}/postclustering_oe_agreement.pdf'.format(pltsavedir), format='pdf')
    plt.clf()

    xaxis, yaxis = [], []
    for k, v in stepwise_report.items():
        xaxis.append(k)
        yaxis.append(v[2])

    plt.bar(xaxis, yaxis)
    plt.xlabel('Number of Labels')
    plt.ylabel("Cohen's Kappa")
    plt.savefig('{}/postclustering_ck.svg'.format(pltsavedir), format='svg')
    plt.savefig('{}/postclustering_ck.pdf'.format(pltsavedir), format='pdf')
    plt.clf()
    return stepwise_report

def report_reproducibility(loci_1, loci_2, pltsavedir):
    """
    get basic reproducibility results for a pair of 
    experiments (two locis)
    """
    if os.path.exists(pltsavedir) == False:
        os.mkdir(pltsavedir)

    if os.path.exists(pltsavedir+"/cc") == False:
        os.mkdir(pltsavedir+"/cc")

    if os.path.exists(pltsavedir+"/agr") == False:
        os.mkdir(pltsavedir+"/agr")
    
    if os.path.exists(pltsavedir+"/snk") == False:
        os.mkdir(pltsavedir+"/snk")

    if os.path.exists(pltsavedir+"/clb") == False:
        os.mkdir(pltsavedir+"/clb")
    
    to_report = {}

    # cc = correspondence_curve(loci_1, loci_2, pltsavedir+"/cc")
    # cc.plot_curve(plot_general=False, merge_plots=False)
    # del cc
    # plt.close("all")
    # plt.style.use('default')

    # calb = posterior_calibration(
    #     loci_1, loci_2, log_transform=False, ignore_overconf=False, filter_nan=True, 
    #     oe_transform=True, savedir=pltsavedir+"/clb")
    # calibrated_loci_1 = calb.perlabel_calibration_function(degree=5, num_bins=25, return_caliberated_matrix=True)
    # plt.close("all")
    # plt.style.use('default')
    
    agr = Agreement(loci_1, loci_2, pltsavedir+"/agr")
    to_report["per-label agreement"] = agr.per_label_agreement()
    to_report["general agreement"] = agr.general_agreement()
    to_report["general log(o/e) agreement"] = agr.general_OE_ratio()
    to_report["general Cohens Kappa score"] = agr.general_cohens_kappa()
    agr.plot_agreement()
    agr.plot_CK()
    agr.plot_OE()
    del agr
    plt.close("all")
    plt.style.use('default')

    vis = sankey(loci_1, loci_2, pltsavedir+"/snk")
    vis.sankey_diag()
    vis.heatmap()
    del vis
    plt.close('all')
    plt.style.use('default')
    return to_report

def full_reproducibility_report(replicate_1_dir, replicate_2_dir, pltsavedir):
    """
    get full reproducibility results
    including stepwise merging of labels through 
    post-clustering.
    """
    loci_1, loci_2 = intersect_parsed_posteriors(
        replicate_1_dir+"/parsed_posterior.csv", 
        replicate_2_dir+"/parsed_posterior.csv")

    num_labels = loci_1.shape[1]-3
    loci_1.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    loci_2.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]

    print('generating confmat 1')
    
    conf_mat = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=False)

    assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')

    loci_1_mnemon = read_mnemonics(replicate_1_dir+"/mnemonics.txt")
    loci_2_mnemon = read_mnemonics(replicate_2_dir+"/mnemonics.txt")
    
    mnemon1_dict = {}
    for i in loci_1_mnemon:
        mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]

    mnemon2_dict = {}
    for i in loci_2_mnemon:
        mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]
    
    for i in range(len(assignment_pairs)):
        assignment_pairs[i] = (mnemon1_dict[str(assignment_pairs[i][0])], mnemon2_dict[str(assignment_pairs[i][1])])
    print(assignment_pairs)

    loci_1, loci_2 = \
        connect_bipartite(loci_1, loci_2, assignment_pairs)
    
    print('connected barpartite')

    print('generating confmat 2')
    conf_mat = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=True)  

    mat_max = conf_mat.max(axis=1).max(axis=0)
    mat_min = conf_mat.min(axis=1).min(axis=0)

    conf_mat = (conf_mat - mat_min) / (mat_max - mat_min)
    distance_matrix = 1 - conf_mat
    linkage = hc.linkage(distance_matrix, method='average')
    
    if os.path.exists(pltsavedir)==False:
        os.mkdir(pltsavedir)

    cmap = sns.cm.rocket_r
    sns.clustermap(
        distance_matrix, row_linkage=linkage, 
        col_linkage=linkage, annot=True, cmap=cmap)

    if os.path.exists(pltsavedir+'/post_clustering/')==False:
        os.mkdir(pltsavedir+'/post_clustering/')

    plt.savefig('{}/post_clustering/clustermap.pdf'.format(pltsavedir), format='pdf')
    plt.savefig('{}/post_clustering/clustermap.svg'.format(pltsavedir), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')
     
    reports = {}
    reports[str(num_labels)] = report_reproducibility(
        loci_1, loci_2, 
        pltsavedir=pltsavedir+"/{}_labels".format(num_labels))

    # merging clusters one at a time
    merged_label_ID = {}
    labels = loci_1.iloc[:, 3:].columns
    for i in range(num_labels):
        merged_label_ID[i] = labels[i]

    for m in range(len(linkage)):
        to_be_merged = [
            merged_label_ID[int(linkage[m, 0])],
            merged_label_ID[int(linkage[m, 1])],
        ]

        merged_label_ID[num_labels + m] = str(
            merged_label_ID[int(linkage[m, 0])] + "+" + merged_label_ID[int(linkage[m, 1])]
        )

        loci_1[merged_label_ID[num_labels + m]] = \
            loci_1[to_be_merged[0]] + loci_1[to_be_merged[1]]
        loci_1 = loci_1.drop(to_be_merged, axis=1)

        loci_2[merged_label_ID[num_labels + m]] = \
            loci_2[to_be_merged[0]] + loci_2[to_be_merged[1]]
        loci_2 = loci_2.drop(to_be_merged, axis=1)

        reports[str((num_labels-1) - m)] = report_reproducibility(
            loci_1, loci_2, 
            pltsavedir=pltsavedir+"/{}_labels".format((num_labels-1) - m))

    nl = list(reports.keys())
    ys =[reports[k]["general agreement"] for k in nl]
    plt.bar(list(nl), list(ys), color="grey")
    plt.title("Post-clustering Progress")
    plt.xlabel("Number of Labels")
    plt.ylabel("General Agreement")
    plt.savefig('{}/post_clustering/Progress.pdf'.format(pltsavedir), format='pdf')
    plt.savefig('{}/post_clustering/Progress.svg'.format(pltsavedir), format='svg')
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

def run_single_reprod_analysis(input_dict):
    print("running type: {}".format(input_dict["runtype"]))
    full_reproducibility_report(input_dict["rep1_dir"], input_dict["rep2_dir"], input_dict["output_dir"])

def RUN_ALL_REPROD_ANALYSIS(runs_dir, CellType_list, output_dir, multi_p=True, type="segway", n_processors=8):
    """Given a directory containing all segway or chromHMM runs, 
    generates all orders of reproducibility analysis including 
    replicates, concatenated, and param_init modes. Stores all results in 
    output_dir"""
    
    ls = os.listdir(runs_dir)
    random_seeds = list(set(
        [s.split("_")[-1] for s in ls if "rs" in s]
        ))

    run_instances = {}
    
    for ct in CellType_list:
        ct_runs = {}
        ct_runs["replicates"] = [
            "{}/{}_rep1".format(runs_dir, ct), 
            "{}/{}_rep2".format(runs_dir, ct), 
            "{}/{}/rep1_vs_rep2/".format(output_dir, ct)]
        
        ct_runs["rep1_paraminit"] = [
            "{}/{}_rep1_{}".format(runs_dir, ct, random_seeds[0]), 
            "{}/{}_rep1_{}".format(runs_dir, ct, random_seeds[1]), 
            "{}/{}/rep1_paraminit/".format(output_dir, ct)]

        ct_runs["rep2_paraminit"] = [
            "{}/{}_rep2_{}".format(runs_dir, ct, random_seeds[0]), 
            "{}/{}_rep2_{}".format(runs_dir, ct, random_seeds[1]), 
            "{}/{}/rep2_paraminit/".format(output_dir, ct)]

        if type == "segway":
            ct_runs["concat"] = [
            "{}/{}_concat_rep1".format(runs_dir, ct), 
            "{}/{}_concat_rep2".format(runs_dir, ct), 
            "{}/{}/concatenated/".format(output_dir, ct)]

        elif type == "chmm":
            ### TO BE COMPLETED
            ct_runs["concat"] = []

        run_instances[ct] = ct_runs
    
    list_of_runs = []
    for k, v in run_instances.items():
        for kk in v.keys():
            list_of_runs.append(
                {"runtype": "{}_{}".format(k, kk),
                "rep1_dir": v[kk][0],
                "rep2_dir": v[kk][1],
                "output_dir": v[kk][2]
                }
            )
    
    if multi_p:
        p_obj = mp.Pool(n_processors)
        p_obj.map(run_single_reprod_analysis, list_of_runs)

    else:
        for r in list_of_runs:
            run_single_reprod_analysis(r)

"""when running the whole script from start to end to generate (and reproduce) results
remember to put label interpretation in try blocks (skippable) to prevent any kind of
dependency issue of segtools to cause issues in reproducibility of results"""

if __name__=="__main__":
    CellType_list = np.array(
        ['K562', 'MCF-7', 'GM12878', 'HeLa-S3', 'CD14-positive monocyte'])

    download_dir = 'files/'
    segway_dir = 'segway_runs/'
    res_dir = 'reprod_results/'

    if os.path.exists(res_dir) == False:
        os.mkdir(res_dir)

    print('list of target celltypes', CellType_list)
    existing_data = np.array(check_if_data_exists(CellType_list, download_dir))
    CellType_list = [CellType_list[i] for i in range(len(CellType_list)) if existing_data[i]==False]

    if len(CellType_list) != 0:
        download_encode_files(CellType_list, download_dir, "GRCh38")
    else:
        print('No download required!')

    CellType_list = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+ct)]

    # clean up potential space characters in directory names to prevent later issues
    for ct in CellType_list:
        if " " in ct:
            os.system("mv {} {}".format(
                ct.replace(' ', '\ '), ct.replace(" ", "_")
            ))

    CellType_list = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+ct)]
    create_trackname_assay_file(download_dir)

    assays = {}
    for ct in CellType_list:
        assays[ct] = read_list_of_assays(download_dir+ct)

    print(assays)

    # convert all bigwigs to bedgraphs (for segway)
    for k, v in assays.items():
        for t in v:
            Convert_all_BW2BG(download_dir+k+'/'+t)

    # metadata = read_metadata(download_dir)
    for c in CellType_list:
        gather_replicates(celltype_dir=download_dir+c)

    
    # download chromosome sizes file for hg38
    if os.path.exists(download_dir+"hg38.chrom.sizes") == False:
        sizes_url = 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes'
        sizes_file_dl_response = requests.get(sizes_url, allow_redirects=True)
        open(download_dir+"hg38.chrom.sizes", 'wb').write(sizes_file_dl_response.content)
        print('downloaded the hg38.chrom.sizes file')

    # check for existence of genomedata files
    gd_exists = []
    for ct in CellType_list:
        if os.path.exists(download_dir+ct+'/rep1.genomedata') == False or \
            os.path.exists(download_dir+ct+'/rep2.genomedata') == False:
            gd_exists.append(False)
        
        else:
            gd_exists.append(True)

    gd_to_create = [CellType_list[i] for i in range(len(CellType_list)) if gd_exists[i]==False]

    if len(gd_to_create) != 0:
        p_obj = mp.Pool(len(gd_to_create))
        p_obj.map(partial(
            create_genomedata, sequence_file=download_dir+"hg38.chrom.sizes"), 
            [download_dir + ct for ct in gd_to_create])
    
    if os.path.exists(segway_dir) == False:
        os.mkdir(segway_dir)

    # Run segway replicates     MP
    partial_runs_i = partial(
        RunParse_segway_replicates, output_dir=segway_dir, random_seed=73)
    p_obj = mp.Pool(int(len(CellType_list)/2))
    p_obj.map(partial_runs_i, [download_dir+ct for ct in CellType_list])

    # parse_posteriors 
    print('Checking for unparsed posteriors...')
    list_of_seg_runs = [
        d for d in os.listdir(segway_dir) if os.path.isdir(segway_dir+'/'+d)]
    print(list_of_seg_runs)
    for d in list_of_seg_runs:
        print('-Checking for {}  ...'.format(segway_dir+'/'+d+'/parsed_posterior.csv'))

        if os.path.exists(segway_dir+'/'+d+'/parsed_posterior.csv') == False:
            parse_posterior_results(segway_dir+'/'+d, 100, mp=False)

        else:
            print('-Exists!')

    print('All parsed!')




