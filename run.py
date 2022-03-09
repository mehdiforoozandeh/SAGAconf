from _utils import *
from _pipeline import *
from _reproducibility import *
from _interpret import *
from _visualize import *
from get_data import *
from matchlabels import *
import pandas as pd
import glob, os 
import multiprocessing as mp
from functools import partial

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

        nav_rep1_chmm[tfmd.loc['assay', 'rep1_spv']] = '{}/{}/{}.bedGraph'.format(
            celltype_dir, tr, tfmd.loc['accession', 'rep1_spv'])

        nav_rep2_chmm[tfmd.loc['assay', 'rep2_spv']] = '{}/{}/{}.bedGraph'.format(
            celltype_dir, tr, tfmd.loc['accession', 'rep2_spv'])

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

def RunParse_segway_replicates(celltype_dir, output_dir, random_seed=73):
    name_sig = celltype_dir.split('/')[-1]
    num_tracks = len(
        [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)])

    params_dict_1 = {
        "random_seed":random_seed, "track_weight":0.01,
        "stws":1, "ruler_scale":100, "prior_strength":1, "resolution":100, 
        "mini_batch_fraction":0.01, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
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
        "mini_batch_fraction":0.01, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
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
        "mini_batch_fraction":0.01, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
        "name_sig":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0]), 
        "genomedata_file":celltype_dir+'/{}.genomedata'.format(replicate_number), 
        "traindir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0])+'_train', 
        "posteriordir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[0])+'_posterior'
    }

    params_dict_2 = {
        "random_seed":random_seeds[1], "track_weight":0.01,
        "stws":1, "ruler_scale":100, "prior_strength":1, "resolution":100, 
        "mini_batch_fraction":0.01, "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
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

def concat_genomedata():
    pass

def RunParse_segway_concat():
        '''
    navigate fcoc files in rep1
    navigate fcoc files in rep2

    rename_chroms()
    virtualize_chroms()

    virtualize_sizes_file()

    make genomedata_concat
    make genomedata_rep1_renamed
    make genomedata_rep2_renamed

    define params_dict
    run concat seg + post clean up

    parse posterior1
    parse posterior2
    '''

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

def report_reproducibility(loci_1, loci_2, pltsavedir, general=True, num_bins=20, merge_cc_curves=False, full_report=True):
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

    try:
        agr = Agreement(loci_1, loci_2, pltsavedir+"/agreement")
        vis = sankey(loci_1, loci_2, pltsavedir+"/sankey")

        print('general agreement:    ', agr.general_agreement())
        print('general o/e ratio:    ', agr.general_OE_ratio(log_transform=False))
        print('general CK:    ', agr.general_cohens_kappa())

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
            repr = Reprodroducibility_vs_posterior(loci_1,loci_2, pltsavedir+"/calib",log_transform=False)

            cc.plot_curve(plot_general=general, merge_plots=merge_cc_curves)
            repr.per_label_count_independent(num_bins=num_bins)

            if general:
                repr.general_count_independent(num_bins=num_bins)
        
    except:
        print('skipped eval')


def post_clustering(loci_1, loci_2, pltsavedir, OE_transform=True):
    num_labels = loci_1.shape[1]-3

    print('generating confmat 1')
    
    conf_mat = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=False)
    
    print('generated confmat 1')

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

    print('generated confmat 1')

    mat_max = conf_mat.max(axis=1).max(axis=0)
    mat_min = conf_mat.min(axis=1).min(axis=0)

    conf_mat = (conf_mat - mat_min) / (mat_max - mat_min)
    distance_matrix = 1 - conf_mat
    linkage = hc.linkage(distance_matrix, method='average')

    if os.path.exists(pltsavedir)==False:
        os.mkdir(pltsavedir)

    if os.path.exists(pltsavedir+'/post_clustering/')==False:
        os.mkdir(pltsavedir+'/post_clustering/')

    # initial reproducibility without cluster merging
    if os.path.exists(pltsavedir+"/{}_labels".format(num_labels)) == False:
            os.mkdir(pltsavedir+"/{}_labels".format(num_labels))
    report_reproducibility(
            corrected_loci_1, corrected_loci_2, pltsavedir+"/{}_labels".format(num_labels), 
            general=False, num_bins=20, merge_cc_curves=False)

    for m in range(len(linkage)):
        # merging clusters one at a time
        to_be_merged = [
            "posterior{}".format(int(linkage[m, 0])), 
            "posterior{}".format(int(linkage[m, 1]))]

        corrected_loci_1['posterior{}'.format(num_labels + m)] = \
            corrected_loci_1[to_be_merged[0]] + corrected_loci_1[to_be_merged[1]]

        corrected_loci_1 = corrected_loci_1.drop(to_be_merged, axis=1)

        corrected_loci_2['posterior{}'.format(num_labels + m)] = \
            corrected_loci_2[to_be_merged[0]] + corrected_loci_2[to_be_merged[1]]

        corrected_loci_2 = corrected_loci_2.drop(to_be_merged, axis=1)

        if os.path.exists(pltsavedir+"/"+str(num_labels-(m+1))+'_labels') == False:
            os.mkdir(pltsavedir+"/"+str(num_labels-(m+1))+'_labels')

        report_reproducibility(
            corrected_loci_1, corrected_loci_2, pltsavedir+"/"+str(num_labels-(m+1))+'_labels', 
            general=False, num_bins=20, merge_cc_curves=False, full_report=False)
    
    c_grid = sns.clustermap(
        distance_matrix, row_linkage=linkage, 
        col_linkage=linkage, annot=True)

    plt.savefig('{}/post_clustering/clustermap.pdf'.format(pltsavedir), format='pdf')
    plt.savefig('{}/post_clustering/clustermap.svg'.format(pltsavedir), format='svg')
    plt.clf()


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
    p_obj = mp.Pool(len(CellType_list))
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

    for ct in CellType_list:
        loci_1, loci_2 = intersect_parsed_posteriors(
            segway_dir+ct+"_rep1/parsed_posterior.csv",
            segway_dir+ct+"_rep2/parsed_posterior.csv")

        loci_1 = loci_1.iloc[500000:1000000, :].reset_index(drop=True)
        loci_2 = loci_2.iloc[500000:1000000, :].reset_index(drop=True)
        print('loaded and intersected parsed posteriors for {}'.format(ct))
        print('starting reproducibility evaluation for {}'.format(ct))
        """EVALUATE REPRODUCIBILITY"""
        post_clustering(loci_1, loci_2, res_dir+ct, OE_transform=True)

    # Run segway param-init test     MP

    # partial_runs_ii = partial(
    #     RunParse_segway_param_init, 
    #     replicate_number = 'rep1', output_dir=segway_dir, random_seeds=[7, 5])

    # p_obj = mp.Pool(len(CellType_list))
    # p_obj.map(partial_runs_ii, [download_dir+ct for ct in CellType_list])

    # partial_runs_iii = partial(
    #     RunParse_segway_param_init, replicate_number = 'rep2', output_dir=segway_dir, 
    #     sizes_file=download_dir+"hg38.chrom.sizes", random_seeds=[7, 5])
        
    # p_obj = mp.Pool(len(CellType_list))
    # p_obj.map(partial_runs_iii, [download_dir+ct for ct in CellType_list])

    # # parse_posteriors 
    # print('Checking for unparsed posteriors...')
    # list_of_seg_runs = [
    #     d for d in os.listdir(segway_dir) if os.path.isdir(segway_dir+'/'+d)]
    # print(list_of_seg_runs)
    # for d in list_of_seg_runs:
    #     print('-Checking for {}  ...'.format(segway_dir+'/'+d+'/parsed_posterior.csv'))

    #     if os.path.exists(segway_dir+'/'+d+'/parsed_posterior.csv') == False:
    #         parse_posterior_results(segway_dir+'/'+d, 100, mp=False)

    #     else:
    #         print('-Exists!')

    # print('All parsed!')