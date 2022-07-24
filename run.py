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
from pseudoreplicate import *
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

def make_pseudo_replicates(celltype_dir, m_p=True):
    to_do_bams = []
    tracks = [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)]
    for tr in tracks:
        bams = [celltype_dir+'/'+tr+"/"+bm for bm in os.listdir(celltype_dir+'/'+tr) if ".bam" in bm]
        to_do_bams = to_do_bams + bams

    print(to_do_bams)
    if m_p:
        p_obj = mp.Pool(8)
        p_obj.map(psdrep_pipeline, to_do_bams)
    else:
        for bm in to_do_bams:
            psdrep_pipeline(bm)


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


def create_psdrep_genomedata(celltype_dir, sequence_file):
    '''
    find the bam name corresponding to rep1/rep2
    find psdrep* files for each track and each replicate

    create genomedata file_rep1 in celltype_dir
    create genomedata file_rep2 in celltype_dir
    '''

    print('creating genomedata files for {}...'.format(celltype_dir))

    assaylist = [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)]# and tr in essential_tracks]
    navigate = {
        "rep1_psdrep1":[],
        "rep1_psdrep2":[],
        "rep2_psdrep1":[],
        "rep2_psdrep2":[]
    }

    for tr in assaylist:
        tfmd = pd.read_csv(celltype_dir+'/'+tr+'/track_files_metadata.csv')
        tfmd.index = list(tfmd['Unnamed: 0'])
        tfmd = tfmd.drop('Unnamed: 0', axis=1) 
        navigate["rep1_psdrep1"].append(
            [tfmd.loc['assay', 'rep1_alig'], "rep1", str(tfmd.loc['accession', 'rep1_alig'])+"_psdrep1.fc.signal.bedGraph"]
        )
        navigate["rep1_psdrep2"].append(
            [tfmd.loc['assay', 'rep1_alig'], "rep1", str(tfmd.loc['accession', 'rep1_alig'])+"_psdrep2.fc.signal.bedGraph"]
        )
        navigate["rep2_psdrep1"].append(
            [tfmd.loc['assay', 'rep2_alig'], "rep2", str(tfmd.loc['accession', 'rep2_alig'])+"_psdrep1.fc.signal.bedGraph"]
        )
        navigate["rep2_psdrep2"].append(
            [tfmd.loc['assay', 'rep2_alig'], "rep2", str(tfmd.loc['accession', 'rep2_alig'])+"_psdrep2.fc.signal.bedGraph"]
        )
    tracklist_rep1psd1 = ''
    for e in navigate["rep1_psdrep1"]:
        tracklist_rep1psd1 = tracklist_rep1psd1 + '-t {}={} '.format(e[0], "{}/{}/{}".format(celltype_dir, e[0], e[2]))

    tracklist_rep1psd2 = ''
    for e in navigate["rep1_psdrep2"]:
        tracklist_rep1psd2 = tracklist_rep1psd2 + '-t {}={} '.format(e[0], "{}/{}/{}".format(celltype_dir, e[0], e[2]))

    tracklist_rep2psd1 = ''
    for e in navigate["rep2_psdrep1"]:
        tracklist_rep2psd1 = tracklist_rep2psd1 + '-t {}={} '.format(e[0], "{}/{}/{}".format(celltype_dir, e[0], e[2]))

    tracklist_rep2psd2 = ''
    for e in navigate["rep2_psdrep2"]:
        tracklist_rep2psd2 = tracklist_rep2psd2 + '-t {}={} '.format(e[0], "{}/{}/{}".format(celltype_dir, e[0], e[2]))

    os.system(
        'genomedata-load -s {} --sizes {} --verbose {}.genomedata'.format(
            sequence_file, tracklist_rep1psd1, celltype_dir+'/rep1_psdrep1'))

    os.system(
        'genomedata-load -s {} --sizes {} --verbose {}.genomedata'.format(
            sequence_file, tracklist_rep1psd2, celltype_dir+'/rep1_psdrep2'))

    os.system(
        'genomedata-load -s {} --sizes {} --verbose {}.genomedata'.format(
            sequence_file, tracklist_rep2psd1, celltype_dir+'/rep2_psdrep1'))

    os.system(
        'genomedata-load -s {} --sizes {} --verbose {}.genomedata'.format(
            sequence_file, tracklist_rep2psd2, celltype_dir+'/rep2_psdrep2'))

def create_concat_psdrep_genomedata(celltype_dir, sequence_file):
    # navigate all psdfiles for all tracks
    print('creating genomedata files for {}...'.format(celltype_dir))

    assaylist = [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)]# and tr in essential_tracks]
    navigate = {
        "rep1_psdrep1":[],
        "rep1_psdrep2":[],
        "rep2_psdrep1":[],
        "rep2_psdrep2":[]
    }

    for tr in assaylist:
        tfmd = pd.read_csv(celltype_dir+'/'+tr+'/track_files_metadata.csv')
        tfmd.index = list(tfmd['Unnamed: 0'])
        tfmd = tfmd.drop('Unnamed: 0', axis=1)

        navigate["rep1_psdrep1"].append(
            [tfmd.loc['assay', 'rep1_alig'], "rep1", str(tfmd.loc['accession', 'rep1_alig'])+"_psdrep1.fc.signal.bedGraph"]
        )
        navigate["rep1_psdrep2"].append(
            [tfmd.loc['assay', 'rep1_alig'], "rep1", str(tfmd.loc['accession', 'rep1_alig'])+"_psdrep2.fc.signal.bedGraph"]
        )
        navigate["rep2_psdrep1"].append(
            [tfmd.loc['assay', 'rep2_alig'], "rep2", str(tfmd.loc['accession', 'rep2_alig'])+"_psdrep1.fc.signal.bedGraph"]
        )
        navigate["rep2_psdrep2"].append(
            [tfmd.loc['assay', 'rep2_alig'], "rep2", str(tfmd.loc['accession', 'rep2_alig'])+"_psdrep2.fc.signal.bedGraph"]
        )

    # chrom-rename for all of them (rename chroms based on psd number)
    for i in range(len(assaylist)):
        if os.path.exists(celltype_dir + "/" + navigate["rep1_psdrep1"][i][0] + "/" + navigate["rep1_psdrep1"][i][2].replace(".bedGraph", "_concat.bedGraph")) == False:
            concat_rename_chroms(
                celltype_dir + "/" + navigate["rep1_psdrep1"][i][0] + "/" + navigate["rep1_psdrep1"][i][2], 
                celltype_dir + "/" + navigate["rep1_psdrep1"][i][0] + "/" + navigate["rep1_psdrep1"][i][2].replace(".bedGraph", "_concat.bedGraph"), 1)

        if os.path.exists(celltype_dir + "/" + navigate["rep1_psdrep2"][i][0] + "/" + navigate["rep1_psdrep2"][i][2].replace(".bedGraph", "_concat.bedGraph")) == False:
            concat_rename_chroms(
                celltype_dir + "/" + navigate["rep1_psdrep2"][i][0] + "/" + navigate["rep1_psdrep2"][i][2], 
                celltype_dir + "/" + navigate["rep1_psdrep2"][i][0] + "/" + navigate["rep1_psdrep2"][i][2].replace(".bedGraph", "_concat.bedGraph"), 2)

        if os.path.exists(celltype_dir + "/" + navigate["rep2_psdrep1"][i][0] + "/" + navigate["rep2_psdrep1"][i][2].replace(".bedGraph", "_concat.bedGraph")) == False:
            concat_rename_chroms(
                celltype_dir + "/" + navigate["rep2_psdrep1"][i][0] + "/" + navigate["rep2_psdrep1"][i][2], 
                celltype_dir + "/" + navigate["rep2_psdrep1"][i][0] + "/" + navigate["rep2_psdrep1"][i][2].replace(".bedGraph", "_concat.bedGraph"), 1)

        if os.path.exists(celltype_dir + "/" + navigate["rep2_psdrep2"][i][0] + "/" + navigate["rep2_psdrep2"][i][2].replace(".bedGraph", "_concat.bedGraph")) == False:
            concat_rename_chroms(
                celltype_dir + "/" + navigate["rep2_psdrep2"][i][0] + "/" + navigate["rep2_psdrep2"][i][2], 
                celltype_dir + "/" + navigate["rep2_psdrep2"][i][0] + "/" + navigate["rep2_psdrep2"][i][2].replace(".bedGraph", "_concat.bedGraph"), 2)

        # virtualize bedgraphs
        if os.path.exists(celltype_dir + "/" + navigate["rep1_psdrep1"][i][0] + "/rep1psd_concatenated.bedGraph") == False:
            concat_virtualize_chroms(
                celltype_dir + "/" + navigate["rep1_psdrep1"][i][0] + "/" + navigate["rep1_psdrep1"][i][2].replace(".bedGraph", "_concat.bedGraph"),
                celltype_dir + "/" + navigate["rep1_psdrep2"][i][0] + "/" + navigate["rep1_psdrep2"][i][2].replace(".bedGraph", "_concat.bedGraph"),
                celltype_dir + "/" + navigate["rep1_psdrep1"][i][0] + "/rep1psd_concatenated.bedGraph")

        if os.path.exists(celltype_dir + "/" + navigate["rep2_psdrep1"][i][0] + "/rep2psd_concatenated.bedGraph") == False:
            concat_virtualize_chroms(
                celltype_dir + "/" + navigate["rep2_psdrep1"][i][0] + "/" + navigate["rep2_psdrep1"][i][2].replace(".bedGraph", "_concat.bedGraph"),
                celltype_dir + "/" + navigate["rep2_psdrep2"][i][0] + "/" + navigate["rep2_psdrep2"][i][2].replace(".bedGraph", "_concat.bedGraph"),
                celltype_dir + "/" + navigate["rep2_psdrep1"][i][0] + "/rep2psd_concatenated.bedGraph")

    # virtualize chrsz file
    if os.path.exists(sequence_file.replace(".sizes", "_concatenated.sizes")) == False:
        concat_virtualize_include_file(
            sequence_file, sequence_file.replace(".sizes", "_concatenated.sizes"), 
            write_separate=True)

    # create 6 genomedata files per celltype:  
    # concatenated_rep1psd
    rep1psd_tracklist = ''

    # concat_rep1psd1
    rep1psd1_tracklist = ''
    
    # concat_rep1psd2
    rep1psd2_tracklist = ''

    # concatenated_rep2psd
    rep2psd_tracklist = ''

    # concat_rep2psd1
    rep2psd1_tracklist = ''
    
    # concat_rep2psd2
    rep2psd2_tracklist = ''

    for i in range(len(assaylist)):
        rep1psd_tracklist = rep1psd_tracklist + '-t {}={} '.format(
            navigate["rep1_psdrep1"][i][0], celltype_dir + "/" + navigate["rep1_psdrep1"][i][0] + "/rep1psd_concatenated.bedGraph")

        rep1psd1_tracklist = rep1psd1_tracklist + '-t {}={} '.format(
            navigate["rep1_psdrep1"][i][0], celltype_dir + "/" + navigate["rep1_psdrep1"][i][0] + "/" + navigate["rep1_psdrep1"][i][2].replace(".bedGraph", "_concat.bedGraph"))

        rep1psd2_tracklist = rep1psd2_tracklist + '-t {}={} '.format(
            navigate["rep1_psdrep1"][i][0], celltype_dir + "/" + navigate["rep1_psdrep2"][i][0] + "/" + navigate["rep1_psdrep2"][i][2].replace(".bedGraph", "_concat.bedGraph"))
        
        # rep2psd_tracklist = rep2psd_tracklist + '-t {}={} '.format(
        #     navigate["rep2_psdrep1"][i][0], celltype_dir + "/" + navigate["rep2_psdrep1"][i][0] + "/rep2psd_concatenated.bedGraph")

        # rep2psd1_tracklist = rep2psd1_tracklist + '-t {}={} '.format(
        #     navigate["rep2_psdrep1"][i][0], celltype_dir + "/" + navigate["rep2_psdrep1"][i][0] + "/" + navigate["rep2_psdrep1"][i][2].replace(".bedGraph", "_concat.bedGraph"))
            
        # rep2psd2_tracklist = rep2psd2_tracklist + '-t {}={} '.format(
        #     navigate["rep2_psdrep2"][i][0], celltype_dir + "/" + navigate["rep2_psdrep2"][i][0] + "/" + navigate["rep2_psdrep2"][i][2].replace(".bedGraph", "_concat.bedGraph"))
        
    if os.path.exists(celltype_dir + "/" + "rep1psd_concatenated.genomedata") == False:
        os.system('genomedata-load -s {} --sizes {} --verbose {}'.format(
            sequence_file.replace(".sizes", "_concatenated.sizes"), 
            rep1psd_tracklist, 
            celltype_dir + "/" + "rep1psd_concatenated.genomedata"))

    if os.path.exists(celltype_dir + "/" + "rep1psd1_concat.genomedata") == False:
        os.system('genomedata-load -s {} --sizes {} --verbose {}'.format(
            sequence_file.replace(".sizes", "_VirtRep1_concatenated.sizes"), 
            rep1psd1_tracklist, 
            celltype_dir + "/" + "rep1psd1_concat.genomedata"))

    if os.path.exists(celltype_dir + "/" + "rep1psd2_concat.genomedata") == False:
        os.system('genomedata-load -s {} --sizes {} --verbose {}'.format(
            sequence_file.replace(".sizes", "_VirtRep2_concatenated.sizes"), 
            rep1psd2_tracklist, 
            celltype_dir + "/" + "rep1psd2_concat.genomedata"))
    
    #========================================================================================#

    # if os.path.exists(celltype_dir + "/" + "rep2psd_concatenated.genomedata") == False:
    #     os.system('genomedata-load -s {} --sizes {} --verbose {}'.format(
    #         sequence_file.replace(".sizes", "_concatenated.sizes"), 
    #         rep2psd_tracklist, 
    #         celltype_dir + "/" + "rep2psd_concatenated.genomedata"))

    # if os.path.exists(celltype_dir + "/" + "rep2psd1_concat.genomedata") == False:
    #     os.system('genomedata-load -s {} --sizes {} --verbose {}'.format(
    #         sequence_file.replace(".sizes", "_VirtRep1_concatenated.sizes"), 
    #         rep2psd1_tracklist, 
    #         celltype_dir + "/" + "rep2psd1_concat.genomedata"))

    # if os.path.exists(celltype_dir + "/" + "rep2psd2_concat.genomedata") == False:
    #     os.system('genomedata-load -s {} --sizes {} --verbose {}'.format(
    #         sequence_file.replace(".sizes", "_VirtRep2_concatenated.sizes"), 
    #         rep2psd2_tracklist, 
    #         celltype_dir + "/" + "rep2psd2_concat.genomedata"))


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

def RunParse_segway_psdreps(celltype_dir, output_dir, random_seed=73):
    name_sig = celltype_dir.split('/')[-1]
    num_tracks = len(
        [tr for tr in os.listdir(celltype_dir) if os.path.isdir(celltype_dir+'/'+tr)])

    concat_param_dict = {
        "random_seed":random_seed, "track_weight":0.01,"stws":1, "ruler_scale":100, 
        "prior_strength":1, "resolution":100, "mini_batch_fraction":0.01,
        "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
        "name_sig_concat":output_dir+name_sig+'_rep1psd_concat', 
        "name_sig_rep1":output_dir+name_sig+'_concat_rep1psd1', 
        "name_sig_rep2":output_dir+name_sig+'_concat_rep1psd2',

        "genomedata_file_concat":celltype_dir+'/rep1psd_concatenated.genomedata',  
        "genomedata_file_rep1":celltype_dir+'/rep1psd1_concat.genomedata', 
        "genomedata_file_rep2":celltype_dir+'/rep1psd2_concat.genomedata',

        "traindir":output_dir+name_sig+'_rep1psd_concat_train', 
        "posteriordir_rep1":output_dir+name_sig+'_rep1psd_concat_posterior_1',
        "posteriordir_rep2":output_dir+name_sig+'_rep1psd_concat_posterior_2'}

    if os.path.exists(concat_param_dict['name_sig_rep1']) == False or \
         os.path.exists(concat_param_dict['name_sig_rep2']) == False:
        print('Running Segway concatenated celltype {}...'.format(celltype_dir))
        concat_segwayrun_and_postprocess(concat_param_dict)

    # concat_param_dict = {
    #     "random_seed":random_seed, "track_weight":0.01,"stws":1, "ruler_scale":100, 
    #     "prior_strength":1, "resolution":100, "mini_batch_fraction":0.03,
    #     "num_labels": 10 + (2 * int(np.sqrt(num_tracks))), 
    #     "name_sig_concat":output_dir+name_sig+'_concat', 
    #     "name_sig_rep1":output_dir+name_sig+'_concat_rep1', 
    #     "name_sig_rep2":output_dir+name_sig+'_concat_rep2',
    #     "genomedata_file_concat":celltype_dir+'/concatenated.genomedata',  
    #     "genomedata_file_rep1":celltype_dir+'/concat_rep1.genomedata', 
    #     "genomedata_file_rep2":celltype_dir+'/concat_rep2.genomedata',
    #     "traindir":output_dir+name_sig+'_concat_train', 
    #     "posteriordir_rep1":output_dir+name_sig+'_concat_posterior_1',
    #     "posteriordir_rep2":output_dir+name_sig+'_concat_posterior_2'}

    # if os.path.exists(concat_param_dict['name_sig_rep1']) == False or \
    #      os.path.exists(concat_param_dict['name_sig_rep2']) == False:
    #     print('Running Segway concatenated celltype {}...'.format(celltype_dir))
    #     concat_segwayrun_and_postprocess(concat_param_dict)


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
        "traindir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[1])+'_train', 
        "posteriordir":output_dir+name_sig+'_{}_rs{}'.format(replicate_number, random_seeds[1])+'_posterior'
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
    df1 = pd.read_csv(parsed_df_dir_1).drop("Unnamed: 0", axis=1)
    df1.iloc[:, 3:] = df1.iloc[:, 3:].astype("float16")
    df2 = pd.read_csv(parsed_df_dir_2).drop("Unnamed: 0", axis=1)
    df2.iloc[:, 3:] = df2.iloc[:, 3:].astype("float16")
    
    # to handle concat indexing
    if "_1" in df1.iloc[0, 0] or "_2" in df1.iloc[0, 0]:
        chrdf1 = list(df1.chr)
        for i in range(len(chrdf1)):
            if "_1" in chrdf1[i]:
                chrdf1[i] = chrdf1[i].replace("_1", "")
            elif "_2" in chrdf1[i]:
                chrdf1[i] = chrdf1[i].replace("_2", "")
        df1.chr = np.array(chrdf1)

    if "_1" in df2.iloc[0, 0] or "_2" in df2.iloc[0, 0]:
        chrdf2 = list(df2.chr)
        for i in range(len(chrdf2)):
            if "_2" in chrdf2[i]:
                chrdf2[i] = chrdf2[i].replace("_2", "")
            elif "_1" in chrdf2[i]:
                chrdf2[i] = chrdf2[i].replace("_1", "")
        df2.chr = np.array(chrdf2)
    
    # if "_1" in df1.iloc[0, 0] or "_2" in df1.iloc[0, 0]:
    #     for i in range(len(df1)):
    #         if "_1" in df1.iloc[i, 0]:
    #             df1.iloc[i, 0] = df1.iloc[i, 0].replace("_1","")
    #         elif "_2" in df1.iloc[i, 0]:
    #             df1.iloc[i, 0] = df1.iloc[i, 0].replace("_2","")

    # if "_1" in df2.iloc[0, 0] or "_2" in df2.iloc[0, 0]:
    #     for i in range(len(df2)):
    #         if "_1" in df2.iloc[i, 0]:
    #             df2.iloc[i, 0] = df2.iloc[i, 0].replace("_1","")
    #         elif "_2" in df2.iloc[i, 0]:
    #             df2.iloc[i, 0] = df2.iloc[i, 0].replace("_2","")

    intersect = pd.merge(
        df1, 
        df2, 
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
        if int(df["old"][0]) == 0:
            mnemon.append(str(df["old"][i])+"_"+df["new"][i])
        elif int(df["old"][0]) == 1:
            mnemon.append(str(int(df["old"][i])-1)+"_"+df["new"][i])
    return mnemon

def get_short_report(replicate_1_dir, replicate_2_dir, outdir, type="chmm"):
    """
    copy dir1 emission to outdir
    copy dir2 emission to outdir

    read parsed_posterior1
    read parsed_posterior2

    get raw overlap CM
    get matched overlap CM

    get agreements rep1 vs rep2
    get agreements rep2 vs rep1
    """

    if os.path.exists(outdir)==False:
        os.mkdir(outdir)

    ls1 = os.listdir("/".join(replicate_1_dir.split("/")[:-1]))
    ls2 = os.listdir("/".join(replicate_2_dir.split("/")[:-1]))

    if type == "chmm":
        for l in ls1:
            if "emissions" in l and ".svg" in l:
                os.system("cp {}/{} {}/emissions_1.svg".format("/".join(replicate_1_dir.split("/")[:-1]), l, outdir))

        for l in ls2:
            if "emissions" in l and ".svg" in l:
                os.system("cp {}/{} {}/emissions_2.svg".format("/".join(replicate_2_dir.split("/")[:-1]), l, outdir))

    elif type == "segw":
        os.system("cp {}/gmtk_parameters/gmtk_parameters.stats.slide.png {}/emissions_1.png".format(
            "/".join(replicate_1_dir.split("/")[:-1]), outdir))
        os.system("cp {}/gmtk_parameters/gmtk_parameters.stats.slide.png {}/emissions_2.png".format(
            "/".join(replicate_2_dir.split("/")[:-1]), outdir))

    loci_1, loci_2 = intersect_parsed_posteriors(
        replicate_1_dir, 
        replicate_2_dir)
    
    num_labels = loci_1.shape[1]-3
    loci_1.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    loci_2.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]

    confmat_raw = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=False, symmetric=False)
    plot_heatmap(confmat_raw, outdir, type="RawOverlapMatrix", columns=[list(loci_1.columns[3:]), list(loci_2.columns[3:])])

    confmat_OE = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=False)    
    plot_heatmap(confmat_OE, outdir, type="log_OE_OverlapMatrix", columns=[list(loci_1.columns[3:]), list(loci_2.columns[3:])])

    assignment_pairs = Hungarian_algorithm(confmat_OE, conf_or_dis='conf')
    new_columns = ["{}|{}".format(c[0], c[1]) for c in assignment_pairs]

    if os.path.exists(
        "/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics_rep1.txt") and os.path.exists(
        "/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics_rep2.txt"):

        print("reading concat mnemonics")
        loci_1_mnemon = read_mnemonics("/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics_rep1.txt")
        loci_2_mnemon = read_mnemonics("/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics_rep2.txt")
    else:
        print("reading mnemonics")
        loci_1_mnemon = read_mnemonics("/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics.txt")
        loci_2_mnemon = read_mnemonics("/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics.txt")

    mnemon1_dict = {}
    for i in loci_1_mnemon:
        mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]

    mnemon2_dict = {}
    for i in loci_2_mnemon:
        mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]

    #handle missing mnemonics
    for i in range(num_labels):
        if str(i) not in mnemon1_dict.keys():
            mnemon1_dict[str(i)] = str(i)
        if str(i) not in mnemon2_dict.keys():
            mnemon2_dict[str(i)] = str(i)
        

    for i in range(len(assignment_pairs)):
        assignment_pairs[i] = (mnemon1_dict[str(assignment_pairs[i][0])], mnemon2_dict[str(assignment_pairs[i][1])])
    print(assignment_pairs)

    loci_1, loci_2 = \
        connect_bipartite(loci_1, loci_2, assignment_pairs)

    print('connected barpartite')

    print('generating confmat 2')
    confmat_OE_matched = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=False) 

    confmat_raw_matched = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=False, symmetric=False)  

    
    plot_heatmap(confmat_raw_matched, outdir, type="RawOverlapMatrix_matched", columns=[new_columns,new_columns])
    plot_heatmap(confmat_OE_matched, outdir, type="logOEOverlapMatrix_matched", columns=[new_columns,new_columns])

    agreement_report_1 = {}
    agreement_report_2 = {}
    
    print("calculating agr 1vs2")
    agr1 = Agreement(loci_1, loci_2, outdir+"/agr1", assume_uniform_coverage=False)
    agreement_report_1["PL_agr"] = agr1.per_label_agreement()
    agreement_report_1["PL_OE"] = agr1.per_label_OE_ratio()
    agreement_report_1["PL_CK"] = agr1.per_label_cohens_kappa()
    agreement_report_1["G_agr"] = agr1.general_agreement()
    agreement_report_1["G_OE"] = agr1.general_OE_ratio()
    agreement_report_1["G_CK"] = agr1.general_cohens_kappa()

    print("calculating agr 2vs1")
    agr2 = Agreement(loci_2, loci_1, outdir+"/agr2", assume_uniform_coverage=False)
    agreement_report_2["PL_agr"] = agr2.per_label_agreement()
    agreement_report_2["PL_OE"] = agr2.per_label_OE_ratio()
    agreement_report_2["PL_CK"] = agr2.per_label_cohens_kappa()
    agreement_report_2["G_agr"] = agr2.general_agreement()
    agreement_report_2["G_OE"] = agr2.general_OE_ratio()
    agreement_report_2["G_CK"] = agr2.general_cohens_kappa()

    print("plotting bidirectional agreements")
    plot_bidir_bar_chart(agreement_report_1["PL_agr"], agreement_report_2["PL_agr"], type="Raw_Agreement", savedir=outdir)
    plot_bidir_bar_chart(agreement_report_1["PL_OE"], agreement_report_2["PL_OE"], type="logOE_Agreement", savedir=outdir)
    plot_bidir_bar_chart(agreement_report_1["PL_CK"], agreement_report_2["PL_CK"], type="Cohens_Kappa", savedir=outdir)


def report_reproducibility(loci_1, loci_2, pltsavedir, cc_calb=True):
    """
    get basic reproducibility results for a pair of 
    experiments (two locis)
    """
    if os.path.exists(pltsavedir) == False:
        os.mkdir(pltsavedir)

    if os.path.exists(pltsavedir+"/agr") == False:
        os.mkdir(pltsavedir+"/agr")
    
    if os.path.exists(pltsavedir+"/snk") == False:
        os.mkdir(pltsavedir+"/snk")
    
    # plot posterior distribution

    to_report = {}

    agr = Agreement(loci_1, loci_2, pltsavedir+"/agr", assume_uniform_coverage=False)
    to_report["per-label agreement"] = agr.per_label_agreement()
    to_report["per-label log(o/e) agreement"] = agr.per_label_OE_ratio()
    to_report["per-label  Cohens Kappa score"] = agr.per_label_cohens_kappa()

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
    
    ##============================================================================================================##

    # if cc_calb:
    #     try:
    #         if os.path.exists(pltsavedir+"/pstdist_rep1") == False:
    #             os.mkdir(pltsavedir+"/pstdist_rep1")
    #         pstdist = Posterior_dist(loci_1, pltsavedir+"/pstdist_rep1")
    #         pstdist.plot_posterior_histogram()

    #         if os.path.exists(pltsavedir+"/pstdist_rep2") == False:
    #             os.mkdir(pltsavedir+"/pstdist_rep2")
    #         pstdist = Posterior_dist(loci_2, pltsavedir+"/pstdist_rep2")
    #         pstdist.plot_posterior_histogram()

    #         del pstdist
            
    #     except:
    #         pass

    #     if os.path.exists(pltsavedir+"/cc") == False:
    #         os.mkdir(pltsavedir+"/cc")

    #     try:
    #         cc = correspondence_curve(loci_1, loci_2, pltsavedir+"/cc")
    #         cc.plot_curve(plot_general=False, merge_plots=False)
    #         del cc
    #         plt.close("all")
    #         plt.style.use('default')

    #     except:
    #         print("could not generate corresp. curve")

    #     try:
    #         if os.path.exists(pltsavedir+"/clb_1") == False:
    #             os.mkdir(pltsavedir+"/clb_1")
    #         calb = posterior_calibration(
    #             loci_1, loci_2, log_transform=False, ignore_overconf=False, filter_nan=True, 
    #             oe_transform=True, savedir=pltsavedir+"/clb_1")
    #         calibrated_loci_1 = calb.perlabel_calibration_function(
    #             degree=5, num_bins=25, return_caliberated_matrix=True, scale_columnwise=True)
            
    #         plt.close("all")
    #         plt.style.use('default')

    #         if os.path.exists(pltsavedir+"/tss_rep1") == False:
    #             os.mkdir(pltsavedir+"/tss_rep1")
    #         TSS_obj = TSS_enrichment(calibrated_loci_1, TSSdir="RefSeqTSS.hg38.txt", savedir=pltsavedir+"/tss_rep1")
    #         TSS_obj.tss_enrich(m_p=False)
    #         TSS_obj.tss_enrich_vs_repr()

    #     except:
    #         print("could not generate calibrations and TSS enrichment")
        
        ##============================================================================================================##
        
        # try:
        #     if os.path.exists(pltsavedir+"/clb_2") == False:
        #         os.mkdir(pltsavedir+"/clb_2")
        #     calb = posterior_calibration(
        #         loci_2, loci_1, log_transform=False, ignore_overconf=False, filter_nan=True, 
        #         oe_transform=True, savedir=pltsavedir+"/clb_2")
        #     calibrated_loci_2 = calb.perlabel_calibration_function(
        #         degree=5, num_bins=25, return_caliberated_matrix=True, scale_columnwise=True)
            
        #     plt.close("all")
        #     plt.style.use('default')

        #     if os.path.exists(pltsavedir+"/tss_rep2") == False:
        #         os.mkdir(pltsavedir+"/tss_rep2")
        #     TSS_obj = TSS_enrichment(calibrated_loci_2, TSSdir="RefSeqTSS.hg38.txt", savedir=pltsavedir+"/tss_rep2")
        #     TSS_obj.tss_enrich(m_p=False)
        #     TSS_obj.tss_enrich_vs_repr()
        # except:
        #     print("could not generate calibrations and TSS enrichment")

    return to_report

def full_reproducibility_report(replicate_1_dir, replicate_2_dir, pltsavedir, run_on_subset=False):
    """
    get full reproducibility results
    including stepwise merging of labels through 
    post-clustering.
    """

    print("loading and intersecting")
    loci_1, loci_2 = intersect_parsed_posteriors(
        replicate_1_dir, 
        replicate_2_dir)

    if run_on_subset:
        loci_1 = loci_1.iloc[[i for i in range(0, len(loci_1), 100)], :].reset_index(drop=True)
        loci_2 = loci_2.iloc[[i for i in range(0, len(loci_2), 100)], :].reset_index(drop=True)

    print("the shapes of the input matrices are: {}, {}".format(str(loci_1.shape), str(loci_2.shape)))

    num_labels = loci_1.shape[1]-3
    loci_1.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    loci_2.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    
    print('generating confmat 1')
    
    conf_mat = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=False)

    assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')

    new_columns = ["{}|{}".format(c[0], c[1]) for c in assignment_pairs]

    if os.path.exists(
        "/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics_rep1.txt") and os.path.exists(
        "/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics_rep2.txt"):

        print("reading concat mnemonics")
        loci_1_mnemon = read_mnemonics("/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics_rep1.txt")
        loci_2_mnemon = read_mnemonics("/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics_rep2.txt")
    else:
        print("reading mnemonics")
        loci_1_mnemon = read_mnemonics("/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics.txt")
        loci_2_mnemon = read_mnemonics("/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics.txt")

    mnemon1_dict = {}
    for i in loci_1_mnemon:
        mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]

    mnemon2_dict = {}
    for i in loci_2_mnemon:
        mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]

    #handle missing mnemonics
    for i in range(num_labels):
        if str(i) not in mnemon1_dict.keys():
            mnemon1_dict[str(i)] = str(i)
        if str(i) not in mnemon2_dict.keys():
            mnemon2_dict[str(i)] = str(i)
        

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

    vis = sankey(loci_1, loci_2, pltsavedir)
    vis.heatmap(new_columns)
     
    reports = {}
    reports[str(num_labels)] = report_reproducibility(
        loci_1, loci_2, 
        pltsavedir=pltsavedir+"/{}_labels".format(num_labels), cc_calb=True)

    coverage1= {}
    coverage2= {}

    merge_track1 = []
    merge_track2 = []

    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    coverage1[str(num_labels)] = MAP1.value_counts()/len(MAP1)

    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
    coverage2[str(num_labels)] = MAP2.value_counts()/len(MAP2)

    merge_track1.append(MAP1)
    merge_track2.append(MAP2)

    # to avoid MAP changes during the merging process, disregard probability granularity and convert 
    # to hard zero vs one matrix

    for c in loci_1.columns[3:]:
        loci_1.loc[MAP1==c, c] = 1
        loci_1.loc[MAP1!=c, c] = 0
    
    for c in loci_2.columns[3:]:
        loci_2.loc[MAP2==c, c] = 1
        loci_2.loc[MAP2!=c, c] = 0

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

        MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
        coverage1[str((num_labels-1) - m)] = MAP1.value_counts()/len(MAP1)
        merge_track1.append(MAP1)

        MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
        coverage2[str((num_labels-1) - m)] = MAP2.value_counts()/len(MAP2)
        merge_track2.append(MAP2)
        
        reports[str((num_labels-1) - m)] = report_reproducibility(
            loci_1, loci_2, 
            pltsavedir=pltsavedir+"/{}_labels".format((num_labels-1) - m), 
            cc_calb=False)

    merge_track1 = pd.concat(merge_track1, axis=1)
    merge_track1.columns = ["{}_labels".format(16-lm) for lm in range(16)] 

    merge_track2 = pd.concat(merge_track2, axis=1)
    merge_track2.columns = ["{}_labels".format(16-lm) for lm in range(16)] 

    merge_track1.to_csv('{}/post_clustering/merged_MAP1.csv'.format(pltsavedir))
    os.system("gzip {}/post_clustering/merged_MAP1.csv".format(pltsavedir))

    merge_track2.to_csv('{}/post_clustering/merged_MAP2.csv'.format(pltsavedir))
    os.system("gzip {}/post_clustering/merged_MAP2.csv".format(pltsavedir))

    nl = list(reports.keys())
    ys =[reports[k]["general agreement"] for k in nl]
    plt.bar(list(nl), list(ys), color="grey")
    plt.title("Post-clustering Progress")
    plt.xlabel("Number of Labels")
    plt.ylabel("General Agreement")
    plt.savefig('{}/post_clustering/agr_Progress.pdf'.format(pltsavedir), format='pdf')
    plt.savefig('{}/post_clustering/agr_Progress.svg'.format(pltsavedir), format='svg')
    plt.clf()

    with open('{}/post_clustering/agr_Progress.txt'.format(pltsavedir), "w") as wrfile:
        # xlabel, ylabel, bars, heights
        wrfile.write(
            "{}\n{}\n{}\n{}".format(
                "Number of Labels", "General Agreement",
                str(nl), str(ys)
            )
        )

    ys =[reports[k]["general log(o/e) agreement"] for k in nl]
    plt.bar(list(nl), list(ys), color="grey")
    plt.title("Post-clustering Progress")
    plt.xlabel("Number of Labels")
    plt.ylabel("log(O/E) Agreement")
    plt.savefig('{}/post_clustering/oe_agr_Progress.pdf'.format(pltsavedir), format='pdf')
    plt.savefig('{}/post_clustering/oe_agr_Progress.svg'.format(pltsavedir), format='svg')
    plt.clf()

    with open('{}/post_clustering/oe_agr_Progress.txt'.format(pltsavedir), "w") as wrfile:
        # xlabel, ylabel, bars, heights
        wrfile.write(
            "{}\n{}\n{}\n{}".format(
                "Number of Labels", "log(O/E) Agreement",
                str(nl), str(ys)
            )
        )

    ys =[reports[k]["general Cohens Kappa score"] for k in nl]
    plt.bar(list(nl), list(ys), color="grey")
    plt.title("Post-clustering Progress")
    plt.xlabel("Number of Labels")
    plt.ylabel("Cohen's Kappa")
    plt.savefig('{}/post_clustering/ck_Progress.pdf'.format(pltsavedir), format='pdf')
    plt.savefig('{}/post_clustering/ck_Progress.svg'.format(pltsavedir), format='svg')
    plt.clf()

    with open('{}/post_clustering/ck_Progress.txt'.format(pltsavedir), "w") as wrfile:
        # xlabel, ylabel, bars, heights
        wrfile.write(
            "{}\n{}\n{}\n{}".format(
                "Number of Labels", "Cohen's Kappa",
                str(nl), str(ys)
            )
        )

    sns.reset_orig
    plt.style.use('default')

def run_single_reprod_analysis(input_dict):
    print("running type: {}".format(input_dict["runtype"]))

    with open(input_dict["output_dir"]+"/run_info.txt", 'w') as fw:
        fw.write(str(input_dict))
    try:
        full_reproducibility_report(input_dict["rep1_dir"], input_dict["rep2_dir"], input_dict["output_dir"], run_on_subset=False)
    except:
        pass
    
def run_single_short_report(input_dict):
    print("running type: {}".format(input_dict["runtype"]))

    with open(input_dict["output_dir"]+"/run_info.txt", 'w') as fw:
        fw.write(str(input_dict))

    if "segway" in input_dict["rep1_dir"]:
        type = "segw"
    elif "chromhmm" in input_dict["rep1_dir"]:
        type = "chmm"

    try:
        get_short_report(input_dict["rep1_dir"], input_dict["rep2_dir"], input_dict["output_dir"], type=type)
    except:
        pass

def RUN_ALL_REPROD_ANALYSIS(runs_dir, CellType_list, output_dir, multi_p=True, type="segway", n_processors=8, run="long"):
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

        if os.path.exists("{}/{}".format(output_dir, ct))==False:
            os.mkdir("{}/{}".format(output_dir, ct))

        ################################################################################
        if os.path.exists("{}/{}/rep1_vs_rep2/".format(output_dir, ct))==False:
            os.mkdir("{}/{}/rep1_vs_rep2/".format(output_dir, ct))
        
        ct_runs["replicates"] = [
            "{}/{}_rep1/parsed_posterior.csv".format(runs_dir, ct), 
            "{}/{}_rep2/parsed_posterior.csv".format(runs_dir, ct), 
            "{}/{}/rep1_vs_rep2/".format(output_dir, ct)]
        
        ################################################################################
        if os.path.exists("{}/{}/rep1_paraminit".format(output_dir, ct))==False:
            os.mkdir("{}/{}/rep1_paraminit".format(output_dir, ct))
        
        if type=="segway":
            ct_runs["rep1_paraminit"] = [
                "{}/{}_rep1/parsed_posterior.csv".format(runs_dir, ct), 
                "{}/{}_rep1_rs5/parsed_posterior.csv".format(runs_dir, ct), 
                "{}/{}/rep1_paraminit/".format(output_dir, ct)]

        elif type == "chmm":
            ct_runs["rep1_paraminit"] = [
                "{}/{}_rep1_{}/parsed_posterior.csv".format(runs_dir, ct, random_seeds[0]), 
                "{}/{}_rep1_{}/parsed_posterior.csv".format(runs_dir, ct, random_seeds[1]), 
                "{}/{}/rep1_paraminit/".format(output_dir, ct)]

        ################################################################################
        if os.path.exists("{}/{}/rep1_pseudoreps".format(output_dir, ct))==False:
            os.mkdir("{}/{}/rep1_pseudoreps".format(output_dir, ct))
        
        if type=="chmm":
            ct_runs["rep1_pseudoreps"] = [
            "{}/{}_rep1psd_concat/parsed_posterior_rep1_psd1.csv".format(runs_dir, ct), 
            "{}/{}_rep1psd_concat/parsed_posterior_rep1_psd2.csv".format(runs_dir, ct), 
            "{}/{}/rep1_pseudoreps/".format(output_dir, ct)]

        elif type == "segway":
            ct_runs["rep1_pseudoreps"] = [
            "{}/{}_concat_rep1psd1/parsed_posterior.csv".format(runs_dir, ct), 
            "{}/{}_concat_rep1psd2/parsed_posterior.csv".format(runs_dir, ct), 
            "{}/{}/rep1_pseudoreps/".format(output_dir, ct)]
        
        ################################################################################
        if os.path.exists("{}/{}/concatenated/".format(output_dir, ct))==False:
            os.mkdir("{}/{}/concatenated/".format(output_dir, ct))

        if type == "segway":
            ct_runs["concat"] = [
            "{}/{}_concat_rep1/parsed_posterior.csv".format(runs_dir, ct), 
            "{}/{}_concat_rep2/parsed_posterior.csv".format(runs_dir, ct), 
            "{}/{}/concatenated/".format(output_dir, ct)]

        elif type == "chmm":
            ct_runs["concat"] = [
            "{}/{}_concat/parsed_posterior_rep1.csv".format(runs_dir, ct), 
            "{}/{}_concat/parsed_posterior_rep2.csv".format(runs_dir, ct), 
            "{}/{}/concatenated/".format(output_dir, ct)]

        ################################################################################
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
    if run=="long":
        if multi_p:
            p_obj = mp.Pool(n_processors)
            p_obj.map(run_single_reprod_analysis, list_of_runs)

        else:
            for r in list_of_runs:
                run_single_reprod_analysis(r)

    elif run == "short":
        if multi_p:
            p_obj = mp.Pool(n_processors)
            p_obj.map(run_single_short_report, list_of_runs)

        else:
            for r in list_of_runs:
                run_single_short_report(r)

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




