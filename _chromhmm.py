import os, sys
import pandas as pd

def intersect_bed(original_bed, include_bed, output_bed):
    cmdline = """bedtools intersect -a {} -b {} > {}"""
    os.system(
        cmdline.format(
           original_bed, include_bed, output_bed 
        )
    )

def create_cellmarkfiletable(celltype_dir, out_filename, suffix_to_look_for="bedGraph"):
    '''tab delimited file each row contains the 
    cell type the all the associated marks'''
    walk_obj = os.walk(celltype_dir)
    dirs_and_subdirs = [x for x in walk_obj]
    with open(out_filename, 'w') as cmft_file:
        for i in range(1, len(dirs_and_subdirs)):
            for j in dirs_and_subdirs[i][2]:
                if suffix_to_look_for in j:
                    cmft_file.write('{}\t{}\t{}\n'.format(
                        celltype_dir, dirs_and_subdirs[i][0].split('/')[1],dirs_and_subdirs[i][0]+'/'+j
                    ))

def binarize_data(inputbeddir, cellmarkfiletable, outputdir, resolution=100, chromlength='CHROMSIZES/hg19.txt'):
    cmdline = "java -mx1600M -jar ChromHMM.jar BinarizeBed -center -b {} {} {} {} {}".format(
        resolution, chromlength, inputbeddir, cellmarkfiletable, outputdir
    )
    os.system(cmdline)

def learnModel(binary_input_dir, output_dir, num_labels='16', assembly='hg19', n_threads='0'):
    learnmodel_cmdline = "java -mx1600M -jar ChromHMM.jar LearnModel -init random -printposterior -p {} {} {} {} {}".format(
        n_threads, binary_input_dir, output_dir, num_labels, assembly
    )
    os.system(learnmodel_cmdline)


def run_chromhmm(
    celltype_dir,
    include_bed='encodePilotRegions.hg19.bed', resolution=100, chromlength='CHROMSIZES/hg19.txt', 
    num_labels='16', assembly='hg19', n_threads='0'):
    '''run the whole chromhmm pipeline for a celltype'''

    cellmarkfiletable_filename =  celltype_dir + '/cmft.txt'
    binarized_files_directory = celltype_dir + '/binarized_data'
    final_output_directory = celltype_dir + '/chromhmm_output'

    # if include!= None: for i in original bedgraphs: trim i
    if include_bed != None:
        walk_obj = os.walk(celltype_dir)
        dirs_and_subdirs = [x for x in walk_obj]
        for i in range(1, len(dirs_and_subdirs)):
            for j in dirs_and_subdirs[i][2]:
                if "bedGraph" in j:
                    original_bed = dirs_and_subdirs[i][0]+'/'+j
                    output_bed = original_bed.replace(".bedGraph", "_trimmed.BED")
                    intersect_bed(original_bed, include_bed, output_bed)
                
    # create cmft
    if include_bed != None:
        create_cellmarkfiletable(celltype_dir, cellmarkfiletable_filename, suffix_to_look_for="bedGraph")
    else:
        create_cellmarkfiletable(celltype_dir, cellmarkfiletable_filename, suffix_to_look_for="_trimmed.BED")

    # binarize data
    binarize_data(
        celltype_dir, cellmarkfiletable_filename, binarized_files_directory, 
        resolution=resolution, chromlength=chromlength)

    # run learn model
    learnModel(
        binarized_files_directory, final_output_directory, 
        num_labels=num_labels, assembly=assembly, n_threads=n_threads)


def read_emissions():
    '''
    return DF
    '''
    pass

def read_posteriors():
    '''
    return DF
    '''
    pass