import os, sys
import pandas as pd

def intersect_bed(original_bed, include_bed, output_bed):
    cmdline = """bedtools intersect -a {} -b {} > {}"""
    os.system(
        cmdline.format(
           original_bed, include_bed, output_bed 
        )
    )

def create_cellmarkfiletable(celltype_dir, out_filename):
    '''tab delimited file each row contains the 
    cell type the all the associated marks'''
    walk_obj = os.walk(celltype_dir)
    dirs_and_subdirs = [x for x in walk_obj]
    with open(out_filename, 'w') as cmft_file:
        for i in range(1, len(dirs_and_subdirs)):
            for j in dirs_and_subdirs[i][2]:
                if "bedGraph" in j:
                    cmft_file.write('{}\t{}\t{}\n'.format(
                        celltype_dir, dirs_and_subdirs[i][0].split('/')[1],dirs_and_subdirs[i][0]+'/'+j
                    ))

def binarize_data(inputbeddir, cellmarkfiletable, outputdir, resolution=100, chromlength='hg19'):
    cmdline = "java -mx1600M -jar ChromHMM.jar BinarizeBed -center -b {} {} {} {} {}".format(
        resolution, chromlength, inputbeddir, cellmarkfiletable, outputdir
    )
    os.system(cmdline)

def prepare_inputdir():
    pass

def run_chromhmm(binary_input_dir, output_dir, num_labels='16', assembly='hg19', n_threads='0'):
    learnmodel_cmdline = "java -mx1600M -jar ChromHMM.jar LearnModel -init random -printposterior -p {} {} {} {} {}".format(
        n_threads, binary_input_dir, output_dir, num_labels, assembly
    )
    os.system(learnmodel_cmdline)

def QC():
    pass

def gather_results():
    pass

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
