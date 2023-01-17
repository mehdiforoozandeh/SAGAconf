from reports import *

def r1vsr2():
    ################### Rep1 vs Rep2 ###################
        
        ######## GM12878 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/GM12878_rep1/", 
        replicate_2_dir="chromhmm_runs/GM12878_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/GM12878_rep1/", 
        replicate_2_dir="segway_runs/GM12878_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")
            
        ######## MCF-7 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/MCF-7_rep1/", 
        replicate_2_dir="chromhmm_runs/MCF-7_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/MCF-7_rep1/", 
        replicate_2_dir="segway_runs/MCF-7_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

        ######## CD14 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/CD14-positive_monocyte_rep1/", 
        replicate_2_dir="chromhmm_runs/CD14-positive_monocyte_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/CD14-positive_monocyte_rep1/", 
        replicate_2_dir="segway_runs/CD14-positive_monocyte_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")
        
        ######## K562 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/K562_rep1/", 
        replicate_2_dir="chromhmm_runs/K562_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/K562_rep1/", 
        replicate_2_dir="segway_runs/K562_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

        ######## HeLa-S3 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/HeLa-S3_rep1/", 
        replicate_2_dir="chromhmm_runs/HeLa-S3_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/HeLa-S3_rep1/", 
        replicate_2_dir="segway_runs/HeLa-S3_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

def concat():
    ################### concatenated ###################

        ######## GM12878 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/GM12878_concat/", 
        replicate_2_dir="chromhmm_runs/GM12878_concat/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/GM12878_concat_rep1/", 
        replicate_2_dir="segway_runs/GM12878_concat_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")
            
        ######## MCF-7 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/MCF-7_concat/", 
        replicate_2_dir="chromhmm_runs/MCF-7_concat/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/MCF-7_concat_rep1/", 
        replicate_2_dir="segway_runs/MCF-7_concat_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

        ######## CD14 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/CD14-positive_monocyte_concat/", 
        replicate_2_dir="chromhmm_runs/CD14-positive_monocyte_concat/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/CD14-positive_monocyte_concat_rep1/", 
        replicate_2_dir="segway_runs/CD14-positive_monocyte_concat_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")
        
        ######## K562 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/K562_concat/", 
        replicate_2_dir="chromhmm_runs/K562_concat/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/K562_concat_rep1/", 
        replicate_2_dir="segway_runs/K562_concat_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

        ######## HeLa-S3 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/HeLa-S3_concat/", 
        replicate_2_dir="chromhmm_runs/HeLa-S3_concat/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/HeLa-S3_concat_rep1/", 
        replicate_2_dir="segway_runs/HeLa-S3_concat_rep2/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

def paraminit():
    ################### param - init ###################

        ######## GM12878 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/GM12878_rep1_rs5/", 
        replicate_2_dir="chromhmm_runs/GM12878_rep1_rs27/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/GM12878_rep1_rs5/", 
        replicate_2_dir="segway_runs/GM12878_rep1_rs7/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")
            
        ######## MCF-7 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/MCF-7_rep1_rs5/", 
        replicate_2_dir="chromhmm_runs/MCF-7_rep1_rs27/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/MCF-7_rep1_rs5/", 
        replicate_2_dir="segway_runs/MCF-7_rep1_rs7/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

        ######## CD14 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/CD14-positive_monocyte_rep1_rs5/", 
        replicate_2_dir="chromhmm_runs/CD14-positive_monocyte_rep1_rs27/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/CD14-positive_monocyte_rep1_rs5/", 
        replicate_2_dir="segway_runs/CD14-positive_monocyte_rep1_rs7/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")
        
        ######## K562 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/K562_rep1_rs5/", 
        replicate_2_dir="chromhmm_runs/K562_rep1_rs27/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/K562_rep1_rs5/", 
        replicate_2_dir="segway_runs/K562_rep1_rs7/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

        ######## HeLa-S3 ########
    GET_ALL(
        replicate_1_dir="chromhmm_runs/HeLa-S3_rep1_rs5/", 
        replicate_2_dir="chromhmm_runs/HeLa-S3_rep1_rs27/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")

    GET_ALL(
        replicate_1_dir="segway_runs/HeLa-S3_rep1_rs5/", 
        replicate_2_dir="segway_runs/HeLa-S3_rep1_rs7/", 
        genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
        rnaseq="")