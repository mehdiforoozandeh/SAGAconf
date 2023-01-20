from reports import *
import multiprocessing as mp

def r1vsr2(maindir="runs012023_subset"):
    ################### Rep1 vs Rep2 ###################
        ######## GM12878 ########
    
    if os.path.exists(maindir)==False:
        os.mkdir(maindir)

    if os.path.exists(maindir+"/r1vsr2")==False:
        os.mkdir(maindir+"/r1vsr2")

    if os.path.exists(maindir+"/r1vsr2/chmm")==False:
        os.mkdir(maindir+"/r1vsr2/chmm")

    if os.path.exists(maindir+"/r1vsr2/segway")==False:
        os.mkdir(maindir+"/r1vsr2/segway")

    listofruns = [
        {"replicate_1_dir":"chromhmm_runs/GM12878_rep1/", 
        "replicate_2_dir":"chromhmm_runs/GM12878_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv",
        "savedir":"{}/r1vsr2/chmm/gm12878/".format(maindir)},

    
        {"replicate_1_dir":"segway_runs/GM12878_rep1/", 
        "replicate_2_dir":"segway_runs/GM12878_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv",
        "savedir":"{}/r1vsr2/segway/gm12878/".format(maindir)},
            
        ######## MCF-7 ########
    
        {"replicate_1_dir":"chromhmm_runs/MCF-7_rep1/", 
        "replicate_2_dir":"chromhmm_runs/MCF-7_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl",
        "savedir":"{}/r1vsr2/chmm/MCF-7/".format(maindir)},

    
        {"replicate_1_dir":"segway_runs/MCF-7_rep1/", 
        "replicate_2_dir":"segway_runs/MCF-7_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl",
        "savedir":"{}/r1vsr2/segway/MCF-7/".format(maindir)},

        ######## CD14 ########
    
        {"replicate_1_dir":"chromhmm_runs/CD14-positive_monocyte_rep1/", 
        "replicate_2_dir":"chromhmm_runs/CD14-positive_monocyte_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode,_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/r1vsr2/chmm/CD14/".format(maindir)},

    
        {"replicate_1_dir":"segway_runs/CD14-positive_monocyte_rep1/", 
        "replicate_2_dir":"segway_runs/CD14-positive_monocyte_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/r1vsr2/segway/CD14/".format(maindir)},
        
        ######## K562 ########
    
        {"replicate_1_dir":"chromhmm_runs/K562_rep1/", 
        "replicate_2_dir":"chromhmm_runs/K562_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv",
        "savedir":"{}/r1vsr2/chmm/K562/".format(maindir)},

    
        {"replicate_1_dir":"segway_runs/K562_rep1/", 
        "replicate_2_dir":"segway_runs/K562_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv",
        "savedir":"{}/r1vsr2/segway/K562/".format(maindir)},

        ######## HeLa-S3 ########
    
        {"replicate_1_dir":"chromhmm_runs/HeLa-S3_rep1/", 
        "replicate_2_dir":"chromhmm_runs/HeLa-S3_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/r1vsr2/chmm/HeLa-S3/".format(maindir)},

        {"replicate_1_dir":"segway_runs/HeLa-S3_rep1/", 
        "replicate_2_dir":"segway_runs/HeLa-S3_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/r1vsr2/segway/HeLa-S3/".format(maindir)}
        ]
    return listofruns

def concat(maindir="runs012023_subset"):
    if os.path.exists(maindir)==False:
        os.mkdir(maindir)
    
    if os.path.exists(maindir+"/concat")==False:
        os.mkdir(maindir+"/concat")
    
    if os.path.exists(maindir+"/concat/chmm")==False:
        os.mkdir(maindir+"/concat/chmm")

    if os.path.exists(maindir+"/concat/segway")==False:
        os.mkdir(maindir+"/concat/segway")

    ################### concatenated ###################
    listofruns = [
        ######## GM12878 ########
        {"replicate_1_dir":"chromhmm_runs/GM12878_concat/", 
        "replicate_2_dir":"chromhmm_runs/GM12878_concat/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv",
        "savedir":"{}/concat/chmm/GM12878/".format(maindir)},

        {"replicate_1_dir":"segway_runs/GM12878_concat_rep1/", 
        "replicate_2_dir":"segway_runs/GM12878_concat_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv",
        "savedir":"{}/concat/segway/GM12878/".format(maindir)},
            
        ######## MCF-7 ########

        {"replicate_1_dir":"chromhmm_runs/MCF-7_concat/", 
        "replicate_2_dir":"chromhmm_runs/MCF-7_concat/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl",
        "savedir":"{}/concat/chmm/MCF-7/".format(maindir)},

        {"replicate_1_dir":"segway_runs/MCF-7_concat_rep1/", 
        "replicate_2_dir":"segway_runs/MCF-7_concat_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl",
        "savedir":"{}/concat/segway/MCF-7/".format(maindir)},

        ######## CD14 ########

        {"replicate_1_dir":"chromhmm_runs/CD14-positive_monocyte_concat/", 
        "replicate_2_dir":"chromhmm_runs/CD14-positive_monocyte_concat/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/concat/chmm/CD14/".format(maindir)},

        {"replicate_1_dir":"segway_runs/CD14-positive_monocyte_concat_rep1/", 
        "replicate_2_dir":"segway_runs/CD14-positive_monocyte_concat_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/concat/segway/CD14/".format(maindir)},
        
        ######## K562 ########

        {"replicate_1_dir":"chromhmm_runs/K562_concat/", 
        "replicate_2_dir":"chromhmm_runs/K562_concat/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv",
        "savedir":"{}/concat/chmm/K562/".format(maindir)},

        {"replicate_1_dir":"segway_runs/K562_concat_rep1/", 
        "replicate_2_dir":"segway_runs/K562_concat_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv",
        "savedir":"{}/concat/segway/K562/".format(maindir)},

        ######## HeLa-S3 ########

        {"replicate_1_dir":"chromhmm_runs/HeLa-S3_concat/", 
        "replicate_2_dir":"chromhmm_runs/HeLa-S3_concat/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/concat/chmm/HeLa-S3/".format(maindir)},

        {"replicate_1_dir":"segway_runs/HeLa-S3_concat_rep1/", 
        "replicate_2_dir":"segway_runs/HeLa-S3_concat_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/concat/segway/HeLa-S3/".format(maindir)}
    ]
    return listofruns

def paraminit(maindir="runs012023_subset"):
    if os.path.exists(maindir)==False:
        os.mkdir(maindir)
    
    if os.path.exists(maindir+"/paraminit")==False:
        os.mkdir(maindir+"/paraminit")
    
    if os.path.exists(maindir+"/paraminit/chmm")==False:
        os.mkdir(maindir+"/paraminit/chmm")

    if os.path.exists(maindir+"/paraminit/segway")==False:
        os.mkdir(maindir+"/paraminit/segway")

    ################### param - init ###################
    listofruns = [
        ######## GM12878 ########

        {"replicate_1_dir":"chromhmm_runs/GM12878_rep1_rs5/", 
        "replicate_2_dir":"chromhmm_runs/GM12878_rep1_rs27/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv",
        "savedir":"{}/paraminit/chmm/GM12878/".format(maindir)},


        {"replicate_1_dir":"segway_runs/GM12878_rep1_rs5/", 
        "replicate_2_dir":"segway_runs/GM12878_rep1_rs7/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv",
        "savedir":"{}/paraminit/segway/GM12878/".format(maindir)},
            
        ######## MCF-7 ########

        {"replicate_1_dir":"chromhmm_runs/MCF-7_rep1_rs5/", 
        "replicate_2_dir":"chromhmm_runs/MCF-7_rep1_rs27/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl",
        "savedir":"{}/paraminit/chmm/MCF-7/".format(maindir)},


        {"replicate_1_dir":"segway_runs/MCF-7_rep1_rs5/", 
        "replicate_2_dir":"segway_runs/MCF-7_rep1_rs7/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/MCF-7/geneExp_dict_ENCFF721BRA.pkl",
        "savedir":"{}/paraminit/segway/MCF-7/".format(maindir)},

        ######## CD14 ########

        {"replicate_1_dir":"chromhmm_runs/CD14-positive_monocyte_rep1_rs5/", 
        "replicate_2_dir":"chromhmm_runs/CD14-positive_monocyte_rep1_rs27/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/paraminit/chmm/CD14/".format(maindir)},


        {"replicate_1_dir":"segway_runs/CD14-positive_monocyte_rep1_rs5/", 
        "replicate_2_dir":"segway_runs/CD14-positive_monocyte_rep1_rs7/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/paraminit/segway/CD14/".format(maindir)},
        
        ######## K562 ########

        {"replicate_1_dir":"chromhmm_runs/K562_rep1_rs5/", 
        "replicate_2_dir":"chromhmm_runs/K562_rep1_rs27/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv",
        "savedir":"{}/paraminit/chmm/K562/".format(maindir)},


        {"replicate_1_dir":"segway_runs/K562_rep1_rs5/", 
        "replicate_2_dir":"segway_runs/K562_rep1_rs7/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv",
        "savedir":"{}/paraminit/segway/K562/".format(maindir)},

        ######## HeLa-S3 ########

        {"replicate_1_dir":"chromhmm_runs/HeLa-S3_rep1_rs5/", 
        "replicate_2_dir":"chromhmm_runs/HeLa-S3_rep1_rs27/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/paraminit/chmm/HeLa-S3/".format(maindir)},


        {"replicate_1_dir":"segway_runs/HeLa-S3_rep1_rs5/", 
        "replicate_2_dir":"segway_runs/HeLa-S3_rep1_rs7/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/paraminit/segway/HeLa-S3/".format(maindir)}
    ]
    return listofruns

def run(param_dict):
    GET_ALL(
        replicate_1_dir=param_dict["replicate_1_dir"], 
        replicate_2_dir=param_dict["replicate_2_dir"], 
        genecode_dir=param_dict["genecode_dir"], 
        savedir=param_dict["savedir"], 
        rnaseq=param_dict["rnaseq"], 
        contour=True
    )
    print("\n")

def m_p(nt=10):
    pool = mp.Pool(nt)
    _ = pool.map(run, r1vsr2())
    _ = pool.map(run, concat())
    _ = pool.map(run, paraminit())

if __name__=="__main__":
    m_p()