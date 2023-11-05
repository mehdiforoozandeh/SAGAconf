# from ..reports import *
import multiprocessing as mp
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging, ast, sys, shutil, os
import matplotlib.patches as mpatches


def r1vsr2(maindir="subset"):
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
        "savedir":"{}/r1vsr2/chmm/GM12878/".format(maindir)},

    
        {"replicate_1_dir":"segway_runs/GM12878_rep1/", 
        "replicate_2_dir":"segway_runs/GM12878_rep2/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv",
        "savedir":"{}/r1vsr2/segway/GM12878/".format(maindir)},
            
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
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
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

def concat(maindir="subset"):
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

def paraminit(maindir="subset"):
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

        {"replicate_1_dir":"chromhmm_runs/GM12878_rep2_rs5/", 
        "replicate_2_dir":"chromhmm_runs/GM12878_rep2_rs27/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv",
        "savedir":"{}/paraminit/chmm/GM12878/".format(maindir)},


        {"replicate_1_dir":"segway_runs/GM12878_rep1/", 
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

        {"replicate_1_dir":"chromhmm_runs/K562_rep2_rs5/", 
        "replicate_2_dir":"chromhmm_runs/K562_rep2_rs27/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv",
        "savedir":"{}/paraminit/chmm/K562/".format(maindir)},


        {"replicate_1_dir":"segway_runs/K562_rep1/", 
        "replicate_2_dir":"segway_runs/K562_rep1_rs7/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":"biovalidation/RNA_seq/K562/preferred_default_ENCFF840UYD.tsv",
        "savedir":"{}/paraminit/segway/K562/".format(maindir)},

        ######## HeLa-S3 ########

        {"replicate_1_dir":"chromhmm_runs/HeLa-S3_rep2_rs5/", 
        "replicate_2_dir":"chromhmm_runs/HeLa-S3_rep2_rs27/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/paraminit/chmm/HeLa-S3/".format(maindir)},


        {"replicate_1_dir":"segway_runs/HeLa-S3_rep1/", 
        "replicate_2_dir":"segway_runs/HeLa-S3_rep1_rs7/", 
        "genecode_dir":"biovalidation/parsed_genecode_data_hg38_release42.csv", 
        "rnaseq":None,
        "savedir":"{}/paraminit/segway/HeLa-S3/".format(maindir)}
    ]
    return listofruns

def run(param_dict):

    print("RUNNING {} VS {}".format(param_dict["replicate_1_dir"], param_dict["replicate_2_dir"]))
    try:
        GET_ALL(
            replicate_1_dir=param_dict["replicate_1_dir"], 
            replicate_2_dir=param_dict["replicate_2_dir"], 
            genecode_dir=param_dict["genecode_dir"], 
            savedir=param_dict["savedir"], 
            rnaseq=param_dict["rnaseq"], 
            contour=False
        )
        with open(param_dict["savedir"]+"/run_info.txt", "w") as f:
            f.write(str(param_dict))

        print("RUNNING {} VS {} is OVER!".format(param_dict["replicate_1_dir"], param_dict["replicate_2_dir"]))
        print("\n")
        
    except Exception as e:
        print("failed at running {} VS {}".format(param_dict["replicate_1_dir"], param_dict["replicate_2_dir"]))

        with open(param_dict["savedir"]+"/run_info.txt", "w") as f:
            f.write(str(param_dict))
            f.write("\nFAILED!")
            f.write("\n\n" + str(e))
            
    finally:
        pass

def m_p(setting, nt=10):
    if setting == "s1":
        with mp.Pool(nt) as pool:
            r_ = pool.map(run, r1vsr2())

    elif setting == "s2":
        with mp.Pool(nt) as pool:
            c = pool.map(run, concat())

    elif setting == "s3":
        with mp.Pool(nt) as pool:
            p = pool.map(run, paraminit())
    
class COMPARATIVE(object):
    def __init__(self, maindir):
        "navigate celltypes, settings, and saga model directories"
        "filter the ones with faulty output"
        self.navigate_results = {}

        for s in [x for x in os.listdir(maindir) if os.path.isdir("{}/{}".format(maindir, x))]:
            self.navigate_results[s] = {}

            for m in [xx for xx in os.listdir("{}/{}".format(maindir, s)) if os.path.isdir("{}/{}/{}".format(maindir, s, xx))]:
                self.navigate_results[s][m] = {}

                for c in [xxx for xxx in os.listdir("{}/{}/{}".format(maindir, s, m)) if os.path.isdir("{}/{}/{}/{}".format(maindir, s, m, xxx))]:

                    self.navigate_results[s][m][c] = ["{}/{}/{}/{}".format(maindir, s, m, c)]
            
        # print(self.navigate_results)
        self.maindir = maindir
        self.var_setting_dict = {
            "concat":"S2: Diff. data, Same model", "r1vsr2": "S1: Diff. data, Diff. model", 
            "paraminit":"S3: Same data, Diff. model"}
        
        self.coverage_quantizer = {(0.0, 0.02):10, (0.02, 0.2):15, (0.2, 0.5):20, (0.5, 1):25}

        self.var_setting_dict_inverse = {v:k for k, v in self.var_setting_dict.items()}
        self.SORT_ORDER = {"Prom": 0, "Prom_fla":1, "Enha":2, "Enha_low":3, "Biva":4, "Tran":5, "Cons":6, "Facu":7, "K9K3":8, "Quie":9}

    # self.navigate_results[s][m][c] = [
    #         dir, [NMI_map, NMI_post], [ovr_ratio_w0, ovr_ratio_w1000], {auc/mauc}, general_rep, {ratio_robust}, {coverage}
    #         {nmi_rec}, {robust_rec}, {avgr_rec}, {robust_rec_entropy}, {NMI_rec_entropy}, {avgr_entropy}, {MI_entropy}, {MI_rec},
    #         {naive_entropy}, {naive_rec}]

    def compare_NMI(self):
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s]:
                for c in self.navigate_results[s][m]:
                    
                    if os.path.exists(self.navigate_results[s][m][c][0] + "/NMI.txt"):
                        dir_NMI_file = self.navigate_results[s][m][c][0] + "/NMI.txt"

                        with open(dir_NMI_file, "r") as nmif:
                            nmilines = nmif.readlines()
                        
                        nmilines[0] = float(nmilines[0].split("=")[1][:7])
                        nmilines[1] = float(nmilines[1].split("=")[2][:7])

                        self.navigate_results[s][m][c].append(nmilines)

                    else:
                        self.navigate_results[s][m][c].append([None, None])

    def compare_ratio_raw_overlap(self):
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s]:
                for c in self.navigate_results[s][m]:
                    
                    if os.path.exists(self.navigate_results[s][m][c][0] + "/raw_conditional_overlap_ratio.txt"):
                        dir_ovr_file = self.navigate_results[s][m][c][0] + "/raw_conditional_overlap_ratio.txt"

                        with open(dir_ovr_file, "r") as ovrf:
                            ovrlines = ovrf.readlines()
                            ovrlines[0] = float(ovrlines[0].split(":")[1][:7])
                            ovrlines[1] = float(ovrlines[1].split(":")[1][:7])

                        self.navigate_results[s][m][c].append([ovrlines[0], ovrlines[1]])

                    else:
                        self.navigate_results[s][m][c].append([None, None])

    def compare_AUC_mAUC(self):
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s]:
                for c in self.navigate_results[s][m]:
                    
                    if os.path.exists(self.navigate_results[s][m][c][0] + "/AUC_mAUC.txt"):
                        dir_auc_file = self.navigate_results[s][m][c][0] + "/AUC_mAUC.txt"
                        auc_mauc = ast.literal_eval(open(dir_auc_file, "r").read())
                        self.navigate_results[s][m][c].append(auc_mauc)

                    else:
                        self.navigate_results[s][m][c].append({})
    
    def load_coverages(self):
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s]:
                for c in self.navigate_results[s][m]:
                    
                    if os.path.exists(self.navigate_results[s][m][c][0] + "/coverages1.txt"):
                        coverage_file = self.navigate_results[s][m][c][0] + "/coverages1.txt"
                        coverage = ast.literal_eval(open(coverage_file, "r").read())
                        self.navigate_results[s][m][c].append(coverage)

                    else:
                        self.navigate_results[s][m][c].append({})                 

    def load_post_clustering_progress(self):
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s]:
                for c in self.navigate_results[s][m]:
                    progress_file = self.navigate_results[s][m][c][0] + "/post_clustering_progress.txt"
                   
                    if os.path.exists(progress_file):
                        with open(progress_file,"r") as f:
                            lines = f.readlines()

                        nmi_rec = ast.literal_eval(lines[0][:-1])
                        robust_rec = ast.literal_eval(lines[1][:-1])
                        avgr_rec = ast.literal_eval(lines[2][:-1])
                        robust_rec_entropy = ast.literal_eval(lines[3][:-1])
                        NMI_rec_entropy = ast.literal_eval(lines[4][:-1])
                        avgr_entropy = ast.literal_eval(lines[5])

                        self.navigate_results[s][m][c].append(nmi_rec)
                        self.navigate_results[s][m][c].append(robust_rec)
                        self.navigate_results[s][m][c].append(avgr_rec)
                        self.navigate_results[s][m][c].append(robust_rec_entropy)
                        self.navigate_results[s][m][c].append(NMI_rec_entropy)
                        self.navigate_results[s][m][c].append(avgr_entropy)

                    else:
                        self.navigate_results[s][m][c].append({})  
                        self.navigate_results[s][m][c].append({})  
                        self.navigate_results[s][m][c].append({})  
                        self.navigate_results[s][m][c].append({})  
                        self.navigate_results[s][m][c].append({})  
                        self.navigate_results[s][m][c].append({})  

                    progress_file = self.navigate_results[s][m][c][0] + "/MI_post_clustering_progress.txt"
                   
                    if os.path.exists(progress_file):
                        with open(progress_file,"r") as f:
                            lines = f.readlines()

                        MI_rec = ast.literal_eval(lines[0][:-1])
                        MI_entropy = ast.literal_eval(lines[1][:-1])

                        self.navigate_results[s][m][c].append(MI_entropy)
                        self.navigate_results[s][m][c].append(MI_rec)

                    else:
                        self.navigate_results[s][m][c].append({})  
                        self.navigate_results[s][m][c].append({})  

                    progress_file = self.navigate_results[s][m][c][0] + "/naive_post_clustering_progress.txt"
                   
                    if os.path.exists(progress_file):
                        with open(progress_file,"r") as f:
                            lines = f.readlines()

                        naive_rec = ast.literal_eval(lines[0][:-1])
                        naive_entropy = ast.literal_eval(lines[1][:-1])

                        self.navigate_results[s][m][c].append(naive_entropy)
                        self.navigate_results[s][m][c].append(naive_rec)

                    else:
                        self.navigate_results[s][m][c].append({})  
                        self.navigate_results[s][m][c].append({})

    def compare_ratio_robust(self):
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s]:
                for c in self.navigate_results[s][m]:
                    
                    if os.path.exists(self.navigate_results[s][m][c][0] + "/general_reproducibility_score_POSTERIOR.txt"):
                        dir_robust_file = self.navigate_results[s][m][c][0] + "/general_reproducibility_score_POSTERIOR.txt"
                        ratio_robust = {}

                        rlines = open(dir_robust_file, "r").readlines()
                        for l in rlines:
                            ratio_robust[l.split(" = ")[0].replace(" reprod score", "")] = float(l.split(" = ")[1][:-1])

                        if "general" in ratio_robust.keys():
                            self.navigate_results[s][m][c].append(ratio_robust["general"])
                            del ratio_robust["general"]
                        
                        self.navigate_results[s][m][c].append(ratio_robust)

                    else:
                        self.navigate_results[s][m][c].append(None)
                        self.navigate_results[s][m][c].append({})

    def compare_ratio_raw_overlap_perlabel(self):
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s]:
                for c in self.navigate_results[s][m]:
                    
                    if os.path.exists(self.navigate_results[s][m][c][0] + "/raw_conditional_overlap_ratio.txt"):
                        dir_ovr_file = self.navigate_results[s][m][c][0] + "/raw_conditional_overlap_ratio.txt"

                        with open(dir_ovr_file, "r") as ovrf:
                            ovrlines = ovrf.readlines()
                            raw_ovr = {}
                            for i in range(2, len(ovrlines)):
                                raw_ovr[str(ovrlines[i].split(" : ")[0])] = float(ovrlines[i].split(":")[1][:7])

                        self.navigate_results[s][m][c].append(raw_ovr)

                    else:
                        self.navigate_results[s][m][c].append({})

    def visualize_NMI(self):
        chmm = {}
        segway = {}

        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():

                    if m == "segway":
                        if s not in segway.keys():
                            segway[s] = {}

                        segway[s][c] = self.navigate_results[s][m][c][1]

                    elif m == "chmm":
                        if s not in chmm.keys():
                            chmm[s] = {}
                        
                        chmm[s][c] = self.navigate_results[s][m][c][1]

        df_chmm = []
        for s in chmm.keys():
            for c in chmm[s].keys():
                df_chmm.append([c, self.var_setting_dict[s], chmm[s][c][0], chmm[s][c][1]])

        df_segway = []
        for s in segway.keys():
            for c in segway[s].keys():
                df_segway.append([c, self.var_setting_dict[s], segway[s][c][0], segway[s][c][1]])

        df_chmm = pd.DataFrame(
            df_chmm, columns = ["CellType", "Setting", "NMI_MAP", "NMI_post"]).sort_values(by=["Setting", "CellType"])
        df_segway = pd.DataFrame(
            df_segway, columns = ["CellType", "Setting", "NMI_MAP", "NMI_post"]).sort_values(by=["Setting", "CellType"])

        sns.set_theme(style="whitegrid")
        sns.reset_orig
        plt.style.use('default')

        sns.set(rc={'figure.figsize':(10,8)})
        sns.set_palette(sns.color_palette("deep"))

        sns.barplot(x="CellType", y="NMI_MAP", hue="Setting", data=df_chmm)
        plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
        plt.ylabel("NMI")
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.tight_layout()

        plt.savefig(self.maindir+"/NMI_Comparison_ChormHMM.pdf", format="pdf")
        plt.savefig(self.maindir+"/NMI_Comparison_ChormHMM.svg", format="svg")
        plt.clf()

        sns.set_theme(style="whitegrid")
        sns.reset_orig
        plt.style.use('default')

        sns.set(rc={'figure.figsize':(10,8)})
        sns.set_palette(sns.color_palette("deep"))

        sns.barplot(x="CellType", y="NMI_MAP", hue="Setting", data=df_segway)
        plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
        plt.ylabel("NMI")
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.tight_layout()
        plt.savefig(self.maindir+"/NMI_Comparison_Segway.pdf", format="pdf")
        plt.savefig(self.maindir+"/NMI_Comparison_Segway.svg", format="svg")
        plt.clf()

        sns.set_theme(style="whitegrid")
        sns.reset_orig
        plt.style.use('default')

        sns.set(rc={'figure.figsize':(10,8)})
        sns.set_palette(sns.color_palette("deep"))

        sns.barplot(x="CellType", y="NMI_post", hue="Setting", data=df_chmm)
        plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
        plt.ylabel("pNMI")
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.tight_layout()

        plt.savefig(self.maindir+"/pNMI_Comparison_ChormHMM.pdf", format="pdf")
        plt.savefig(self.maindir+"/pNMI_Comparison_ChormHMM.svg", format="svg")
        plt.clf()

        sns.set_theme(style="whitegrid")
        sns.reset_orig
        plt.style.use('default')

        sns.set(rc={'figure.figsize':(10,8)})
        sns.set_palette(sns.color_palette("deep"))

        sns.barplot(x="CellType", y="NMI_post", hue="Setting", data=df_segway)
        plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
        plt.ylabel("pNMI")
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.tight_layout()
        plt.savefig(self.maindir+"/pNMI_Comparison_Segway.pdf", format="pdf")
        plt.savefig(self.maindir+"/pNMI_Comparison_Segway.svg", format="svg")
        plt.clf()
    
    def visualize_ovr(self):
        chmm = {}
        segway = {}

        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():

                    if m == "segway":
                        if s not in segway.keys():
                            segway[s] = {}

                        segway[s][c] = self.navigate_results[s][m][c][2]

                    elif m == "chmm":
                        if s not in chmm.keys():
                            chmm[s] = {}
                        
                        chmm[s][c] = self.navigate_results[s][m][c][2]

        df_chmm = []
        for s in chmm.keys():
            for c in chmm[s].keys():
                df_chmm.append([c, self.var_setting_dict[s], chmm[s][c][0], chmm[s][c][1]])

        df_segway = []
        for s in segway.keys():
            for c in segway[s].keys():
                df_segway.append([c, self.var_setting_dict[s], segway[s][c][0], segway[s][c][1]])

        df_chmm = pd.DataFrame(
            df_chmm, columns = ["CellType", "Setting", "ovr_ratio_w0", "ovr_ratio_w1000"]).sort_values(by=["Setting", "CellType"])
        df_segway = pd.DataFrame(
            df_segway, columns = ["CellType", "Setting", "ovr_ratio_w0", "ovr_ratio_w1000"]).sort_values(by=["Setting", "CellType"])

        sns.set_theme(style="whitegrid")
        sns.reset_orig
        plt.style.use('default')

        sns.set(rc={'figure.figsize':(10,8)})
        sns.set_palette(sns.color_palette("deep"))

        sns.barplot(x="CellType", y="ovr_ratio_w0", hue="Setting", data=df_chmm)
        plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
        plt.ylabel("Overlap Ratio | w=0")
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.tight_layout()

        plt.savefig(self.maindir+"/ovrw0_Comparison_ChormHMM.pdf", format="pdf")
        plt.savefig(self.maindir+"/ovrw0_Comparison_ChormHMM.svg", format="svg")
        plt.clf()

        ######################################################################################################
        sns.set_theme(style="whitegrid")
        sns.reset_orig
        plt.style.use('default')

        sns.set(rc={'figure.figsize':(10,8)})
        sns.set_palette(sns.color_palette("deep"))

        sns.barplot(x="CellType", y="ovr_ratio_w0", hue="Setting", data=df_segway)
        plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
        plt.ylabel("Overlap Ratio | w=0")
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.tight_layout()
        plt.savefig(self.maindir+"/ovrw0_Comparison_Segway.pdf", format="pdf")
        plt.savefig(self.maindir+"/ovrw0_Comparison_Segway.svg", format="svg")
        plt.clf()
        ######################################################################################################
        sns.set_theme(style="whitegrid")
        sns.reset_orig
        plt.style.use('default')

        sns.set(rc={'figure.figsize':(10,8)})
        sns.set_palette(sns.color_palette("deep"))

        sns.barplot(x="CellType", y="ovr_ratio_w1000", hue="Setting", data=df_chmm)
        plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
        plt.ylabel("Overlap Ratio | w=1000bp")
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.tight_layout()

        plt.savefig(self.maindir+"/ovrw1000_Comparison_ChormHMM.pdf", format="pdf")
        plt.savefig(self.maindir+"/ovrw1000_Comparison_ChormHMM.svg", format="svg")
        plt.clf()
        ######################################################################################################
        sns.set_theme(style="whitegrid")
        sns.reset_orig
        plt.style.use('default')

        sns.set(rc={'figure.figsize':(10,8)})
        sns.set_palette(sns.color_palette("deep"))

        sns.barplot(x="CellType", y="ovr_ratio_w1000", hue="Setting", data=df_segway)
        plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
        plt.ylabel("Overlap Ratio | w=1000bp")
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.tight_layout()
        plt.savefig(self.maindir+"/ovrw1000_Comparison_Segway.pdf", format="pdf")
        plt.savefig(self.maindir+"/ovrw1000_Comparison_Segway.svg", format="svg")
        plt.clf()

    def visualize_AUC_mAUC(self):
        chmm = {}
        segway = {}

        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():

                    if m == "segway":
                        if s not in segway.keys():
                            segway[s] = {}

                        smc_dict = self.navigate_results[s][m][c][3]
                        smc_coverage = self.navigate_results[s][m][c][6]


                        smc_list = []
                        for k in smc_dict.keys():
                            new_k = "_".join(k.split("_")[1:])

                            for q in self.coverage_quantizer.keys():
                                if q[0] <= smc_coverage[k] and smc_coverage[k] < q[1]:
                                    smc_list.append([new_k, smc_dict[k], self.coverage_quantizer[q]])

                        segway[s][c] = smc_list

                    elif m == "chmm":
                        if s not in chmm.keys():
                            chmm[s] = {}
                        
                        smc_dict = self.navigate_results[s][m][c][3]
                        smc_coverage = self.navigate_results[s][m][c][6]

                        smc_list = []
                        for k in smc_dict.keys():
                            new_k = "_".join(k.split("_")[1:])
                            for q in self.coverage_quantizer.keys():
                                if q[0] <= smc_coverage[k] and smc_coverage[k] < q[1]:
                                    smc_list.append([new_k, smc_dict[k], self.coverage_quantizer[q]])

                        chmm[s][c] = smc_list

        df_chmm = []
        for s in chmm.keys():
            for c in chmm[s].keys():
                for l in chmm[s][c]:
                    if l[1] > 0:
                        df_chmm.append([c, self.var_setting_dict[s], l[0], l[1], l[2]])
        
        df_chmm = pd.DataFrame(
            df_chmm, columns = ["CellType", "Setting", "Genomic_function", "AUC/maxAUC", "coverage"]).sort_values(
                by=["Setting", "CellType"])

        df_segway = []
        for s in segway.keys():
            for c in segway[s].keys():
                for l in segway[s][c]:
                    if l[1] > 0:
                        df_segway.append([c, self.var_setting_dict[s], l[0], l[1], l[2]])
        
        df_segway = pd.DataFrame(
            df_segway, columns = ["CellType", "Setting", "Genomic_function", "AUC/maxAUC", "coverage"]).sort_values(
                by=["Setting", "CellType"])

        ####################################################################################################################################
        for s in np.unique(df_chmm["Setting"]):
            plt.clf()
            sns.reset_orig
            plt.style.use('default')

            fig, ax = plt.subplots(figsize=(10, 8))

            chmmdata = df_chmm.loc[df_chmm["Setting"]==s,:]

            order = chmmdata['Genomic_function'].unique()
            unique_coverages = chmmdata['coverage'].unique()
            unique_celltypes = chmmdata['CellType'].unique()
            palette = sns.color_palette("deep", len(unique_celltypes))

            for coverage in unique_coverages:
                mask = chmmdata['coverage'] == coverage
                size = coverage
                for i, celltype in enumerate(unique_celltypes):
                    celltype_mask = chmmdata['CellType'] == celltype
                    data = chmmdata[mask & celltype_mask]
                    if not data.empty:
                        sns.stripplot(
                            x='Genomic_function', y="AUC/maxAUC", data=data, ax=ax, size=size, 
                            color=palette[i], order=order, alpha=0.85)

            n_categories = len(ax.get_xticks())

            for i in range(n_categories):
                ax.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

            # ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)
            handles = [plt.Line2D([], [], marker='o', color=palette[i], linestyle='None') for i in range(len(unique_celltypes))]
            ax.legend(
                handles, unique_celltypes, bbox_to_anchor=(0,1.02,1,0.2), 
                loc="lower left", mode="expand", borderaxespad=0, ncol=5, fontsize=12)

            # Show the plot
            plt.yticks(np.arange(0.5, 1.05, step=0.1), fontsize=12)
            plt.xticks(fontsize=12)
            plt.tight_layout()

            plt.savefig("{}/{}/{}/AUC_mAUC_stripplot.pdf".format(self.maindir, self.var_setting_dict_inverse[s], "chmm"), format="pdf")
            plt.savefig("{}/{}/{}/AUC_mAUC_stripplot.svg".format(self.maindir, self.var_setting_dict_inverse[s], "chmm"), format="svg")

            plt.clf()
            sns.reset_orig
            plt.style.use('default')

            
        ####################################################################################################################################
        for s in np.unique(df_segway["Setting"]):
            fig, ax = plt.subplots(figsize=(10, 8))

            segdata = df_segway.loc[df_segway["Setting"]==s,:]

            order = segdata['Genomic_function'].unique()
            unique_coverages = segdata['coverage'].unique()
            unique_celltypes = segdata['CellType'].unique()
            palette = sns.color_palette("deep", len(unique_celltypes))

            for coverage in unique_coverages:
                mask = segdata['coverage'] == coverage
                size = coverage
                for i, celltype in enumerate(unique_celltypes):
                    celltype_mask = segdata['CellType'] == celltype
                    data = segdata[mask & celltype_mask]
                    if not data.empty:
                        sns.stripplot(
                            x='Genomic_function', y="AUC/maxAUC", data=data, ax=ax, size=size, 
                            color=palette[i], order=order, alpha=0.85)

            n_categories = len(ax.get_xticks())

            for i in range(n_categories):
                ax.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

            # ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)
            handles = [plt.Line2D([], [], marker='o', color=palette[i], linestyle='None') for i in range(len(unique_celltypes))]
            ax.legend(handles, unique_celltypes, bbox_to_anchor=(0,1.02,1,0.2), 
            loc="lower left", mode="expand", borderaxespad=0, ncol=5, fontsize=12)

            # Show the plot
            plt.yticks(np.arange(0.5, 1.05, step=0.1), fontsize=12)
            plt.xticks(fontsize=12)
            plt.tight_layout()

            plt.savefig("{}/{}/{}/AUC_mAUC_stripplot.pdf".format(self.maindir, self.var_setting_dict_inverse[s], "segway"), format="pdf")
            plt.savefig("{}/{}/{}/AUC_mAUC_stripplot.svg".format(self.maindir, self.var_setting_dict_inverse[s], "segway"), format="svg")

            plt.clf()
            sns.reset_orig
            plt.style.use('default')
        ####################################################################################################################################

    def visualize_robust(self):
        chmm = {}
        segway = {}

        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():

                    if m == "segway":
                        if s not in segway.keys():
                            segway[s] = {}

                        smc_dict = self.navigate_results[s][m][c][5]
                        smc_coverage = self.navigate_results[s][m][c][6]


                        smc_list = []
                        for k in smc_dict.keys():
                            new_k = "_".join(k.split("_")[1:])

                            for q in self.coverage_quantizer.keys():
                                if q[0] <= smc_coverage[k] and smc_coverage[k] <= q[1]:
                                    smc_list.append([new_k, smc_dict[k], self.coverage_quantizer[q]])

                        segway[s][c] = smc_list

                    elif m == "chmm":
                        if s not in chmm.keys():
                            chmm[s] = {}
                        
                        smc_dict = self.navigate_results[s][m][c][5]
                        smc_coverage = self.navigate_results[s][m][c][6]

                        smc_list = []
                        for k in smc_dict.keys():
                            new_k = "_".join(k.split("_")[1:])

                            for q in self.coverage_quantizer.keys():
                                if q[0] <= smc_coverage[k] and smc_coverage[k] < q[1]:
                                    smc_list.append([new_k, smc_dict[k], self.coverage_quantizer[q]])
                                    # if smc_dict[k] == 0:
                                        # print(s, m, c, k, smc_dict[k], smc_coverage[k])
                                
                        chmm[s][c] = smc_list

        df_chmm = []
        for s in chmm.keys():
            for c in chmm[s].keys():
                for l in chmm[s][c]:
                    df_chmm.append([c, self.var_setting_dict[s], l[0], l[1], l[2]])
        
        df_chmm = pd.DataFrame(
            df_chmm, columns = ["CellType", "Setting", "Genomic_function", "ratio_robust", "coverage"]).sort_values(
                by=["Setting", "CellType"])
        
        df_chmm['sort_col'] = df_chmm["Genomic_function"].map(self.SORT_ORDER)
        df_chmm.sort_values("sort_col", inplace=True)
        df_chmm.drop("sort_col", axis=1, inplace=True)

        df_segway = []
        for s in segway.keys():
            for c in segway[s].keys():
                for l in segway[s][c]:
                    df_segway.append([c, self.var_setting_dict[s], l[0], l[1], l[2]])
        
        df_segway = pd.DataFrame(
            df_segway, columns = ["CellType", "Setting", "Genomic_function", "ratio_robust", "coverage"]).sort_values(
                by=["Setting", "CellType"])
        
        df_segway['sort_col'] = df_segway["Genomic_function"].map(self.SORT_ORDER)
        df_segway.sort_values("sort_col", inplace=True)
        df_segway.drop("sort_col", axis=1, inplace=True)
        ####################################################################################################################################
        for s in np.unique(df_chmm["Setting"]):
            fig, ax = plt.subplots(figsize=(10, 8))

            chmmdata = df_chmm.loc[df_chmm["Setting"]==s,:]
            

            order = chmmdata['Genomic_function'].unique()
            unique_coverages = chmmdata['coverage'].unique()
            unique_celltypes = chmmdata['CellType'].unique()
            palette = sns.color_palette("deep", len(unique_celltypes))

                
            for coverage in unique_coverages:
                mask = chmmdata['coverage'] == coverage
                size = coverage
                for i, celltype in enumerate(unique_celltypes):
                    celltype_mask = chmmdata['CellType'] == celltype
                    data = chmmdata[mask & celltype_mask]
                    if not data.empty:
                        sns.stripplot(
                            x='Genomic_function', y="ratio_robust", data=data, ax=ax, size=size, 
                            color=palette[i], order=order, alpha=0.85)

            n_categories = len(ax.get_xticks())

            for i in range(n_categories):
                ax.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

            # ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)
            handles = [plt.Line2D([], [], marker='o', color=palette[i], linestyle='None') for i in range(len(unique_celltypes))]
            ax.legend(handles, unique_celltypes, bbox_to_anchor=(0,1.02,1,0.2), 
            loc="lower left", mode="expand", borderaxespad=0, ncol=5, fontsize=12)
            
            # Show the plot
            plt.yticks(np.arange(0, 1.05, step=0.1), fontsize=12)
            plt.xticks(fontsize=12)
            plt.tight_layout()
            plt.savefig("{}/{}/{}/ratio_robust_stripplot.pdf".format(self.maindir, self.var_setting_dict_inverse[s], "chmm"), format="pdf")
            plt.savefig("{}/{}/{}/ratio_robust_stripplot.svg".format(self.maindir, self.var_setting_dict_inverse[s], "chmm"), format="svg")

            plt.clf()
            sns.reset_orig
            plt.style.use('default')
        ####################################################################################################################################
        for s in np.unique(df_segway["Setting"]):
            fig, ax = plt.subplots(figsize=(10, 8))
            
            segdata = df_segway.loc[df_segway["Setting"]==s,:]

            order = segdata['Genomic_function'].unique()
            unique_coverages = segdata['coverage'].unique()
            unique_celltypes = segdata['CellType'].unique()
            palette = sns.color_palette("deep", len(unique_celltypes))

            for coverage in unique_coverages:
                mask = segdata['coverage'] == coverage
                size = coverage
                for i, celltype in enumerate(unique_celltypes):
                    celltype_mask = segdata['CellType'] == celltype
                    data = segdata[mask & celltype_mask]
                    if not data.empty:
                        sns.stripplot(
                            x='Genomic_function', y="ratio_robust", data=data, ax=ax, size=size, 
                            color=palette[i], order=order, alpha=0.85)

            n_categories = len(ax.get_xticks())

            for i in range(n_categories):
                ax.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

            # ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)
            handles = [plt.Line2D([], [], marker='o', color=palette[i], linestyle='None') for i in range(len(unique_celltypes))]
            ax.legend(handles, unique_celltypes, bbox_to_anchor=(0,1.02,1,0.2), 
            loc="lower left", mode="expand", borderaxespad=0, ncol=5, fontsize=12)

            # Show the plot
            plt.yticks(np.arange(0, 1.05, step=0.1), fontsize=12)
            plt.xticks(fontsize=12)
            plt.tight_layout()
            plt.savefig("{}/{}/{}/ratio_robust_stripplot.pdf".format(self.maindir, self.var_setting_dict_inverse[s], "segway"), format="pdf")
            plt.savefig("{}/{}/{}/ratio_robust_stripplot.svg".format(self.maindir, self.var_setting_dict_inverse[s], "segway"), format="svg")

            plt.clf()
            sns.reset_orig
            plt.style.use('default')

    def visualize_raw_ovr(self):
        chmm = {}
        segway = {}

        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():

                    if m == "segway":
                        if s not in segway.keys():
                            segway[s] = {}

                        smc_dict = self.navigate_results[s][m][c][15]
                        smc_coverage = self.navigate_results[s][m][c][6]


                        smc_list = []
                        for k in smc_dict.keys():
                            new_k = "_".join(k.split("_")[1:])

                            for q in self.coverage_quantizer.keys():
                                if q[0] <= smc_coverage[k] and smc_coverage[k] < q[1]:
                                    smc_list.append([new_k, smc_dict[k], self.coverage_quantizer[q]])

                        segway[s][c] = smc_list

                    elif m == "chmm":
                        if s not in chmm.keys():
                            chmm[s] = {}
                        
                        smc_dict = self.navigate_results[s][m][c][15]
                        smc_coverage = self.navigate_results[s][m][c][6]

                        smc_list = []
                        for k in smc_dict.keys():
                            new_k = "_".join(k.split("_")[1:])
                            for q in self.coverage_quantizer.keys():
                                if q[0] <= smc_coverage[k] and smc_coverage[k] < q[1]:
                                    smc_list.append([new_k, smc_dict[k], self.coverage_quantizer[q]])

                        chmm[s][c] = smc_list

        df_chmm = []
        for s in chmm.keys():
            for c in chmm[s].keys():
                for l in chmm[s][c]:
                    if l[1] > 0:
                        df_chmm.append([c, self.var_setting_dict[s], l[0], l[1], l[2]])
        
        df_chmm = pd.DataFrame(
            df_chmm, columns = ["CellType", "Setting", "Genomic_function", "raw_ovr", "coverage"]).sort_values(
                by=["Setting", "CellType"])

        df_segway = []
        for s in segway.keys():
            for c in segway[s].keys():
                for l in segway[s][c]:
                    if l[1] > 0:
                        df_segway.append([c, self.var_setting_dict[s], l[0], l[1], l[2]])
        
        df_segway = pd.DataFrame(
            df_segway, columns = ["CellType", "Setting", "Genomic_function", "raw_ovr", "coverage"]).sort_values(
                by=["Setting", "CellType"])

        ####################################################################################################################################
        for s in np.unique(df_chmm["Setting"]):
            plt.clf()
            sns.reset_orig
            plt.style.use('default')

            fig, ax = plt.subplots(figsize=(10, 8))

            chmmdata = df_chmm.loc[df_chmm["Setting"]==s,:]

            order = chmmdata['Genomic_function'].unique()
            unique_coverages = chmmdata['coverage'].unique()
            unique_celltypes = chmmdata['CellType'].unique()
            palette = sns.color_palette("deep", len(unique_celltypes))

            for coverage in unique_coverages:
                mask = chmmdata['coverage'] == coverage
                size = coverage
                for i, celltype in enumerate(unique_celltypes):
                    celltype_mask = chmmdata['CellType'] == celltype
                    data = chmmdata[mask & celltype_mask]
                    if not data.empty:
                        sns.stripplot(
                            x='Genomic_function', y="raw_ovr", data=data, ax=ax, size=size, 
                            color=palette[i], order=order, alpha=0.85)

            n_categories = len(ax.get_xticks())

            for i in range(n_categories):
                ax.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

            # ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)
            handles = [plt.Line2D([], [], marker='o', color=palette[i], linestyle='None') for i in range(len(unique_celltypes))]
            ax.legend(
                handles, unique_celltypes, bbox_to_anchor=(0,1.02,1,0.2), 
                loc="lower left", mode="expand", borderaxespad=0, ncol=5, fontsize=12)

            # Show the plot
            plt.yticks(np.arange(0, 1.05, step=0.2), fontsize=12)
            plt.xticks(fontsize=12)
            plt.tight_layout()

            plt.savefig("{}/{}/{}/raw_ovr_stripplot.pdf".format(self.maindir, self.var_setting_dict_inverse[s], "chmm"), format="pdf")
            plt.savefig("{}/{}/{}/raw_ovr_stripplot.svg".format(self.maindir, self.var_setting_dict_inverse[s], "chmm"), format="svg")

            plt.clf()
            sns.reset_orig
            plt.style.use('default')

            
        ####################################################################################################################################
        for s in np.unique(df_segway["Setting"]):
            fig, ax = plt.subplots(figsize=(10, 8))

            segdata = df_segway.loc[df_segway["Setting"]==s,:]

            order = segdata['Genomic_function'].unique()
            unique_coverages = segdata['coverage'].unique()
            unique_celltypes = segdata['CellType'].unique()
            palette = sns.color_palette("deep", len(unique_celltypes))

            for coverage in unique_coverages:
                mask = segdata['coverage'] == coverage
                size = coverage
                for i, celltype in enumerate(unique_celltypes):
                    celltype_mask = segdata['CellType'] == celltype
                    data = segdata[mask & celltype_mask]
                    if not data.empty:
                        sns.stripplot(
                            x='Genomic_function', y="raw_ovr", data=data, ax=ax, size=size, 
                            color=palette[i], order=order, alpha=0.85)

            n_categories = len(ax.get_xticks())

            for i in range(n_categories):
                ax.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

            # ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)
            handles = [plt.Line2D([], [], marker='o', color=palette[i], linestyle='None') for i in range(len(unique_celltypes))]
            ax.legend(handles, unique_celltypes, bbox_to_anchor=(0,1.02,1,0.2), 
            loc="lower left", mode="expand", borderaxespad=0, ncol=5, fontsize=12)

            # Show the plot
            plt.yticks(np.arange(0, 1.05, step=0.2), fontsize=12)
            plt.xticks(fontsize=12)
            plt.tight_layout()

            plt.savefig("{}/{}/{}/raw_ovr_stripplot.pdf".format(self.maindir, self.var_setting_dict_inverse[s], "segway"), format="pdf")
            plt.savefig("{}/{}/{}/raw_ovr_stripplot.svg".format(self.maindir, self.var_setting_dict_inverse[s], "segway"), format="svg")

            plt.clf()
            sns.reset_orig
            plt.style.use('default')
        ####################################################################################################################################

    def post_clustering_comparatives(self, drop_decrease=True):
        saga_translate = {"chmm":"ChromHMM","segway":"Segway"}
        numlabels_DF = []
        entropy_DF = []

        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():
                    
                    for k in self.navigate_results[s][m][c][7].keys():
                        smc_nlabels_list = [self.var_setting_dict[s], saga_translate[m], c]
                        smc_nlabels_list.append(k)
                        smc_nlabels_list.append(self.navigate_results[s][m][c][7][k])
                        smc_nlabels_list.append(self.navigate_results[s][m][c][8][k])
                        smc_nlabels_list.append(self.navigate_results[s][m][c][9][k])
                        # smc_nlabels_list.append(self.navigate_results[s][m][c][14][k])

                        numlabels_DF.append(smc_nlabels_list)
                    
                    for k in self.navigate_results[s][m][c][10].keys():
                        smc_entropy_list = [self.var_setting_dict[s], saga_translate[m], c]
                        smc_entropy_list.append(k)
                        smc_entropy_list.append(self.navigate_results[s][m][c][10][k])
                        smc_entropy_list.append(self.navigate_results[s][m][c][11][k])
                        smc_entropy_list.append(self.navigate_results[s][m][c][12][k])
                        # smc_entropy_list.append(self.navigate_results[s][m][c][13][k])

                        entropy_DF.append(smc_entropy_list)

        
        numlabels_DF = pd.DataFrame(numlabels_DF, columns= [
            "setting", "saga", "celltype", "num_labels", "nmi", "ratio_robust", "avgr"]).sort_values(by=["setting", "celltype"]) 
        numlabels_DF['num_labels'] = numlabels_DF['num_labels'].astype(int)

        entropy_DF = pd.DataFrame(entropy_DF, columns= [
            "setting", "saga", "celltype", "entropy", "ratio_robust", "nmi", "avgr"]).sort_values(by=["setting", "celltype"]) 
    

        markers = {entropy_DF.setting.unique()[0]: 'o', entropy_DF.setting.unique()[1]: '^', entropy_DF.setting.unique()[2]: 's'}
        linestyles = {
            entropy_DF.setting.unique()[0]: 'dashed', 
            entropy_DF.setting.unique()[1]: "solid", 
            entropy_DF.setting.unique()[2]: 'dotted'}

        colors = sns.color_palette('colorblind', n_colors=len(entropy_DF.celltype.unique()))
        colors = dict(zip(entropy_DF.celltype.unique(), colors))

        methods = {"Segway":1, "ChromHMM":0}
        settings = {
            entropy_DF.setting.unique()[0]: 0, 
            entropy_DF.setting.unique()[1]: 1, 
            entropy_DF.setting.unique()[2]: 2}
        ########################################################################################################
        def dd(x, y):
            xs = []
            ys = []
            max_sofar = y[0]
            for i in range(len(x)):
                if y[i] >= max_sofar:
                    xs.append(x[i])
                    ys.append(y[i])
                    max_sofar = y[i]

            return np.array(xs), np.array(ys)
        ########################################################################################################

        fig, axs = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True)

        for s in entropy_DF.setting.unique():
            for m in entropy_DF.saga.unique():
                for c in entropy_DF.celltype.unique():
                    
                    data = entropy_DF.loc[
                        (entropy_DF["setting"] == s) &  
                        (entropy_DF["saga"] == m) &
                        (entropy_DF["celltype"] == c), :
                    ]
                    
                    data = data.sort_values(by=["entropy"], ascending=False) 

                    xs = np.array(data.entropy)
                    ys = np.array(data.ratio_robust)

                    if drop_decrease:
                        xs, ys = dd(xs, ys)
                    
                    try:
                        axs[methods[m], settings[s]].scatter(xs, ys, marker=markers[s], color=colors[c], label=c)
                        axs[methods[m], settings[s]].plot(xs, ys, linestyle=linestyles[s], color=colors[c])
                        axs[methods[m], settings[s]].set_title(f"{m} | {s}")
                        axs[methods[m], settings[s]].set_xlabel("Entorpy")
                        axs[methods[m], settings[s]].set_ylabel("Ratio Confident")
                        axs[methods[m], settings[s]].legend()
                    except:
                        pass
        
        plt.tight_layout()
        plt.savefig(self.maindir+"/ent_robust.pdf", format="pdf")
        plt.savefig(self.maindir+"/ent_robust.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')

        ########################################################################################################
        fig, axs = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True)

        for s in entropy_DF.setting.unique():
            for m in entropy_DF.saga.unique():
                for c in entropy_DF.celltype.unique():
                    
                    data = entropy_DF.loc[
                        (entropy_DF["setting"] == s) &  
                        (entropy_DF["saga"] == m) &
                        (entropy_DF["celltype"] == c), :
                    ]
                    
                    data = data.sort_values(by=["entropy"], ascending=False) 

                    xs = np.array(data.entropy)
                    ys = np.array(data.avgr)

                    if drop_decrease:
                        xs, ys = dd(xs, ys)
                    
                    try:
                        axs[methods[m], settings[s]].scatter(xs, ys, marker=markers[s], color=colors[c], label=c)
                        axs[methods[m], settings[s]].plot(xs, ys, linestyle=linestyles[s], color=colors[c])
                        axs[methods[m], settings[s]].set_title(f"{m} | {s}")
                        axs[methods[m], settings[s]].set_xlabel("Entorpy")
                        axs[methods[m], settings[s]].set_ylabel("Average r_value")
                        axs[methods[m], settings[s]].legend()
                    except:
                        pass
        
        plt.tight_layout()
        plt.savefig(self.maindir+"/ent_avgr.pdf", format="pdf")
        plt.savefig(self.maindir+"/ent_avgr.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')

        ########################################################################################################
        fig, axs = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True)

        for s in numlabels_DF.setting.unique():
            for m in numlabels_DF.saga.unique():
                for c in numlabels_DF.celltype.unique():
                    
                    data = numlabels_DF.loc[
                        (numlabels_DF["setting"] == s) &  
                        (numlabels_DF["saga"] == m) &
                        (numlabels_DF["celltype"] == c), :
                    ]
                    
                    data = data.sort_values(by=["num_labels"], ascending=False) 

                    xs = np.array(data.num_labels)
                    ys = np.array(data.ratio_robust)

                    if drop_decrease:
                        xs, ys = dd(xs, ys)
                    
                    try:
                        axs[methods[m], settings[s]].scatter(xs, ys, marker=markers[s], color=colors[c], label=c)
                        axs[methods[m], settings[s]].plot(xs, ys, linestyle=linestyles[s], color=colors[c])
                        axs[methods[m], settings[s]].set_title(f"{m} | {s}")
                        axs[methods[m], settings[s]].set_xlabel("Number of Chromatin States")
                        axs[methods[m], settings[s]].set_ylabel("Ratio Confident")
                        axs[methods[m], settings[s]].legend()
                    except:
                        pass
        
        plt.tight_layout()
        plt.savefig(self.maindir+"/nlabels_robust.pdf", format="pdf")
        plt.savefig(self.maindir+"/nlabels_robust.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')

        #######################################################################################################
        fig, axs = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True)

        for s in numlabels_DF.setting.unique():
            for m in numlabels_DF.saga.unique():
                for c in numlabels_DF.celltype.unique():
                    
                    data = numlabels_DF.loc[
                        (numlabels_DF["setting"] == s) &  
                        (numlabels_DF["saga"] == m) &
                        (numlabels_DF["celltype"] == c), :
                    ]
                    
                    data = data.sort_values(by=["num_labels"], ascending=False) 

                    xs = np.array(data.num_labels)
                    ys = np.array(data.avgr)

                    if drop_decrease:
                        xs, ys = dd(xs, ys)
                    
                    try:
                        axs[methods[m], settings[s]].scatter(xs, ys, marker=markers[s], color=colors[c], label=c)
                        axs[methods[m], settings[s]].plot(xs, ys, linestyle=linestyles[s], color=colors[c])
                        axs[methods[m], settings[s]].set_title(f"{m} | {s}")
                        axs[methods[m], settings[s]].set_xlabel("Number of Chromatin States")
                        axs[methods[m], settings[s]].set_ylabel("Average r_value")
                        axs[methods[m], settings[s]].legend()
                    except:
                        pass
        
        plt.tight_layout()
        plt.savefig(self.maindir+"/nlabels_avgr.pdf", format="pdf")
        plt.savefig(self.maindir+"/nlabels_avgr.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')      

    def post_clustering_compressed(self, drop_decrease=True):
        saga_translate = {"chmm":"ChromHMM","segway":"Segway"}
        numlabels_DF = []
        entropy_DF = []

        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():
                    
                    for k in self.navigate_results[s][m][c][7].keys():
                        smc_nlabels_list = [self.var_setting_dict[s], saga_translate[m], c]
                        smc_nlabels_list.append(k)
                        smc_nlabels_list.append(self.navigate_results[s][m][c][7][k])
                        smc_nlabels_list.append(self.navigate_results[s][m][c][8][k])
                        smc_nlabels_list.append(self.navigate_results[s][m][c][9][k])
                        # smc_nlabels_list.append(self.navigate_results[s][m][c][14][k])

                        numlabels_DF.append(smc_nlabels_list)
                    
                    for k in self.navigate_results[s][m][c][10].keys():
                        smc_entropy_list = [self.var_setting_dict[s], saga_translate[m], c]
                        smc_entropy_list.append(k)
                        smc_entropy_list.append(self.navigate_results[s][m][c][10][k])
                        smc_entropy_list.append(self.navigate_results[s][m][c][11][k])
                        smc_entropy_list.append(self.navigate_results[s][m][c][12][k])
                        # smc_entropy_list.append(self.navigate_results[s][m][c][13][k])

                        entropy_DF.append(smc_entropy_list)

        
        numlabels_DF = pd.DataFrame(numlabels_DF, columns= [
            "setting", "saga", "celltype", "num_labels", "nmi", "ratio_robust", "avgr"]).sort_values(by=["setting", "celltype"]) 
        numlabels_DF['num_labels'] = numlabels_DF['num_labels'].astype(int)

        entropy_DF = pd.DataFrame(entropy_DF, columns= [
            "setting", "saga", "celltype", "entropy", "ratio_robust", "nmi", "avgr"]).sort_values(by=["setting", "celltype"]) 
    

        markers = {entropy_DF.setting.unique()[0]: 'o', entropy_DF.setting.unique()[1]: '^', entropy_DF.setting.unique()[2]: 's'}
        linestyles = {
            entropy_DF.saga.unique()[0]: 'solid', 
            entropy_DF.saga.unique()[1]: 'dotted'}

        colors = sns.color_palette('colorblind', n_colors=len(entropy_DF.saga.unique()))
        colors = dict(zip(entropy_DF.saga.unique(), colors))

        methods = {"Segway":1, "ChromHMM":0}
        settings = {
            entropy_DF.setting.unique()[0]: 0, 
            entropy_DF.setting.unique()[1]: 1, 
            entropy_DF.setting.unique()[2]: 2}
        ########################################################################################################
        def dd(x, y):
            xs = []
            ys = []
            max_sofar = y[0]
            for i in range(len(x)):
                if y[i] >= max_sofar:
                    xs.append(x[i])
                    ys.append(y[i])
                    max_sofar = y[i]
                    
            return np.array(xs), np.array(ys)
        ########################################################################################################
        fig, axs = plt.subplots(3, 1, figsize=(4, 8), sharex=True, sharey=True)

        for s in entropy_DF.setting.unique():
            for m in entropy_DF.saga.unique():
                for c in entropy_DF.celltype.unique():
                    
                    data = entropy_DF.loc[
                        (entropy_DF["celltype"] == c) &
                        (entropy_DF["setting"] == s) &  
                        (entropy_DF["saga"] == m), :
                    ]
                    
                    data = data.sort_values(by=["entropy"], ascending=False) 

                    xs = np.array(data.entropy)
                    ys = np.array(data.avgr)
                    if drop_decrease:
                        xs, ys = dd(xs, ys)
                    
                    axs[settings[s]].plot(xs, ys, linestyle=linestyles[m], color=colors[m], label=m)
                    axs[settings[s]].set_title(f"{s}")
                    axs[settings[s]].set_xlabel("Entorpy")
                    axs[settings[s]].set_ylabel("Average r_value")

        lines, labels = axs[settings[s]].get_legend_handles_labels()

        # Create a dictionary to store unique entries
        unique_entries = {}
        # Iterate over the lines and labels
        for line, label in zip(lines, labels):
            # Only add unique entries to the dictionary
            if label not in unique_entries:
                unique_entries[label] = line

        # Create a custom legend using the unique entries
        plt.legend(unique_entries.values(), unique_entries.keys())
        plt.tight_layout()
        plt.savefig(self.maindir+"/compressed_ent_avgr.pdf", format="pdf")
        plt.savefig(self.maindir+"/compressed_ent_avgr.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')
        ########################################################################################################
        fig, axs = plt.subplots(3, 1, figsize=(4, 8), sharex=True, sharey=True)

        for s in entropy_DF.setting.unique():
            for m in entropy_DF.saga.unique():
                for c in entropy_DF.celltype.unique():
                    
                    data = entropy_DF.loc[
                        (entropy_DF["celltype"] == c) &
                        (entropy_DF["setting"] == s) &  
                        (entropy_DF["saga"] == m), :
                    ]
                    
                    data = data.sort_values(by=["entropy"], ascending=False) 

                    xs = np.array(data.entropy)
                    ys = np.array(data.ratio_robust)
                    if drop_decrease:
                        xs, ys = dd(xs, ys)
                    
                    axs[settings[s]].plot(xs, ys, linestyle=linestyles[m], color=colors[m], label=m)
                    axs[settings[s]].set_title(f"{s}")
                    axs[settings[s]].set_xlabel("Entorpy")
                    axs[settings[s]].set_ylabel("Ratio Confident")

        lines, labels = axs[settings[s]].get_legend_handles_labels()

        # Create a dictionary to store unique entries
        unique_entries = {}
        # Iterate over the lines and labels
        for line, label in zip(lines, labels):
            # Only add unique entries to the dictionary
            if label not in unique_entries:
                unique_entries[label] = line

        # Create a custom legend using the unique entries
        plt.legend(unique_entries.values(), unique_entries.keys())
        plt.tight_layout()
        plt.savefig(self.maindir+"/compressed_ent_robust.pdf", format="pdf")
        plt.savefig(self.maindir+"/compressed_ent_robust.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')
        ########################################################################################################
        fig, axs = plt.subplots(3, 1, figsize=(4, 8), sharex=True, sharey=True)

        for s in numlabels_DF.setting.unique():
            for m in numlabels_DF.saga.unique():
                for c in numlabels_DF.celltype.unique():
                    
                    data = numlabels_DF.loc[
                        (numlabels_DF["celltype"] == c) &
                        (numlabels_DF["setting"] == s) &  
                        (numlabels_DF["saga"] == m), :
                    ]
                    
                    data = data.sort_values(by=["num_labels"], ascending=False) 

                    xs = np.array(data.num_labels)
                    ys = np.array(data.avgr)
                    if drop_decrease:
                        xs, ys = dd(xs, ys)
                    
                    axs[settings[s]].plot(xs, ys, linestyle=linestyles[m], color=colors[m], label=m)
                    axs[settings[s]].set_title(f"{s}")
                    axs[settings[s]].set_xlabel("Number of Chromatin states")
                    axs[settings[s]].set_ylabel("Average r_value")

        lines, labels = axs[settings[s]].get_legend_handles_labels()

        # Create a dictionary to store unique entries
        unique_entries = {}
        # Iterate over the lines and labels
        for line, label in zip(lines, labels):
            # Only add unique entries to the dictionary
            if label not in unique_entries:
                unique_entries[label] = line

        # Create a custom legend using the unique entries
        plt.legend(unique_entries.values(), unique_entries.keys())
        plt.tight_layout()
        plt.savefig(self.maindir+"/compressed_nlabels_avgr.pdf", format="pdf")
        plt.savefig(self.maindir+"/compressed_nlabels_avgr.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')
        ########################################################################################################
        fig, axs = plt.subplots(3, 1, figsize=(4, 8), sharex=True, sharey=True)

        for s in numlabels_DF.setting.unique():
            for m in numlabels_DF.saga.unique():
                for c in numlabels_DF.celltype.unique():
                    
                    data = numlabels_DF.loc[
                        (numlabels_DF["celltype"] == c) &
                        (numlabels_DF["setting"] == s) &  
                        (numlabels_DF["saga"] == m), :
                    ]
                    
                    data = data.sort_values(by=["num_labels"], ascending=False) 

                    xs = np.array(data.num_labels)
                    ys = np.array(data.ratio_robust)
                    if drop_decrease:
                        xs, ys = dd(xs, ys)
                    
                    axs[settings[s]].plot(xs, ys, linestyle=linestyles[m], color=colors[m], label=m)
                    axs[settings[s]].set_title(f"{s}")
                    axs[settings[s]].set_xlabel("Number of Chromatin states")
                    axs[settings[s]].set_ylabel("Ratio Confident")

        lines, labels = axs[settings[s]].get_legend_handles_labels()

        # Create a dictionary to store unique entries
        unique_entries = {}
        # Iterate over the lines and labels
        for line, label in zip(lines, labels):
            # Only add unique entries to the dictionary
            if label not in unique_entries:
                unique_entries[label] = line

        # Create a custom legend using the unique entries
        plt.legend(unique_entries.values(), unique_entries.keys())
        plt.tight_layout()
        plt.savefig(self.maindir+"/compressed_nlabels_robust.pdf", format="pdf")
        plt.savefig(self.maindir+"/compressed_nlabels_robust.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')
        
    def post_clustering_MI(self):
        saga_translate = {"chmm":"ChromHMM","segway":"Segway"}
        numlabels_DF = []
        entropy_DF = []

        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():
                    
                    for k in self.navigate_results[s][m][c][14].keys():
                        smc_nlabels_list = [self.var_setting_dict[s], saga_translate[m], c]
                        smc_nlabels_list.append(k)
                        smc_nlabels_list.append(self.navigate_results[s][m][c][14][k])

                        numlabels_DF.append(smc_nlabels_list)
                    
                    for k in self.navigate_results[s][m][c][13].keys():
                        smc_entropy_list = [self.var_setting_dict[s], saga_translate[m], c]
                        smc_entropy_list.append(k)
                        smc_entropy_list.append(self.navigate_results[s][m][c][13][k])

                        entropy_DF.append(smc_entropy_list)

        numlabels_DF = pd.DataFrame(numlabels_DF, columns= [
            "setting", "saga", "celltype", "num_labels", "MI"]).sort_values(by=["setting", "celltype"]) 
        numlabels_DF['num_labels'] = numlabels_DF['num_labels'].astype(int)

        entropy_DF = pd.DataFrame(entropy_DF, columns= [
            "setting", "saga", "celltype", "entropy", "MI"]).sort_values(by=["setting", "celltype"]) 
    

        markers = {entropy_DF.setting.unique()[0]: 'o', entropy_DF.setting.unique()[1]: '^', entropy_DF.setting.unique()[2]: 's'}
        linestyles = {
            entropy_DF.setting.unique()[0]: 'dashed', 
            entropy_DF.setting.unique()[1]: "solid", 
            entropy_DF.setting.unique()[2]: 'dotted'}

        colors = sns.color_palette('colorblind', n_colors=len(entropy_DF.celltype.unique()))
        colors = dict(zip(entropy_DF.celltype.unique(), colors))

        methods = {"Segway":1, "ChromHMM":0}
        settings = {
            entropy_DF.setting.unique()[0]: 0, 
            entropy_DF.setting.unique()[1]: 1, 
            entropy_DF.setting.unique()[2]: 2}

        ########################################################################################################

        fig, axs = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True)

        for s in entropy_DF.setting.unique():
            for m in entropy_DF.saga.unique():
                
                axs[methods[m], settings[s]].plot(
                    [entropy_DF.entropy.min(), entropy_DF.entropy.max()], [
                        entropy_DF.entropy.min(), entropy_DF.entropy.max()], 
                        linestyle="dashed", color="grey")

                for c in entropy_DF.celltype.unique():
                    
                    data = entropy_DF.loc[
                        (entropy_DF["setting"] == s) &  
                        (entropy_DF["saga"] == m) &
                        (entropy_DF["celltype"] == c), :
                    ]
                    
                    data = data.sort_values(by=["entropy"], ascending=False) 

                    xs = np.array(data.entropy)
                    ys = np.array(data.MI)

                    axs[methods[m], settings[s]].scatter(xs, ys, marker=markers[s], color=colors[c], label=c)
                    axs[methods[m], settings[s]].plot(xs, ys, linestyle=linestyles[s], color=colors[c])
                    axs[methods[m], settings[s]].set_title(f"{m} | {s}")
                    axs[methods[m], settings[s]].set_xlabel("Entropy")
                    axs[methods[m], settings[s]].set_ylabel("Mutual Information")
                    axs[methods[m], settings[s]].legend()
    
        plt.tight_layout()
        plt.savefig(self.maindir+"/ent_MI.pdf", format="pdf")
        plt.savefig(self.maindir+"/ent_MI.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')
        ########################################################################################################

        fig, axs = plt.subplots(2, 3, figsize=(15, 10), sharex=True, sharey=True)

        for s in numlabels_DF.setting.unique():
            for m in numlabels_DF.saga.unique():
                for c in numlabels_DF.celltype.unique():
                    
                    data = numlabels_DF.loc[
                        (numlabels_DF["setting"] == s) &  
                        (numlabels_DF["saga"] == m) &
                        (numlabels_DF["celltype"] == c), :
                    ]
                    
                    data = data.sort_values(by=["num_labels"], ascending=False) 

                    xs = np.array(data.num_labels)
                    ys = np.array(data.MI)
                    
                    axs[methods[m], settings[s]].scatter(xs, ys, marker=markers[s], color=colors[c], label=c)
                    axs[methods[m], settings[s]].plot(xs, ys, linestyle=linestyles[s], color=colors[c])
                    axs[methods[m], settings[s]].set_title(f"{m} | {s}")
                    axs[methods[m], settings[s]].set_xlabel("Entropy")
                    axs[methods[m], settings[s]].set_ylabel("Mutual Information")
                    axs[methods[m], settings[s]].legend()
        
        plt.tight_layout()
        plt.savefig(self.maindir+"/num_labels_MI.pdf", format="pdf")
        plt.savefig(self.maindir+"/num_labels_MI.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')

        ########################################################################################################
        linestyles = {
            entropy_DF.saga.unique()[0]: 'solid', 
            entropy_DF.saga.unique()[1]: 'solid'}

        colors = sns.color_palette('colorblind', n_colors=len(entropy_DF.saga.unique()))
        colors = dict(zip(entropy_DF.saga.unique(), colors))

        methods = {"Segway":1, "ChromHMM":0}
        settings = {
            entropy_DF.setting.unique()[0]: 0, 
            entropy_DF.setting.unique()[1]: 1, 
            entropy_DF.setting.unique()[2]: 2}
        ########################################################################################################
        fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)

        for s in entropy_DF.setting.unique():
            for m in entropy_DF.saga.unique():

                axs[settings[s]].plot(
                    [entropy_DF.entropy.min(), entropy_DF.entropy.max()], [
                        entropy_DF.entropy.min(), entropy_DF.entropy.max()], 
                        linestyle="dashed", color="grey")

                for c in entropy_DF.celltype.unique():
                    
                    data = entropy_DF.loc[
                        (entropy_DF["celltype"] == c) &
                        (entropy_DF["setting"] == s) &  
                        (entropy_DF["saga"] == m), :
                    ]
                    
                    data = data.sort_values(by=["entropy"], ascending=False) 

                    xs = np.array(data.entropy)
                    ys = np.array(data.MI)
                    
                    axs[settings[s]].plot(xs, ys, linestyle=linestyles[m], color=colors[m], label=m)
                    axs[settings[s]].set_title(f"{s}")
                    axs[settings[s]].set_xlabel("Entropy")
                    axs[settings[s]].set_ylabel("Mutual Information")

        lines, labels = axs[settings[s]].get_legend_handles_labels()

        # Create a dictionary to store unique entries
        unique_entries = {}
        # Iterate over the lines and labels
        for line, label in zip(lines, labels):
            # Only add unique entries to the dictionary
            if label not in unique_entries:
                unique_entries[label] = line

        # Create a custom legend using the unique entries
        plt.legend(unique_entries.values(), unique_entries.keys())
        plt.tight_layout()
        plt.savefig(self.maindir+"/compressed_ent_MI.pdf", format="pdf")
        plt.savefig(self.maindir+"/compressed_ent_MI.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')
        ########################################################################################################
        fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)

        for s in numlabels_DF.setting.unique():
            for m in numlabels_DF.saga.unique():

                for c in numlabels_DF.celltype.unique():
                    
                    data = numlabels_DF.loc[
                        (numlabels_DF["celltype"] == c) &
                        (numlabels_DF["setting"] == s) &  
                        (numlabels_DF["saga"] == m), :
                    ]
                    
                    data = data.sort_values(by=["num_labels"], ascending=False) 

                    xs = np.array(data.num_labels)
                    ys = np.array(data.MI)
                    
                    axs[settings[s]].plot(xs, ys, linestyle=linestyles[m], color=colors[m], label=m)
                    axs[settings[s]].set_title(f"{s}")
                    axs[settings[s]].set_xlabel("Number of Chromatin states")
                    axs[settings[s]].set_ylabel("Mutual Information")

        lines, labels = axs[settings[s]].get_legend_handles_labels()

        # Create a dictionary to store unique entries
        unique_entries = {}
        # Iterate over the lines and labels
        for line, label in zip(lines, labels):
            # Only add unique entries to the dictionary
            if label not in unique_entries:
                unique_entries[label] = line

        # Create a custom legend using the unique entries
        plt.legend(unique_entries.values(), unique_entries.keys())
        plt.tight_layout()
        plt.savefig(self.maindir+"/compressed_nlabels_MI.pdf", format="pdf")
        plt.savefig(self.maindir+"/compressed_nlabels_MI.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')
        ########################################################################################################
        fig, axs = plt.subplots(1, 1, figsize=(5, 5), sharex=True, sharey=True)

        data = numlabels_DF.loc[
            (numlabels_DF["celltype"] == "GM12878") &
            (numlabels_DF["setting"] == numlabels_DF.setting.unique()[0]) &  
            (numlabels_DF["saga"] == "ChromHMM"), :
        ]
        s= numlabels_DF.setting.unique()[0]
        data = data.sort_values(by=["num_labels"], ascending=False) 

        xs = np.array(data.num_labels)
        ys = np.array(data.MI)
        
        axs.scatter(xs, ys, color=colors["ChromHMM"])
        axs.plot(xs, ys, linestyle=linestyles["ChromHMM"], color=colors["ChromHMM"], label="ChromHMM")
        axs.set_title(f"{s} | GM12878")
        axs.set_xlabel("Number of Chromatin states")
        axs.set_ylabel("Mutual Information")

        lines, labels = axs.get_legend_handles_labels()

        # Create a dictionary to store unique entries
        unique_entries = {}
        # Iterate over the lines and labels
        for line, label in zip(lines, labels):
            # Only add unique entries to the dictionary
            if label not in unique_entries:
                unique_entries[label] = line

        # Create a custom legend using the unique entries
        plt.legend(unique_entries.values(), unique_entries.keys())
        plt.tight_layout()
        plt.savefig(self.maindir+"/compressed_nlabels_MI_EXEMPLAR.pdf", format="pdf")
        plt.savefig(self.maindir+"/compressed_nlabels_MI_EXEMPLAR.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')

    def MPF1(self):
        saga_translate = {"chmm":"ChromHMM","segway":"Segway"}
        general_DF = []
        
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():
                    general_DF.append(
                        [self.var_setting_dict[s], saga_translate[m], c, 
                        self.navigate_results[s][m][c][1][0], self.navigate_results[s][m][c][1][1], 
                        self.navigate_results[s][m][c][2][0], self.navigate_results[s][m][c][2][1]
                        ]
                    )
        
        general_DF = pd.DataFrame(general_DF, columns= [
            "setting", "saga", "celltype", "NMI", "pNMI", "overlap | w=0", "overlap | w=1000"]).sort_values(by=["setting", "celltype"])

        # Create a 2x2 grid of subplots
        fig, axs = plt.subplots(2, 2, figsize=(15, 10), sharex=False, sharey=False)

        # Flatten the axs array to make it easier to iterate over
        axs = axs.flatten()

        # Define the marker styles you want to use for each category of the 'saga' variable
        markers = {general_DF.saga.unique()[0]: 'o', general_DF.saga.unique()[1]: '^'}

        # Plot the data in each subplot
        for saga, marker in markers.items():
            df = general_DF[general_DF['saga'] == saga]
            sns.stripplot(x="setting", y="NMI", hue="celltype", data=df, ax=axs[0], marker=marker, size=10, alpha=0.8)
            sns.stripplot(x="setting", y="pNMI", hue="celltype", data=df, ax=axs[1], marker=marker, size=10, alpha=0.8)
            sns.stripplot(x="setting", y="overlap | w=0", hue="celltype", data=df, ax=axs[2], marker=marker, size=10, alpha=0.8)
            sns.stripplot(x="setting", y="overlap | w=1000", hue="celltype", data=df, ax=axs[3], marker=marker, size=10, alpha=0.8)

        for i in range(len(general_DF.setting.unique())):
            axs[0].axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)
            axs[1].axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)
            axs[2].axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)
            axs[3].axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

        for ax in axs:
            ax.set_xticklabels(ax.get_xticklabels(), fontsize=9)
            ax.set_xlabel('')

        # Get the color palette used by seaborn
        palette = sns.color_palette()

        # Create a custom legend
        legend_elements = []
        for i, celltype in enumerate(general_DF.celltype.unique()):
            legend_elements.append(plt.Line2D([0], [0], color=palette[i], label=celltype))
            
        for i, saga in enumerate(general_DF.saga.unique()):
            legend_elements.append(plt.Line2D([0], [0], color='black', marker=markers[saga], label=saga, linestyle='None'))

        # Add the custom legend to the first subplot
        axs[0].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=7)

        axs[1].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=7)

        axs[2].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=7)

        axs[3].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=7)

        plt.savefig(self.maindir+"/MPF1.pdf", format="pdf")
        plt.savefig(self.maindir+"/MPF1.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')

    def MPF2(self):
        saga_translate = {"chmm":"ChromHMM","segway":"Segway"}
        general_DF = []

        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():
                    general_DF.append(
                        [self.var_setting_dict[s], saga_translate[m], c, 
                        self.navigate_results[s][m][c][1][0], self.navigate_results[s][m][c][1][1], 
                        self.navigate_results[s][m][c][2][0], self.navigate_results[s][m][c][2][1]
                        ]
                    )

        general_DF = pd.DataFrame(general_DF, columns= [
            "setting", "saga", "celltype", "NMI", "pNMI", "overlap | w=0", "overlap | w=1000"]).sort_values(by=["setting", "celltype"])

        # Create a 2x2 grid of subplots
        fig, axs = plt.subplots(2, 2, figsize=(15, 10), sharex=False, sharey=False)

        # Flatten the axs array to make it easier to iterate over
        axs = axs.flatten()

        # Define the marker styles you want to use for each category of the 'celltype' variable
        markers = {general_DF.celltype.unique()[0]: 'o', general_DF.celltype.unique()[1]: '^', general_DF.celltype.unique()[2]: 's', general_DF.celltype.unique()[3]: 'D', general_DF.celltype.unique()[4]: 'v'}

        # Plot the data in each subplot
        for celltype, marker in markers.items():
            df = general_DF[general_DF['celltype'] == celltype]
            sns.stripplot(x="setting", y="NMI", hue="saga", data=df, ax=axs[0], marker=marker, size=11, alpha=0.8)
            sns.stripplot(x="setting", y="pNMI", hue="saga", data=df, ax=axs[1], marker=marker, size=11, alpha=0.8)
            sns.stripplot(x="setting", y="overlap | w=0", hue="saga", data=df, ax=axs[2], marker=marker, size=11, alpha=0.8)
            sns.stripplot(x="setting", y="overlap | w=1000", hue="saga", data=df, ax=axs[3], marker=marker, size=11, alpha=0.8)

        for i in range(len(general_DF.setting.unique())):
            axs[0].axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)
            axs[1].axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)
            axs[2].axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)
            axs[3].axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

        for ax in axs:
            ax.set_xticklabels(ax.get_xticklabels(), fontsize=9)
            ax.set_xlabel('')
        
        axs[0].set_yticks(np.arange(0.2, 1, 0.1))
        axs[0].tick_params(axis='y', labelsize=12)
        axs[1].set_yticks(np.arange(0.2, 1, 0.1))
        axs[1].tick_params(axis='y', labelsize=12)
        axs[2].set_yticks(np.arange(0.3, 1, 0.1))
        axs[2].tick_params(axis='y', labelsize=12)
        axs[3].set_yticks(np.arange(0.3, 1, 0.1))
        axs[3].tick_params(axis='y', labelsize=12)

        # Get the color palette used by seaborn
        palette = sns.color_palette()

        # Create a custom legend
        legend_elements = []
        for i, saga in enumerate(general_DF.saga.unique()):
            legend_elements.append(plt.Line2D([0], [0], color=palette[i], label=saga))
            
        for i, celltype in enumerate(general_DF.celltype.unique()):
            legend_elements.append(plt.Line2D([0], [0], color='black', marker=markers[celltype], label=celltype, linestyle='None'))
    
         # Add the custom legend to the first subplot
        axs[0].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=9)

        axs[1].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=9)

        axs[2].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=9)

        axs[3].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=9)

        plt.savefig(self.maindir+"/MPF2.pdf", format="pdf")
        plt.savefig(self.maindir+"/MPF2.svg", format="svg")

        plt.clf()
        sns.reset_orig
        plt.style.use('default')
    
    def integ_ratio_robust(self):
        saga_translate = {"chmm":"ChromHMM","segway":"Segway"}
        general_DF = []
        
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():
                    general_DF.append(
                        [self.var_setting_dict[s], saga_translate[m], c, 
                        self.navigate_results[s][m][c][4]])
        
        general_DF = pd.DataFrame(general_DF, columns= [
            "setting", "saga", "celltype", "ratio_confident"]).sort_values(by=["setting", "celltype"])

        fig, axs = plt.subplots(1, 1, figsize=(8, 6), sharex=False, sharey=False)

        # Flatten the axs array to make it easier to iterate over
        # axs = axs.flatten()

        # Define the marker styles you want to use for each category of the 'celltype' variable
        markers = {
            general_DF.celltype.unique()[0]: 'o', general_DF.celltype.unique()[1]: '^', 
            general_DF.celltype.unique()[2]: 's', general_DF.celltype.unique()[3]: 'D', 
            general_DF.celltype.unique()[4]: 'v'}

        # Plot the data in each subplot
        for celltype, marker in markers.items():
            df = general_DF[general_DF['celltype'] == celltype]
            sns.stripplot(x="setting", y="ratio_confident", hue="saga", data=df, ax=axs, marker=marker, size=11, alpha=0.8)

        for i in range(len(general_DF.setting.unique())):
            axs.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

        axs.set_xticklabels(axs.get_xticklabels(), fontsize=9)
        axs.set_xlabel('')

        axs.tick_params(axis='y', labelsize=12)

        # Get the color palette used by seaborn
        palette = sns.color_palette()

        # Create a custom legend
        legend_elements = []
        for i, saga in enumerate(general_DF.saga.unique()):
            legend_elements.append(plt.Line2D([0], [0], color=palette[i], label=saga))
            
        for i, celltype in enumerate(general_DF.celltype.unique()):
            legend_elements.append(plt.Line2D([0], [0], color='black', marker=markers[celltype], label=celltype, linestyle='None'))
    
         # Add the custom legend to the first subplot
        axs.legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=9)
        
        plt.savefig(self.maindir+"/integ_ratiorobust.pdf", format="pdf")
        plt.savefig(self.maindir+"/integ_ratiorobust.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')

    def integ_auc_mauc(self):
        saga_translate = {"chmm":"ChromHMM","segway":"Segway"}
        general_DF = []
        
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s].keys():
                for c in self.navigate_results[s][m].keys():
                    avg_auc = 0
                    for k in self.navigate_results[s][m][c][3].keys():
                        avg_auc += self.navigate_results[s][m][c][3][k] * self.navigate_results[s][m][c][6][k]

                    general_DF.append(
                        [self.var_setting_dict[s], saga_translate[m], c, 
                        avg_auc])
        
        general_DF = pd.DataFrame(general_DF, columns= [
            "setting", "saga", "celltype", "average AUC/mAUC"]).sort_values(by=["setting", "celltype"])

        fig, axs = plt.subplots(1, 1, figsize=(8, 6), sharex=False, sharey=False)

        # Flatten the axs array to make it easier to iterate over
        # axs = axs.flatten()

        # Define the marker styles you want to use for each category of the 'celltype' variable
        markers = {
            general_DF.celltype.unique()[0]: 'o', general_DF.celltype.unique()[1]: '^', 
            general_DF.celltype.unique()[2]: 's', general_DF.celltype.unique()[3]: 'D', 
            general_DF.celltype.unique()[4]: 'v'}

        # Plot the data in each subplot
        for celltype, marker in markers.items():
            df = general_DF[general_DF['celltype'] == celltype]
            sns.stripplot(x="setting", y="average AUC/mAUC", hue="saga", data=df, ax=axs, marker=marker, size=10, alpha=0.8)

        for i in range(len(general_DF.setting.unique())):
            axs.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

        axs.set_xticklabels(axs.get_xticklabels(), fontsize=9)
        axs.set_xlabel('')

        axs.tick_params(axis='y', labelsize=12)

        # Get the color palette used by seaborn
        palette = sns.color_palette()

        # Create a custom legend
        legend_elements = []
        for i, saga in enumerate(general_DF.saga.unique()):
            legend_elements.append(plt.Line2D([0], [0], color=palette[i], label=saga))
            
        for i, celltype in enumerate(general_DF.celltype.unique()):
            legend_elements.append(plt.Line2D([0], [0], color='black', marker=markers[celltype], label=celltype, linestyle='None'))
    
         # Add the custom legend to the first subplot
        axs.legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=9)

        plt.savefig(self.maindir+"/integ_auc_mauc.pdf", format="pdf")
        plt.savefig(self.maindir+"/integ_auc_mauc.svg", format="svg")
        plt.clf()
        sns.reset_orig
        plt.style.use('default')

    def ALL(self):
        try:
            self.compare_NMI()
        except:
            print("failed at compare_NMI")
        try:
            self.compare_ratio_raw_overlap()
        except:
            print("failed at compare_ratio_raw_overlap")
        try:
            self.compare_AUC_mAUC()
        except:
            print("failed at compare_AUC_mAUC")
        try:
            self.compare_ratio_robust()
        except:
            print("failed at compare_ratio_robust")
        try:
            self.load_coverages()
        except:
            print("failed at load_coverages")
        try:
            self.load_post_clustering_progress()
        except:
            print("failed at load_post_clustering_progress")
        try:
            self.compare_ratio_raw_overlap_perlabel()
        except:
            print("failed at compare_ratio_raw_overlap_perlabel")
        # try:
        #     self.post_clustering_compressed(drop_decrease=True)
        # except:
        #     print("failed at post_clustering_compressed")
        # try:
        #     self.post_clustering_comparatives()    
        # except:
        #     print("failed at post_clustering_comparatives")
        # try:
        #     self.integ_ratio_robust()
        # except:
        #     print("failed at integ_ratio_robust")
        # try:
        #     self.integ_auc_mauc()
        # except:
        #     print("failed at integ_auc_mauc")
        # try:
        #     self.visualize_raw_ovr()
        # except:
        #     print("failed at visualize_raw_ovr")
        # try:
        #     self.MPF2()
        # except:
        #     print("failed at MPF2")
        # try:
        #     self.visualize_NMI()
        # except:
        #     print("failed at visualize_NMI")
        # try:
        #     self.visualize_ovr()
        # except:
        #     print("failed at visualize_ovr")
        # try:
        #     self.visualize_AUC_mAUC()
        # except:
        #     print("failed at visualize_AUC_mAUC")
        try:
            self.visualize_robust()
        except:
            print("failed at visualize_robust")
        # try:
        #     self.post_clustering_MI()
        # except:
        #     print("failed at post_clustering_MI")

def merge_WG_subset(dir1, dir2, log_file):
    """
    this function fills the missing analysis from WG runs with their corresponding results
    from subset runs.
    """
    if os.path.exists(log_file):
        prev_log = open(log_file, 'r').readlines()
    else:
        prev_log = []

    # Open the log file for writing
    with open(log_file, 'w') as f:
        for l in prev_log:
            f.write(l)

        # Recursively iterate over all subdirectories and files in dir1
        for root, dirs, files in os.walk(dir1):
            for file in files:
                # Construct the absolute path of the current file
                file_path = os.path.join(root, file)
                # Construct the relative path of the current file with respect to dir1
                rel_path = os.path.relpath(file_path, dir1)
                # Construct the target path in dir2
                target_path = os.path.join(dir2, rel_path)
                # Check if the target path exists
                if not os.path.exists(target_path):
                    # Create any necessary subdirectories in dir2
                    os.makedirs(os.path.dirname(target_path), exist_ok=True)
                    # Copy the file to dir2
                    shutil.copy2(file_path, target_path)
                    # Write a log entry to the log file
                    f.write(f'Copied {file_path} to {target_path}\n')

def exemplar_naive(maindir):

    progress_file = maindir + "/naive_post_clustering_progress.txt"
                   
    if os.path.exists(progress_file):
        with open(progress_file,"r") as f:
            lines = f.readlines()

        naive_rec = ast.literal_eval(lines[0][:-1])
        naive_entropy = ast.literal_eval(lines[1][:-1])

    ####################################################################################
    fig, axs = plt.subplots(1, 1, figsize=(5, 5), sharex=True, sharey=True)

    xs = np.array([int(k) for k in naive_rec.keys()])
    ys = np.array(list(naive_rec.values()))

    axs.scatter(xs, ys, color="blue", alpha=0.7)
    axs.plot(xs, ys, linestyle="solid", color="blue", label="ChromHMM", alpha=0.7)
    axs.set_title("S1: diff. data, diff models | GM12878")
    axs.set_xlabel("Number of Chromatin states")
    axs.set_ylabel("Naive overlap")

    lines, labels = axs.get_legend_handles_labels()

    # Create a dictionary to store unique entries
    unique_entries = {}
    # Iterate over the lines and labels
    for line, label in zip(lines, labels):
        # Only add unique entries to the dictionary
        if label not in unique_entries:
            unique_entries[label] = line

    # Create a custom legend using the unique entries
    plt.legend(unique_entries.values(), unique_entries.keys())
    plt.tight_layout()
    plt.savefig(maindir + "/nlabels_naive_EXEMPLAR.pdf", format="pdf")
    plt.savefig(maindir + "/nlabels_naive_EXEMPLAR.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')
    ####################################################################################
    fig, axs = plt.subplots(1, 1, figsize=(5, 5), sharex=True, sharey=True)

    xs = np.array([float(k) for k in naive_entropy.keys()])
    ys = np.array(list(naive_entropy.values()))

    axs.scatter(xs, ys, color="blue", alpha=0.7)
    axs.plot(xs, ys, linestyle="solid", color="blue", label="ChromHMM", alpha=0.7)
    axs.set_title("S1: diff. data, diff models | GM12878")
    axs.set_xlabel("Entropy")
    axs.set_ylabel("Naive overlap")

    lines, labels = axs.get_legend_handles_labels()

    # Create a dictionary to store unique entries
    unique_entries = {}
    # Iterate over the lines and labels
    for line, label in zip(lines, labels):
        # Only add unique entries to the dictionary
        if label not in unique_entries:
            unique_entries[label] = line

    # Create a custom legend using the unique entries
    plt.legend(unique_entries.values(), unique_entries.keys())
    plt.tight_layout()
    plt.savefig(maindir + "/ent_naive_EXEMPLAR.pdf", format="pdf")
    plt.savefig(maindir + "/ent_naive_EXEMPLAR.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

def coverage_legend():
    import seaborn as sns
    import matplotlib.pyplot as plt

    # example data
    data = [0.01, 0.03, 0.1, 0.4, 0.6, 0.8]
    coverage_quantizer = {(0.0, 0.02):10, (0.02, 0.2):15, (0.2, 0.5):20, (0.5, 1):25}

    # create a new list for dot sizes
    dot_sizes = []
    for value in data:
        for key in coverage_quantizer:
            if key[0] <= value < key[1]:
                dot_sizes.append(coverage_quantizer[key])
                break

    # create the strip plot
    fig, ax = plt.subplots(figsize=(24,18))
    ax.scatter(data, [1] * len(data), s=dot_sizes, color='black', alpha=0.8)

    # create the legend
    legend_elements = []
    for key in coverage_quantizer:
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', label=f'{key[0]} - {key[1]}', markerfacecolor='black', markersize=coverage_quantizer[key], alpha=0.8))

    ax.legend(handles=legend_elements, title='Coverage Range', bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)

    plt.show()

def reprod_vs_exp(maindir):
    "navigate celltypes, settings, and saga model directories"
    "filter the ones with faulty output"
    navigate_results = {}

    for s in [x for x in os.listdir(maindir) if os.path.isdir("{}/{}".format(maindir, x))]:
        navigate_results[s] = {}

        for m in [xx for xx in os.listdir("{}/{}".format(maindir, s)) if os.path.isdir("{}/{}/{}".format(maindir, s, xx))]:
            navigate_results[s][m] = {}

            for c in [xxx for xxx in os.listdir("{}/{}/{}".format(maindir, s, m)) if os.path.isdir("{}/{}/{}/{}".format(maindir, s, m, xxx))]:
                if c in ["GM12878", "K562", "MCF-7"]:
                    navigate_results[s][m][c] = ["{}/{}/{}/{}".format(maindir, s, m, c)]
        
    # print(navigate_results)

    var_setting_dict = {
        "concat":"S2: Diff. data, Same model", "r1vsr2": "S1: Diff. data, Diff. model", 
        "paraminit":"S3: Same data, Diff. model"}

    def load_conf_vs_nonconf(file):
        with open(file, "r") as f:
            lines = f.readlines()
            zp_val = float(lines[-1].split("=")[1])
        return zp_val

    def load_r_v_exp(file):
        pearson = {}
        r2 = {}
        with open(file, "r") as f:
            lines = f.readlines()
            for i in range(len(lines)):
                l = lines[i]
                l = l.split(" = ")
                if "Tran" in l[0]:
                    label, metric, score = l[0].split("|")[0], l[0].split("|")[1], float(l[1])

                    if metric == "Pearson_r":
                        pearson[label] = score
                    elif metric == "R2":
                        r2[label] = score
        
        r2 = max(list(r2.values()))
        pearson = max(list(pearson.values()))
        return r2, pearson

    for s in navigate_results.keys():
        for m in navigate_results[s]:
            for c in navigate_results[s][m]:

                if os.path.exists(navigate_results[s][m][c][0] + "/16_states"):
                    s16_dir = navigate_results[s][m][c][0] + "/16_states/"
                    try:
                        con_v_conf16 = -1 * np.log(
                            load_conf_vs_nonconf(s16_dir + "conf_vs_nonconf_meanEXP_metrics.txt") + 1e-19
                        )
                    except:
                        con_v_conf16 = None

                    try:
                        rve_r2_16, rve_pcc_16 = load_r_v_exp(s16_dir + "mean_expression_v_r_metrics.txt")
                    except:
                        rve_r2_16, rve_pcc_16 = None, None

                    try:
                        rve_r2_16_GB, rve_pcc_16_GB = load_r_v_exp(s16_dir + "mean_expression_v_r_genebody_metrics.txt")
                    except:
                        rve_r2_16_GB, rve_pcc_16_GB = None, None
                    
                    ################################################################################################

                    s14_dir = navigate_results[s][m][c][0] + "/14_states/"
                    try:
                        con_v_conf14 = -1 * np.log(
                            load_conf_vs_nonconf(s14_dir + "conf_vs_nonconf_meanEXP_metrics.txt") + 1e-19
                        )
                    except:
                        con_v_conf14 = None

                    try:
                        rve_r2_14, rve_pcc_14 = load_r_v_exp(s14_dir + "mean_expression_v_r_metrics.txt")
                    except:
                        rve_r2_14, rve_pcc_14 = None, None

                    try:
                        rve_r2_14_GB, rve_pcc_14_GB = load_r_v_exp(s14_dir + "mean_expression_v_r_genebody_metrics.txt")
                    except:
                        rve_r2_14_GB, rve_pcc_14_GB = None, None

                    ################################################################################################

                    s12_dir = navigate_results[s][m][c][0] + "/12_states/"
                    try:
                        con_v_conf12 = -1 * np.log(
                            load_conf_vs_nonconf(s12_dir + "conf_vs_nonconf_meanEXP_metrics.txt") + 1e-19
                        )
                    except:
                        con_v_conf12 = None

                    try:
                        rve_r2_12, rve_pcc_12 = load_r_v_exp(s12_dir + "mean_expression_v_r_metrics.txt")
                    except:
                        rve_r2_12, rve_pcc_12 = None, None

                    try:
                        rve_r2_12_GB, rve_pcc_12_GB = load_r_v_exp(s12_dir + "mean_expression_v_r_genebody_metrics.txt")
                    except:
                        rve_r2_12_GB, rve_pcc_12_GB = None, None

                    ################################################################################################

                    s10_dir = navigate_results[s][m][c][0] + "/10_states/"
                    try:
                        con_v_conf10 = -1 * np.log(
                            load_conf_vs_nonconf(s10_dir + "conf_vs_nonconf_meanEXP_metrics.txt") + 1e-19
                        )
                    except:
                        con_v_conf10 = None

                    try:
                        rve_r2_10, rve_pcc_10 = load_r_v_exp(s10_dir + "mean_expression_v_r_metrics.txt")
                    except:
                        rve_r2_10, rve_pcc_10 = None, None

                    try:
                        rve_r2_10_GB, rve_pcc_10_GB = load_r_v_exp(s10_dir + "mean_expression_v_r_genebody_metrics.txt")
                    except:
                        rve_r2_10_GB, rve_pcc_10_GB = None, None
                    
                    ################################################################################################

                    conf_v_conf = [con_v_conf16, con_v_conf14, con_v_conf12, con_v_conf10]

                    rve_r2 = [rve_r2_16, rve_r2_14, rve_r2_12, rve_r2_10]
                    rve_pcc = [rve_pcc_16, rve_pcc_14, rve_pcc_12, rve_pcc_10]

                    GB_rve_r2 = [rve_r2_16_GB, rve_r2_14_GB, rve_r2_12_GB, rve_r2_10_GB]
                    GB_rve_pcc = [rve_pcc_16_GB, rve_pcc_14_GB, rve_pcc_12_GB, rve_pcc_10_GB]

                    navigate_results[s][m][c].append([conf_v_conf, rve_r2, rve_pcc, GB_rve_r2, GB_rve_pcc])

                else:
                    pass

    # print(navigate_results)
    chmm = []
    segway = []
    for s in navigate_results.keys():
        for m in navigate_results[s]:
            for c in navigate_results[s][m]:

                if m == "chmm":
                    # setting, ct, n_labels, con_v_conf, rve_r2, rve_pcc, GB_rve_r2, GB_rve_pcc
                    chmm.append(
                        [
                            var_setting_dict[s], c, 16, navigate_results[s][m][c][1][0][0], navigate_results[s][m][c][1][1][0], 
                            navigate_results[s][m][c][1][2][0], navigate_results[s][m][c][1][3][0], navigate_results[s][m][c][1][4][0]
                        ]
                    )

                    chmm.append(
                        [
                            var_setting_dict[s], c, 14, navigate_results[s][m][c][1][0][1], navigate_results[s][m][c][1][1][1], 
                            navigate_results[s][m][c][1][2][1], navigate_results[s][m][c][1][3][1], navigate_results[s][m][c][1][4][1]
                        ]
                    )

                    chmm.append(
                        [
                            var_setting_dict[s], c, 12, navigate_results[s][m][c][1][0][2], navigate_results[s][m][c][1][1][2], 
                            navigate_results[s][m][c][1][2][2], navigate_results[s][m][c][1][3][2], navigate_results[s][m][c][1][4][2]
                        ]
                    )

                    chmm.append(
                        [
                            var_setting_dict[s], c, 10, navigate_results[s][m][c][1][0][3], navigate_results[s][m][c][1][1][3], 
                            navigate_results[s][m][c][1][2][3], navigate_results[s][m][c][1][3][3], navigate_results[s][m][c][1][4][3]
                        ]
                    )

                elif m == "segway":
                    segway.append(
                        [
                            var_setting_dict[s], c, 16, navigate_results[s][m][c][1][0][0], navigate_results[s][m][c][1][1][0], 
                            navigate_results[s][m][c][1][2][0], navigate_results[s][m][c][1][3][0], navigate_results[s][m][c][1][4][0]
                        ]
                    )

                    segway.append(
                        [
                            var_setting_dict[s], c, 14, navigate_results[s][m][c][1][0][1], navigate_results[s][m][c][1][1][1], 
                            navigate_results[s][m][c][1][2][1], navigate_results[s][m][c][1][3][1], navigate_results[s][m][c][1][4][1]
                        ]
                    )

                    segway.append(
                        [
                            var_setting_dict[s], c, 12, navigate_results[s][m][c][1][0][2], navigate_results[s][m][c][1][1][2], 
                            navigate_results[s][m][c][1][2][2], navigate_results[s][m][c][1][3][2], navigate_results[s][m][c][1][4][2]
                        ]
                    )

                    segway.append(
                        [
                            var_setting_dict[s], c, 10, navigate_results[s][m][c][1][0][3], navigate_results[s][m][c][1][1][3], 
                            navigate_results[s][m][c][1][2][3], navigate_results[s][m][c][1][3][3], navigate_results[s][m][c][1][4][3]
                        ]
                    )
    
    chmm = pd.DataFrame(
        chmm, 
        columns = ["setting", "CellType", "n_labels", "con_v_conf", "rve_r2", "rve_pcc", "GB_rve_r2", "GB_rve_pcc"]).sort_values(
            by=["setting", "CellType"])

    segway = pd.DataFrame(
        segway, 
        columns = ["setting", "CellType", "n_labels", "con_v_conf", "rve_r2", "rve_pcc", "GB_rve_r2", "GB_rve_pcc"]).sort_values(
            by=["setting", "CellType"])

    def visualize(df, metric):
        fig, axs = plt.subplots(1, 1, figsize=(12, 9), sharex=False, sharey=False)

        # Flatten the axs array to make it easier to iterate over
        # axs = axs.flatten()

        # Define the marker styles you want to use for each category of the 'celltype' variable
        markers = {
            'S3: Same data, Diff. model': 'o', 
            'S2: Diff. data, Same model': 's',
            "S1: Diff. data, Diff. model": "v"} 

        # Plot the data in each subplot
        for setting, marker in markers.items():
            df_s = df[df['setting'] == setting]
            sns.stripplot(x="n_labels", y=metric, hue="CellType", data=df_s, ax=axs, marker=marker, size=10, alpha=0.8)

        for i in range(len(df.n_labels.unique())):
            axs.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

        # axs.set_xticklabels(axs.get_xticklabels(), fontsize=9)
        axs.set_xlabel('Number of States', fontsize=14)

        axs.tick_params(axis='y', labelsize=14)
        axs.tick_params(axis='x', labelsize=14)

        # Get the color palette used by seaborn
        palette = sns.color_palette()

        # # Create a custom legend
        legend_elements = []
        for i, ct in enumerate(df.CellType.unique()):
            legend_elements.append(plt.Line2D([0], [0], color=palette[i], label=ct, marker='o', markersize=10, linestyle='None'))
            
        for i, setting in enumerate(df.setting.unique()):
            legend_elements.append(plt.Line2D([0], [0], color='black', marker=markers[setting], label=setting, linestyle='None'))
    
        #  # Add the custom legend to the first subplot
        axs.legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=9)

    ####################################################################################

    visualize(chmm, "rve_r2")
    plt.ylabel("R2", fontsize=14)
    plt.savefig(maindir+"/chmm_rve_r2.pdf", format="pdf")
    plt.savefig(maindir+"/chmm_rve_r2.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    visualize(chmm, "con_v_conf")
    plt.ylabel("Z-test -log(p)", fontsize=14)
    plt.savefig(maindir+"/chmm_rve_con_v_conf.pdf", format="pdf")
    plt.savefig(maindir+"/chmm_rve_con_v_conf.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    visualize(chmm, "rve_pcc")
    plt.ylabel("PCC", fontsize=14)
    plt.savefig(maindir+"/chmm_rve_pcc.pdf", format="pdf")
    plt.savefig(maindir+"/chmm_rve_pcc.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    ####################################################################################
    visualize(segway, "rve_r2")
    plt.ylabel("R2", fontsize=14)
    plt.savefig(maindir+"/segway_rve_r2.pdf", format="pdf")
    plt.savefig(maindir+"/segway_rve_r2.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    visualize(segway, "con_v_conf")
    plt.ylabel("Z-test -log(p)", fontsize=14)
    plt.savefig(maindir+"/segway_rve_con_v_conf.pdf", format="pdf")
    plt.savefig(maindir+"/segway_rve_con_v_conf.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    visualize(segway, "rve_pcc")
    plt.ylabel("PCC", fontsize=14)
    plt.savefig(maindir+"/segway_rve_pcc.pdf", format="pdf")
    plt.savefig(maindir+"/segway_rve_pcc.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

def reprod_vs_exp_per_label(maindir):
    "navigate celltypes, settings, and saga model directories"
    "filter the ones with faulty output"
    navigate_results = {}

    for s in [x for x in os.listdir(maindir) if os.path.isdir("{}/{}".format(maindir, x))]:
        navigate_results[s] = {}

        for m in [xx for xx in os.listdir("{}/{}".format(maindir, s)) if os.path.isdir("{}/{}/{}".format(maindir, s, xx))]:
            navigate_results[s][m] = {}

            for c in [xxx for xxx in os.listdir("{}/{}/{}".format(maindir, s, m)) if os.path.isdir("{}/{}/{}/{}".format(maindir, s, m, xxx))]:
                if c in ["GM12878", "K562", "MCF-7"]:
                    navigate_results[s][m][c] = ["{}/{}/{}/{}".format(maindir, s, m, c)]
        
    # print(navigate_results)

    var_setting_dict = {
        "concat":"S2: Diff. data, Same model", "r1vsr2": "S1: Diff. data, Diff. model", 
        "paraminit":"S3: Same data, Diff. model"}

    setting_inverse = {v:k for k, v in var_setting_dict.items()}

    def load_r_v_exp(file):
        pearson = {}
        r2 = {}
        with open(file, "r") as f:
            lines = f.readlines()
            for i in range(len(lines)):
                l = lines[i]
                l = l.split(" = ")
                label, metric, score = l[0].split("|")[0], l[0].split("|")[1], float(l[1])

                if metric == "Pearson_r":
                    pearson[label] = score
                elif metric == "R2":
                    r2[label] = score
                    
        return r2, pearson

    for s in navigate_results.keys():
        for m in navigate_results[s]:
            for c in navigate_results[s][m]:

                if os.path.exists(navigate_results[s][m][c][0] + "/16_states"):
                    s16_dir = navigate_results[s][m][c][0] + "/16_states/"

                    try:
                        rve_r2_16, rve_pcc_16 = load_r_v_exp(s16_dir + "mean_expression_v_r_metrics.txt")
                    except:
                        rve_r2_16, rve_pcc_16 = {}, {}

                    try:
                        rve_r2_16_GB, rve_pcc_16_GB = load_r_v_exp(s16_dir + "mean_expression_v_r_genebody_metrics.txt")
                    except:
                        rve_r2_16_GB, rve_pcc_16_GB = {}, {}
                    
                    ################################################################################################

                    s14_dir = navigate_results[s][m][c][0] + "/14_states/"

                    try:
                        rve_r2_14, rve_pcc_14 = load_r_v_exp(s14_dir + "mean_expression_v_r_metrics.txt")
                    except:
                        rve_r2_14, rve_pcc_14 = {}, {}

                    try:
                        rve_r2_14_GB, rve_pcc_14_GB = load_r_v_exp(s14_dir + "mean_expression_v_r_genebody_metrics.txt")
                    except:
                        rve_r2_14_GB, rve_pcc_14_GB = {}, {}

                    ################################################################################################

                    s12_dir = navigate_results[s][m][c][0] + "/12_states/"

                    try:
                        rve_r2_12, rve_pcc_12 = load_r_v_exp(s12_dir + "mean_expression_v_r_metrics.txt")
                    except:
                        rve_r2_12, rve_pcc_12 = {}, {}

                    try:
                        rve_r2_12_GB, rve_pcc_12_GB = load_r_v_exp(s12_dir + "mean_expression_v_r_genebody_metrics.txt")
                    except:
                        rve_r2_12_GB, rve_pcc_12_GB = {}, {}

                    ################################################################################################

                    s10_dir = navigate_results[s][m][c][0] + "/10_states/"

                    try:
                        rve_r2_10, rve_pcc_10 = load_r_v_exp(s10_dir + "mean_expression_v_r_metrics.txt")
                    except:
                        rve_r2_10, rve_pcc_10 = {}, {}

                    try:
                        rve_r2_10_GB, rve_pcc_10_GB = load_r_v_exp(s10_dir + "mean_expression_v_r_genebody_metrics.txt")
                    except:
                        rve_r2_10_GB, rve_pcc_10_GB = {}, {}
                    
                    ################################################################################################

                    rve_r2 = [rve_r2_16, rve_r2_14, rve_r2_12, rve_r2_10]
                    rve_pcc = [rve_pcc_16, rve_pcc_14, rve_pcc_12, rve_pcc_10]

                    GB_rve_r2 = [rve_r2_16_GB, rve_r2_14_GB, rve_r2_12_GB, rve_r2_10_GB]
                    GB_rve_pcc = [rve_pcc_16_GB, rve_pcc_14_GB, rve_pcc_12_GB, rve_pcc_10_GB]

                    navigate_results[s][m][c].append([rve_r2, rve_pcc, GB_rve_r2, GB_rve_pcc])

                else:
                    pass

    chmm = []
    segway = []
    for s in navigate_results.keys():
        for m in navigate_results[s]:
            for c in navigate_results[s][m]:

                if m == "chmm":
                    # setting, ct, n_labels, state, con_v_conf, rve_r2, rve_pcc, GB_rve_r2, GB_rve_pcc
                    for k in navigate_results[s][m][c][1][0][0].keys():
                        if "+" in k:
                            raise
                        try:
                            chmm.append(
                                [
                                    var_setting_dict[s], c, 16, k, navigate_results[s][m][c][1][0][0][k], navigate_results[s][m][c][1][1][0][k], 
                                    navigate_results[s][m][c][1][2][0][k], navigate_results[s][m][c][1][3][0][k]
                                ]
                            )
                        except:
                            pass
                    
                    for k in navigate_results[s][m][c][1][0][1].keys():
                        try:
                            chmm.append(
                                [
                                    var_setting_dict[s], c, 14, k, navigate_results[s][m][c][1][0][1][k], navigate_results[s][m][c][1][1][1][k], 
                                    navigate_results[s][m][c][1][2][1][k], navigate_results[s][m][c][1][3][1][k]
                                ]
                            )
                        except:
                            pass
                    
                    for k in navigate_results[s][m][c][1][0][2].keys():
                        try:
                            chmm.append(
                                [
                                    var_setting_dict[s], c, 12, k, navigate_results[s][m][c][1][0][2][k], navigate_results[s][m][c][1][1][2][k], 
                                    navigate_results[s][m][c][1][2][2][k], navigate_results[s][m][c][1][3][2][k]
                                ]
                            )
                        except:
                            pass
                    
                    for k in navigate_results[s][m][c][1][0][3].keys():
                        try:
                            chmm.append(
                                [
                                    var_setting_dict[s], c, 10, k, navigate_results[s][m][c][1][0][3][k], navigate_results[s][m][c][1][1][3][k], 
                                    navigate_results[s][m][c][1][2][3][k], navigate_results[s][m][c][1][3][3][k]
                                ]
                            )
                        except:
                            pass
                
                #############################################################################################################################

                elif m == "segway":
                    for k in navigate_results[s][m][c][1][0][0].keys():
                        
                        try:
                            segway.append(
                                [
                                    var_setting_dict[s], c, 16, k, navigate_results[s][m][c][1][0][0][k], navigate_results[s][m][c][1][1][0][k], 
                                    navigate_results[s][m][c][1][2][0][k], navigate_results[s][m][c][1][3][0][k]
                                ]
                            )
                        except:
                            pass
                    
                    for k in navigate_results[s][m][c][1][0][1].keys():
                        try:
                            segway.append(
                                [
                                    var_setting_dict[s], c, 14, k, navigate_results[s][m][c][1][0][1][k], navigate_results[s][m][c][1][1][1][k], 
                                    navigate_results[s][m][c][1][2][1][k], navigate_results[s][m][c][1][3][1][k]
                                ]
                            )
                        except:
                            pass
                    
                    for k in navigate_results[s][m][c][1][0][2].keys():
                        try:
                            segway.append(
                                [
                                    var_setting_dict[s], c, 12, k, navigate_results[s][m][c][1][0][2][k], navigate_results[s][m][c][1][1][2][k], 
                                    navigate_results[s][m][c][1][2][2][k], navigate_results[s][m][c][1][3][2][k]
                                ]
                            )
                        except:
                            pass
                    
                    for k in navigate_results[s][m][c][1][0][3].keys():
                        try:
                            segway.append(
                                [
                                    var_setting_dict[s], c, 10, k, navigate_results[s][m][c][1][0][3][k], navigate_results[s][m][c][1][1][3][k], 
                                    navigate_results[s][m][c][1][2][3][k], navigate_results[s][m][c][1][3][3][k]
                                ]
                            )
                        except:
                            pass
   
    chmm = pd.DataFrame(
        chmm, 
        columns = ["Setting", "CellType", "n_labels", "State", "rve_r2", "rve_pcc", "GB_rve_r2", "GB_rve_pcc"]).sort_values(
            by=["Setting", "CellType", "n_labels", "State"]).reset_index(drop=True)

    segway = pd.DataFrame(
        segway, 
        columns = ["Setting", "CellType", "n_labels", "State", "rve_r2", "rve_pcc", "GB_rve_r2", "GB_rve_pcc"]).sort_values(
            by=["Setting", "CellType", "n_labels", "State"]).reset_index(drop=True)

    # print(chmm)
    # print(segway)
    
    def visualize(df, score, cell_type, setting):
        df_cell = df[(df['CellType'] == cell_type) & (df['Setting'] == setting)]
        # Create a color palette based on the 'State' column
        palette = ['red' if 'Tran' in state else 'grey' for state in df_cell['State']]

        # Create a grid of subplots based on 'setting' and 'n_labels'
        g = sns.FacetGrid(df_cell, row="Setting", col="n_labels", hue="CellType", margin_titles=True, sharex=False, sharey=True, height=6, aspect=0.6)

        # Map the data to a bar plot with 'State' on the x-axis and the score on the y-axis
        g.map(lambda x, y, **kwargs: sns.barplot(x=x, y=y, palette=['red' if 'Tran' in state else 'grey' for state in x], **kwargs), "State", score)

        # Rotate x-axis labels
        g.set_xticklabels(rotation=90, fontsize=9)

        # Add a legend
        # g.add_legend()

        g.tight_layout()

    cell_types = chmm['CellType'].unique()
    settings = chmm['Setting'].unique()

    for s in settings:
        for c in cell_types:

            visualize(chmm, "rve_r2", c, s)
            plt.ylabel("R2", fontsize=10)
            plt.savefig(maindir + "/" + setting_inverse[s] + "/chmm/" + c + "/rve_r2.pdf", format="pdf")
            plt.savefig(maindir + "/" + setting_inverse[s] + "/chmm/" + c + "/rve_r2.svg", format="svg")
            plt.clf()
            sns.reset_orig
            plt.style.use('default')


            visualize(segway, "rve_r2", c, s)
            plt.ylabel("R2", fontsize=10)
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/rve_r2.pdf", format="pdf")
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/rve_r2.svg", format="svg")
            plt.clf()
            sns.reset_orig
            plt.style.use('default')

            ####################################################################################

            visualize(segway, "rve_pcc", c, s)
            plt.ylabel("PCC", fontsize=10)
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/rve_pcc.pdf", format="pdf")
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/rve_pcc.svg", format="svg")
            plt.clf()
            sns.reset_orig
            plt.style.use('default')
            

            visualize(segway, "rve_pcc", c, s)
            plt.ylabel("PCC", fontsize=10)
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/rve_pcc.pdf", format="pdf")
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/rve_pcc.svg", format="svg")
            plt.clf()
            sns.reset_orig
            plt.style.use('default')

            ####################################################################################

            visualize(segway, "GB_rve_r2", c, s)
            plt.ylabel("R2", fontsize=10)
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/GB_rve_r2.pdf", format="pdf")
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/GB_rve_r2.svg", format="svg")
            plt.clf()
            sns.reset_orig
            plt.style.use('default')

            visualize(segway, "GB_rve_r2", c, s)
            plt.ylabel("R2", fontsize=10)
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/GB_rve_r2.pdf", format="pdf")
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/GB_rve_r2.svg", format="svg")
            plt.clf()
            sns.reset_orig
            plt.style.use('default')

            ####################################################################################

            visualize(segway, "GB_rve_pcc", c, s)
            plt.ylabel("PCC", fontsize=10)
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/GB_rve_pcc.pdf", format="pdf")
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/GB_rve_pcc.svg", format="svg")
            plt.clf()
            sns.reset_orig
            plt.style.use('default')

            visualize(segway, "GB_rve_pcc", c, s)
            plt.ylabel("PCC", fontsize=10)
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/GB_rve_pcc.pdf", format="pdf")
            plt.savefig(maindir + "/" + setting_inverse[s] + "/segway/" + c + "/GB_rve_pcc.svg", format="svg")
            plt.clf()
            sns.reset_orig
            plt.style.use('default')


    exit()
    ####################################################################################
    
    visualize(chmm, "rve_r2")
    plt.ylabel("R2", fontsize=14)
    plt.savefig(maindir+"/chmm_rve_r2.pdf", format="pdf")
    plt.savefig(maindir+"/chmm_rve_r2.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    visualize(chmm, "con_v_conf")
    plt.ylabel("Z-test -log(p)", fontsize=14)
    plt.savefig(maindir+"/chmm_rve_con_v_conf.pdf", format="pdf")
    plt.savefig(maindir+"/chmm_rve_con_v_conf.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    visualize(chmm, "rve_pcc")
    plt.ylabel("PCC", fontsize=14)
    plt.savefig(maindir+"/chmm_rve_pcc.pdf", format="pdf")
    plt.savefig(maindir+"/chmm_rve_pcc.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    ####################################################################################
    visualize(segway, "rve_r2")
    plt.ylabel("R2", fontsize=14)
    plt.savefig(maindir+"/segway_rve_r2.pdf", format="pdf")
    plt.savefig(maindir+"/segway_rve_r2.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    visualize(segway, "con_v_conf")
    plt.ylabel("Z-test -log(p)", fontsize=14)
    plt.savefig(maindir+"/segway_rve_con_v_conf.pdf", format="pdf")
    plt.savefig(maindir+"/segway_rve_con_v_conf.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')

    visualize(segway, "rve_pcc")
    plt.ylabel("PCC", fontsize=14)
    plt.savefig(maindir+"/segway_rve_pcc.pdf", format="pdf")
    plt.savefig(maindir+"/segway_rve_pcc.svg", format="svg")
    plt.clf()
    sns.reset_orig
    plt.style.use('default')


if __name__=="__main__":
    # exemplar_naive("tests/cedar_runs/chmm/GM12878_R1/")
    # exit()

    # merge_WG_subset("tests/subset", "tests/WG", "tests/copy_log.txt")
    # comp = COMPARATIVE("tests/subset")
    # comp.ALL()
    reprod_vs_exp_per_label("/Users/mforooz/Desktop/research/libbrechteam@sfu/segwayconf/gitrepo/SAGAconf/tests/Rebuttal_runs/rebuttal_subset")
    exit()
    comp = COMPARATIVE("tests/WG")
    comp.ALL()
    

    exit()

    if sys.argv[1] == "all":
        m_p("s1")
        m_p("s2")
        m_p("s3")
        # comp = COMPARATIVE("Subset")
        # comp.ALL()

    else:
        setting = sys.argv[1]
        m_p(setting, nt=10)
        # try:
        #     # comp = COMPARATIVE("runs052023_WG")
        #     # comp.ALL()
        # except:
        #     pass
    