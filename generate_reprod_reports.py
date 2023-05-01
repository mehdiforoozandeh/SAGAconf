from reports import *
import multiprocessing as mp
import logging, ast, sys

def r1vsr2(maindir="runs052023_subset"):
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

def concat(maindir="runs052023_subset"):
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

def paraminit(maindir="runs052023_subset"):
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
            contour=True
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
        
        self.var_setting_dict_inverse = {v:k for k, v in self.var_setting_dict.items()}
        self.SORT_ORDER = {"Prom": 0, "Prom_fla":1, "Enha":2, "Enha_low":3, "Biva":4, "Tran":5, "Cons":6, "Facu":7, "K9K3":8, "Quie":9}

    # self.navigate_results[s][m][c] = [dir, [NMI_map, NMI_post], [ovr_ratio_w0, ovr_ratio_w1000],  {auc/mauc}, {ratio_robust}]
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
    
    def compare_ratio_robust(self):
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s]:
                for c in self.navigate_results[s][m]:
                    
                    if os.path.exists(self.navigate_results[s][m][c][0] + "/general_reproducibility_score_POSTERIOR.txt"):
                        dir_robust_file = self.navigate_results[s][m][c][0] + "/general_reproducibility_score_POSTERIOR.txt"
                        ratio_robust = {}

                        rlines = open(dir_robust_file, "r").readlines()
                        for l in rlines:
                            ratio_robust[l.split(" = ")[0].replace(" reprod score", "")] = float(l.split(" = ")[1][:7])

                        if "general" in ratio_robust.keys():
                            del ratio_robust["general"]
                        
                        self.navigate_results[s][m][c].append(ratio_robust)

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
                        smc_list = []
                        for k in smc_dict.keys():
                            new_k = "_".join(k.split("_")[1:])
                            smc_list.append([new_k, smc_dict[k]])

                        segway[s][c] = smc_list

                    elif m == "chmm":
                        if s not in chmm.keys():
                            chmm[s] = {}
                        
                        smc_dict = self.navigate_results[s][m][c][3]
                        smc_list = []
                        for k in smc_dict.keys():
                            new_k = "_".join(k.split("_")[1:])
                            smc_list.append([new_k, smc_dict[k]])

                        chmm[s][c] = smc_list

        df_chmm = []
        for s in chmm.keys():
            for c in chmm[s].keys():
                for l in chmm[s][c]:
                    if l[1] > 0:
                        df_chmm.append([c, self.var_setting_dict[s], l[0], l[1]])
        
        df_chmm = pd.DataFrame(
            df_chmm, columns = ["CellType", "Setting", "Genomic_function", "AUC/maxAUC"]).sort_values(
                by=["Setting", "CellType"])

        df_segway = []
        for s in segway.keys():
            for c in segway[s].keys():
                for l in segway[s][c]:
                    if l[1] > 0:
                        df_segway.append([c, self.var_setting_dict[s], l[0], l[1]])
        
        df_segway = pd.DataFrame(
            df_segway, columns = ["CellType", "Setting", "Genomic_function", "AUC/maxAUC"]).sort_values(
                by=["Setting", "CellType"])


        ####################################################################################################################################
        for s in np.unique(df_chmm["Setting"]):
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.set_palette(sns.color_palette("deep"))
            # sns.set_theme(style="whitegrid")
            
            sns.stripplot(x='Genomic_function', y="AUC/maxAUC", hue="CellType", data=df_chmm.loc[df_chmm["Setting"]==s,:], ax=ax, size=10)
            # sns.boxplot(x='Genomic_function', y="AUC/maxAUC", hue="CellType", data=df_chmm.loc[df_chmm["Setting"]==s,:], ax=ax)
            n_categories = len(ax.get_xticks())

            for i in range(n_categories):
                ax.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

            ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)

            # Show the plot
            plt.yticks(np.arange(0.5, 1.05, step=0.1))
            plt.tight_layout()
            plt.savefig("{}/{}/{}/AUC_mAUC_stripplot.pdf".format(self.maindir, self.var_setting_dict_inverse[s], "chmm"), format="pdf")
            plt.savefig("{}/{}/{}/AUC_mAUC_stripplot.svg".format(self.maindir, self.var_setting_dict_inverse[s], "chmm"), format="svg")

            plt.clf()
            sns.reset_orig
            plt.style.use('default')


        ####################################################################################################################################
        for s in np.unique(df_segway["Setting"]):
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.set_palette(sns.color_palette("deep"))
            # sns.set_theme(style="whitegrid")
            
            sns.stripplot(x='Genomic_function', y="AUC/maxAUC", hue="CellType", data=df_segway.loc[df_segway["Setting"]==s,:], ax=ax, size=10)
            # sns.boxplot(x='Genomic_function', y="AUC/maxAUC", hue="CellType", data=df_chmm.loc[df_chmm["Setting"]==s,:], ax=ax)
            n_categories = len(ax.get_xticks())

            for i in range(n_categories):
                ax.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

            ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)

            # Show the plot
            plt.yticks(np.arange(0.5, 1.05, step=0.1))
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

                        smc_dict = self.navigate_results[s][m][c][4]
                        smc_list = []
                        for k in smc_dict.keys():
                            new_k = "_".join(k.split("_")[1:])
                            smc_list.append([new_k, smc_dict[k]])

                        segway[s][c] = smc_list

                    elif m == "chmm":
                        if s not in chmm.keys():
                            chmm[s] = {}
                        
                        smc_dict = self.navigate_results[s][m][c][4]
                        smc_list = []
                        for k in smc_dict.keys():
                            new_k = "_".join(k.split("_")[1:])
                            smc_list.append([new_k, smc_dict[k]])

                        chmm[s][c] = smc_list

        df_chmm = []
        for s in chmm.keys():
            for c in chmm[s].keys():
                for l in chmm[s][c]:
                    if l[1] > 0:
                        df_chmm.append([c, self.var_setting_dict[s], l[0], l[1]])
        
        df_chmm = pd.DataFrame(
            df_chmm, columns = ["CellType", "Setting", "Genomic_function", "ratio_robust"]).sort_values(
                by=["Setting", "CellType"])
        
        df_chmm['sort_col'] = df_chmm["Genomic_function"].map(self.SORT_ORDER)
        df_chmm.sort_values("sort_col", inplace=True)
        df_chmm.drop("sort_col", axis=1, inplace=True)

        df_segway = []
        for s in segway.keys():
            for c in segway[s].keys():
                for l in segway[s][c]:
                    if l[1] > 0:
                        df_segway.append([c, self.var_setting_dict[s], l[0], l[1]])
        
        df_segway = pd.DataFrame(
            df_segway, columns = ["CellType", "Setting", "Genomic_function", "ratio_robust"]).sort_values(
                by=["Setting", "CellType"])
        
        df_segway['sort_col'] = df_segway["Genomic_function"].map(self.SORT_ORDER)
        df_segway.sort_values("sort_col", inplace=True)
        df_segway.drop("sort_col", axis=1, inplace=True)

        ####################################################################################################################################
        for s in np.unique(df_chmm["Setting"]):
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.set_palette(sns.color_palette("deep"))
            # sns.set_theme(style="whitegrid")
            
            sns.stripplot(x='Genomic_function', y="ratio_robust", hue="CellType", data=df_chmm.loc[df_chmm["Setting"]==s,:], ax=ax, size=10)
            # sns.boxplot(x='Genomic_function', y="AUC/maxAUC", hue="CellType", data=df_chmm.loc[df_chmm["Setting"]==s,:], ax=ax)
            n_categories = len(ax.get_xticks())

            for i in range(n_categories):
                ax.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

            ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)

            # Show the plot
            plt.yticks(np.arange(0, 1.05, step=0.1))
            plt.tight_layout()
            plt.savefig("{}/{}/{}/ratio_robust_stripplot.pdf".format(self.maindir, self.var_setting_dict_inverse[s], "chmm"), format="pdf")
            plt.savefig("{}/{}/{}/ratio_robust_stripplot.svg".format(self.maindir, self.var_setting_dict_inverse[s], "chmm"), format="svg")

            plt.clf()
            sns.reset_orig
            plt.style.use('default')


        ####################################################################################################################################
        for s in np.unique(df_segway["Setting"]):
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.set_palette(sns.color_palette("deep"))
            # sns.set_theme(style="whitegrid")
            
            sns.stripplot(x='Genomic_function', y="ratio_robust", hue="CellType", data=df_segway.loc[df_segway["Setting"]==s,:], ax=ax, size=10)
            # sns.boxplot(x='Genomic_function', y="AUC/maxAUC", hue="CellType", data=df_chmm.loc[df_chmm["Setting"]==s,:], ax=ax)
            n_categories = len(ax.get_xticks())

            for i in range(n_categories):
                ax.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

            ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)

            # Show the plot
            plt.yticks(np.arange(0, 1.05, step=0.1))
            plt.tight_layout()
            plt.savefig("{}/{}/{}/ratio_robust_stripplot.pdf".format(self.maindir, self.var_setting_dict_inverse[s], "segway"), format="pdf")
            plt.savefig("{}/{}/{}/ratio_robust_stripplot.svg".format(self.maindir, self.var_setting_dict_inverse[s], "segway"), format="svg")

            plt.clf()
            sns.reset_orig
            plt.style.use('default')

    def ALL(self):
        self.compare_NMI()
        self.compare_ratio_raw_overlap()
        self.compare_AUC_mAUC()
        self.compare_ratio_robust()
        self.visualize_NMI()
        self.visualize_ovr()
        self.visualize_AUC_mAUC()
        self.visualize_robust()
        
if __name__=="__main__":
    setting = sys.argv[1]
    m_p(setting)
    try:
        comp = COMPARATIVE("runs052023_subset")
        comp.ALL()
    except:
        pass
    