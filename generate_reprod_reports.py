from reports import *
import multiprocessing as mp
import logging, ast, sys

def r1vsr2(maindir="WG"):
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

def concat(maindir="WG"):
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

def paraminit(maindir="WG"):
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
        
        self.coverage_quantizer = {(0.0, 0.02):7, (0.02, 0.2):11, (0.2, 0.5):15, (0.5, 1):20}

        self.var_setting_dict_inverse = {v:k for k, v in self.var_setting_dict.items()}
        self.SORT_ORDER = {"Prom": 0, "Prom_fla":1, "Enha":2, "Enha_low":3, "Biva":4, "Tran":5, "Cons":6, "Facu":7, "K9K3":8, "Quie":9}

     # self.navigate_results[s][m][c] = [dir, [NMI_map, NMI_post], [ovr_ratio_w0, ovr_ratio_w1000], {auc/mauc}, general_rep, {ratio_robust}, {coverage}]

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
                            self.navigate_results[s][m][c].append(ratio_robust["general"])
                            del ratio_robust["general"]
                        
                        self.navigate_results[s][m][c].append(ratio_robust)

                    else:
                        self.navigate_results[s][m][c].append(None)
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
            ax.legend(handles, unique_celltypes, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)

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
            ax.legend(handles, unique_celltypes, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)

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

                        smc_dict = self.navigate_results[s][m][c][5]
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
                        
                        smc_dict = self.navigate_results[s][m][c][5]
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
            df_chmm, columns = ["CellType", "Setting", "Genomic_function", "ratio_robust", "coverage"]).sort_values(
                by=["Setting", "CellType"])
        
        df_chmm['sort_col'] = df_chmm["Genomic_function"].map(self.SORT_ORDER)
        df_chmm.sort_values("sort_col", inplace=True)
        df_chmm.drop("sort_col", axis=1, inplace=True)

        df_segway = []
        for s in segway.keys():
            for c in segway[s].keys():
                for l in segway[s][c]:
                    if l[1] > 0:
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
            ax.legend(handles, unique_celltypes, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)

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
            ax.legend(handles, unique_celltypes, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5)

            # Show the plot
            plt.yticks(np.arange(0, 1.05, step=0.1))
            plt.tight_layout()
            plt.savefig("{}/{}/{}/ratio_robust_stripplot.pdf".format(self.maindir, self.var_setting_dict_inverse[s], "segway"), format="pdf")
            plt.savefig("{}/{}/{}/ratio_robust_stripplot.svg".format(self.maindir, self.var_setting_dict_inverse[s], "segway"), format="svg")

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
            sns.stripplot(x="setting", y="NMI", hue="celltype", data=df, ax=axs[0], marker=marker, size=8, alpha=0.8)
            sns.stripplot(x="setting", y="pNMI", hue="celltype", data=df, ax=axs[1], marker=marker, size=8, alpha=0.8)
            sns.stripplot(x="setting", y="overlap | w=0", hue="celltype", data=df, ax=axs[2], marker=marker, size=8, alpha=0.8)
            sns.stripplot(x="setting", y="overlap | w=1000", hue="celltype", data=df, ax=axs[3], marker=marker, size=8, alpha=0.8)

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
            sns.stripplot(x="setting", y="NMI", hue="saga", data=df, ax=axs[0], marker=marker, size=8, alpha=0.8)
            sns.stripplot(x="setting", y="pNMI", hue="saga", data=df, ax=axs[1], marker=marker, size=8, alpha=0.8)
            sns.stripplot(x="setting", y="overlap | w=0", hue="saga", data=df, ax=axs[2], marker=marker, size=8, alpha=0.8)
            sns.stripplot(x="setting", y="overlap | w=1000", hue="saga", data=df, ax=axs[3], marker=marker, size=8, alpha=0.8)

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
        for i, saga in enumerate(general_DF.saga.unique()):
            legend_elements.append(plt.Line2D([0], [0], color=palette[i], label=saga))
            
        for i, celltype in enumerate(general_DF.celltype.unique()):
            legend_elements.append(plt.Line2D([0], [0], color='black', marker=markers[celltype], label=celltype, linestyle='None'))
    
         # Add the custom legend to the first subplot
        axs[0].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=7)

        axs[1].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=7)

        axs[2].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=7)

        axs[3].legend(handles=legend_elements, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=7)

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
            sns.stripplot(x="setting", y="ratio_confident", hue="saga", data=df, ax=axs, marker=marker, size=8, alpha=0.8)

        for i in range(len(general_DF.setting.unique())):
            axs.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

        axs.set_xticklabels(axs.get_xticklabels(), fontsize=9)
        axs.set_xlabel('')

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
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=7)
        
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
            sns.stripplot(x="setting", y="average AUC/mAUC", hue="saga", data=df, ax=axs, marker=marker, size=8, alpha=0.8)

        for i in range(len(general_DF.setting.unique())):
            axs.axvline(i - 0.5, color="gray", linestyle="--", alpha=0.5)

        axs.set_xticklabels(axs.get_xticklabels(), fontsize=9)
        axs.set_xlabel('')

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
                mode="expand", borderaxespad=0, ncol=len(legend_elements), fontsize=7)

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
            self.integ_ratio_robust()
        except:
            print("failed at integ_ratio_robust")
        try:
            self.integ_auc_mauc()
        except:
            print("failed at integ_auc_mauc")
        try:
            self.MPF1()
        except:
            print("failed at MPF1")
        try:
            self.MPF2()
        except:
            print("failed at MPF2")
        try:
            self.visualize_NMI()
        except:
            print("failed at visualize_NMI")
        try:
            self.visualize_ovr()
        except:
            print("failed at visualize_ovr")
        try:
            self.visualize_AUC_mAUC()
        except:
            print("failed at visualize_AUC_mAUC")
        try:
            self.visualize_robust()
        except:
            print("failed at visualize_robust")

if __name__=="__main__":
    # comp = COMPARATIVE("tests/WG")
    # comp.ALL()
    # comp = COMPARATIVE("tests/runs062023_subset")
    # comp.ALL()

    # exit()

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
    