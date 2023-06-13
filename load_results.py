import os, ast, sys
import pandas as pd

class COMPARATIVE(object):
    def __init__(self, maindir):
        "navigate celltypes, settings, and saga model directories"
        "filter the ones with faulty output"
        self.maindir = maindir
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
            "concat":"S2: Diff. data, Same model", 
            "r1vsr2": "S1: Diff. data, Diff. model", 
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
    
    def load_coverages(self):
        for s in self.navigate_results.keys():
            for m in self.navigate_results[s]:
                for c in self.navigate_results[s][m]:
                    
                    if os.path.exists(self.navigate_results[s][m][c][0] + "/coverage1.txt"):
                        coverage_file = self.navigate_results[s][m][c][0] + "/coverage1.txt"
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
                            del ratio_robust["general"]
                        
                        self.navigate_results[s][m][c].append(ratio_robust)

                    else:
                        self.navigate_results[s][m][c].append({})

    def load_all(self):
        self.compare_NMI()
        self.compare_ratio_raw_overlap()
        self.compare_AUC_mAUC()
        self.compare_ratio_robust()
    
    def save_all(self):

        settings = ['paraminit', 'concat', 'r1vsr2']
        sagas = ['chmm', 'segway']
        cts = ['HeLa-S3', 'GM12878', 'MCF-7', 'CD14', 'K562']


        general_df = [] # settings, saga, ct, NMI_map, NMI_post, ovr_ratio_w0, ovr_ratio_w1000
        perlabel_df = [] # settings, saga, ct, label, auc/mauc, ratio_robust

        for i in settings:
            for j in sagas:
                for k in cts:
                    general_df.append([
                        self.var_setting_dict[i], j, str(k), 
                        res.navigate_results[i][j][k][1][0],
                        res.navigate_results[i][j][k][1][1],
                        res.navigate_results[i][j][k][2][0],
                        res.navigate_results[i][j][k][2][1]])

                    for l in res.navigate_results[i][j][k][3].keys():
                        if l in res.navigate_results[i][j][k][4].keys():
                            perlabel_df.append([
                                self.var_setting_dict[i], j, str(k), l, "_".join(l.split("_")[1:]),
                                res.navigate_results[i][j][k][3][l],
                                res.navigate_results[i][j][k][4][l]
                            ])

        general_df = pd.DataFrame(general_df, columns = ["settings", "saga", "ct", "NMI_map", "NMI_post", "ovr_ratio_w0", "ovr_ratio_w1000"])
        perlabel_df = pd.DataFrame(perlabel_df, columns = ["settings", "saga", "ct", "label", "func", "auc/mauc", "ratio_robust"])

        general_df.to_csv(self.maindir  + "/general_df.csv")
        perlabel_df.to_csv(self.maindir + "/perlabel_df.csv")

if __name__ == "__main__":
    res = COMPARATIVE(sys.argv[1])
    res.load_all()
    res.save_all()