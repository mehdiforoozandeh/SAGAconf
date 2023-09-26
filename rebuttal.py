import os

def get_listofruns(maindir="rebuttal"):
    listofruns = [
        {"replicate_1_dir":"segway_runs/GM12878_rep1/", 
        "replicate_2_dir":"segway_runs/GM12878_rep2/", 
        "savedir":"{}/r1vsr2/segway/GM12878/".format(maindir)},
        ######## MCF-7 ########
    
        {"replicate_1_dir":"chromhmm_runs/MCF-7_rep1/", 
        "replicate_2_dir":"chromhmm_runs/MCF-7_rep2/", 
        "savedir":"{}/r1vsr2/chmm/MCF-7/".format(maindir)},

        {"replicate_1_dir":"segway_runs/MCF-7_rep1/", 
        "replicate_2_dir":"segway_runs/MCF-7_rep2/", 
        "savedir":"{}/r1vsr2/segway/MCF-7/".format(maindir)},

        ######## CD14 ########
    
        {"replicate_1_dir":"chromhmm_runs/CD14-positive_monocyte_rep1/", 
        "replicate_2_dir":"chromhmm_runs/CD14-positive_monocyte_rep2/", 
        "savedir":"{}/r1vsr2/chmm/CD14/".format(maindir)},
    
        {"replicate_1_dir":"segway_runs/CD14-positive_monocyte_rep1/", 
        "replicate_2_dir":"segway_runs/CD14-positive_monocyte_rep2/", 
        "savedir":"{}/r1vsr2/segway/CD14/".format(maindir)},
        
        ######## K562 ########
    
        {"replicate_1_dir":"chromhmm_runs/K562_rep1/", 
        "replicate_2_dir":"chromhmm_runs/K562_rep2/", 
        "savedir":"{}/r1vsr2/chmm/K562/".format(maindir)},
    
        {"replicate_1_dir":"segway_runs/K562_rep1/", 
        "replicate_2_dir":"segway_runs/K562_rep2/", 
        "savedir":"{}/r1vsr2/segway/K562/".format(maindir)},

        ######## HeLa-S3 ########
    
        {"replicate_1_dir":"chromhmm_runs/HeLa-S3_rep1/", 
        "replicate_2_dir":"chromhmm_runs/HeLa-S3_rep2/", 
        "savedir":"{}/r1vsr2/chmm/HeLa-S3/".format(maindir)},

        {"replicate_1_dir":"segway_runs/HeLa-S3_rep1/", 
        "replicate_2_dir":"segway_runs/HeLa-S3_rep2/", 
        "savedir":"{}/r1vsr2/segway/HeLa-S3/".format(maindir)},

        ######## GM12878 ########
        {"replicate_1_dir":"chromhmm_runs/GM12878_concat/", 
        "replicate_2_dir":"chromhmm_runs/GM12878_concat/", 
        "savedir":"{}/concat/chmm/GM12878/".format(maindir)},

        {"replicate_1_dir":"segway_runs/GM12878_concat_rep1/", 
        "replicate_2_dir":"segway_runs/GM12878_concat_rep2/", 
        "savedir":"{}/concat/segway/GM12878/".format(maindir)},
            
        ######## MCF-7 ########

        {"replicate_1_dir":"chromhmm_runs/MCF-7_concat/", 
        "replicate_2_dir":"chromhmm_runs/MCF-7_concat/", 
        "savedir":"{}/concat/chmm/MCF-7/".format(maindir)},

        {"replicate_1_dir":"segway_runs/MCF-7_concat_rep1/", 
        "replicate_2_dir":"segway_runs/MCF-7_concat_rep2/", 
        "savedir":"{}/concat/segway/MCF-7/".format(maindir)},

        ######## CD14 ########

        {"replicate_1_dir":"chromhmm_runs/CD14-positive_monocyte_concat/", 
        "replicate_2_dir":"chromhmm_runs/CD14-positive_monocyte_concat/", 
        "savedir":"{}/concat/chmm/CD14/".format(maindir)},

        {"replicate_1_dir":"segway_runs/CD14-positive_monocyte_concat_rep1/", 
        "replicate_2_dir":"segway_runs/CD14-positive_monocyte_concat_rep2/", 
        "savedir":"{}/concat/segway/CD14/".format(maindir)},
        
        ######## K562 ########

        {"replicate_1_dir":"chromhmm_runs/K562_concat/", 
        "replicate_2_dir":"chromhmm_runs/K562_concat/", 
        "savedir":"{}/concat/chmm/K562/".format(maindir)},

        {"replicate_1_dir":"segway_runs/K562_concat_rep1/", 
        "replicate_2_dir":"segway_runs/K562_concat_rep2/", 
        "savedir":"{}/concat/segway/K562/".format(maindir)},

        ######## HeLa-S3 ########

        {"replicate_1_dir":"chromhmm_runs/HeLa-S3_concat/", 
        "replicate_2_dir":"chromhmm_runs/HeLa-S3_concat/", 
        "savedir":"{}/concat/chmm/HeLa-S3/".format(maindir)},

        {"replicate_1_dir":"segway_runs/HeLa-S3_concat_rep1/", 
        "replicate_2_dir":"segway_runs/HeLa-S3_concat_rep2/", 
        "savedir":"{}/concat/segway/HeLa-S3/".format(maindir)},
        ######## GM12878 ########

        {"replicate_1_dir":"chromhmm_runs/GM12878_rep2_rs5/", 
        "replicate_2_dir":"chromhmm_runs/GM12878_rep2_rs27/", 
        "savedir":"{}/paraminit/chmm/GM12878/".format(maindir)},


        {"replicate_1_dir":"segway_runs/GM12878_rep1/", 
        "replicate_2_dir":"segway_runs/GM12878_rep1_rs7/", 
        "savedir":"{}/paraminit/segway/GM12878/".format(maindir)},
            
        ######## MCF-7 ########

        {"replicate_1_dir":"chromhmm_runs/MCF-7_rep1_rs5/", 
        "replicate_2_dir":"chromhmm_runs/MCF-7_rep1_rs27/", 
        "savedir":"{}/paraminit/chmm/MCF-7/".format(maindir)},


        {"replicate_1_dir":"segway_runs/MCF-7_rep1_rs5/", 
        "replicate_2_dir":"segway_runs/MCF-7_rep1_rs7/", 
        "savedir":"{}/paraminit/segway/MCF-7/".format(maindir)},

        ######## CD14 ########

        {"replicate_1_dir":"chromhmm_runs/CD14-positive_monocyte_rep1_rs5/", 
        "replicate_2_dir":"chromhmm_runs/CD14-positive_monocyte_rep1_rs27/", 
        "savedir":"{}/paraminit/chmm/CD14/".format(maindir)},


        {"replicate_1_dir":"segway_runs/CD14-positive_monocyte_rep1_rs5/", 
        "replicate_2_dir":"segway_runs/CD14-positive_monocyte_rep1_rs7/", 
        "savedir":"{}/paraminit/segway/CD14/".format(maindir)},
        
        ######## K562 ########

        {"replicate_1_dir":"chromhmm_runs/K562_rep2_rs5/", 
        "replicate_2_dir":"chromhmm_runs/K562_rep2_rs27/", 
        "savedir":"{}/paraminit/chmm/K562/".format(maindir)},

        {"replicate_1_dir":"segway_runs/K562_rep1/", 
        "replicate_2_dir":"segway_runs/K562_rep1_rs7/", 
        "savedir":"{}/paraminit/segway/K562/".format(maindir)},

        ######## HeLa-S3 ########

        {"replicate_1_dir":"chromhmm_runs/HeLa-S3_rep2_rs5/", 
        "replicate_2_dir":"chromhmm_runs/HeLa-S3_rep2_rs27/", 
        "savedir":"{}/paraminit/chmm/HeLa-S3/".format(maindir)},

        {"replicate_1_dir":"segway_runs/HeLa-S3_rep1/", 
        "replicate_2_dir":"segway_runs/HeLa-S3_rep1_rs7/", 
        "savedir":"{}/paraminit/segway/HeLa-S3/".format(maindir)}
        ]
    return listofruns

def get_runs(maindir = "rebuttal"):
    list_of_runs = get_listofruns(maindir)

    if os.path.exists(maindir)==False:
        os.mkdir(maindir)

    if os.path.exists(maindir+"/r1vsr2")==False:
        os.mkdir(maindir+"/r1vsr2")

    if os.path.exists(maindir+"/r1vsr2/chmm")==False:
        os.mkdir(maindir+"/r1vsr2/chmm")

    if os.path.exists(maindir+"/r1vsr2/segway")==False:
        os.mkdir(maindir+"/r1vsr2/segway")
    
    if os.path.exists(maindir+"/concat")==False:
        os.mkdir(maindir+"/concat")
    
    if os.path.exists(maindir+"/concat/chmm")==False:
        os.mkdir(maindir+"/concat/chmm")

    if os.path.exists(maindir+"/concat/segway")==False:
        os.mkdir(maindir+"/concat/segway")

    if os.path.exists(maindir+"/paraminit")==False:
        os.mkdir(maindir+"/paraminit")
    
    if os.path.exists(maindir+"/paraminit/chmm")==False:
        os.mkdir(maindir+"/paraminit/chmm")

    if os.path.exists(maindir+"/paraminit/segway")==False:
        os.mkdir(maindir+"/paraminit/segway")


    for r in list_of_runs:
        savedir = r["savedir"]

        if "chromhmm_runs" in r["replicate_1_dir"] and "concat" in r["replicate_1_dir"]:
            base_mnemonics = r["replicate_1_dir"] + "/mnemonics_rep1.txt"
            verif_mnemonics = r["replicate_2_dir"] + "/mnemonics_rep2.txt"

            replicate_1_dir = r["replicate_1_dir"] + "/parsed_posterior_rep1.csv"
            replicate_2_dir = r["replicate_2_dir"] + "/parsed_posterior_rep2.csv"
        else:
            base_mnemonics = r["replicate_1_dir"] + "/mnemonics.txt"
            verif_mnemonics = r["replicate_2_dir"] + "/mnemonics.txt"

            replicate_1_dir = r["replicate_1_dir"] + "/parsed_posterior.csv"
            replicate_2_dir = r["replicate_2_dir"] + "/parsed_posterior.csv"

        os.system(f"python SAGAconf.py --r_only -s -tr 0.9 -bm {base_mnemonics} -vm {verif_mnemonics} {replicate_1_dir} {replicate_2_dir} {savedir}")

if __name__ == "__main__":
    get_runs(maindir = "rebuttal")