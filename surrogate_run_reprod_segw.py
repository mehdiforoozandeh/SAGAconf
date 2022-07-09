from run import *
import sys

CellType_list = np.array(
    ['K562', 'MCF-7', 'GM12878', 'HeLa-S3', 'CD14-positive_monocyte'])

# download_dir = 'files/'
segway_dir = 'segway_runs/'
chmm_dir = "chromhmm_runs"
res_dir = 'reprod_results/'
short_res_dir = 'short_reports/'

################################    SHORT RUNS   ###################################

if sys.argv[1] == "short":
    if os.path.exists(short_res_dir) == False:
        os.mkdir(short_res_dir)

    if os.path.exists(short_res_dir+"chmm/") == False:
        os.mkdir(short_res_dir+"chmm/")

    RUN_ALL_REPROD_ANALYSIS(
        chmm_dir,  CellType_list, short_res_dir+"chmm/", type="chmm", multi_p=True, n_processors=5, run="short")
    
    if os.path.exists(short_res_dir+"segway/") == False:
        os.mkdir(short_res_dir+"segway/")

    RUN_ALL_REPROD_ANALYSIS(
        segway_dir,  CellType_list, short_res_dir+"segway/", type="segway", multi_p=True, n_processors=5, run="short")



################################    LONG RUNS   ###################################

elif sys.argv[1] == "long":
    if os.path.exists(res_dir) == False:
        os.mkdir(res_dir)

    if os.path.exists(res_dir+"segway/") == False:
        os.mkdir(res_dir+"segway/")

    RUN_ALL_REPROD_ANALYSIS(
        segway_dir,  CellType_list, res_dir+"segway/", type="segway", multi_p=True, n_processors=10)

    if os.path.exists(res_dir+"chmm/") == False:
        os.mkdir(res_dir+"chmm/")

    RUN_ALL_REPROD_ANALYSIS(
        chmm_dir,  CellType_list, res_dir+"chmm/", type="chmm", multi_p=True, n_processors=10)