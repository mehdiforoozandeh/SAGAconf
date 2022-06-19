from run import *

download_dir = 'protect_files_/'
segway_dir = 'segway_runs/'
res_dir = 'reprod_results/'

CellType_list = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+ct)]

for ct in CellType_list:
    make_pseudo_replicates(ct, m_p=False)