from run import *

CellType_list = np.array(
    ['K562', 'MCF-7', 'GM12878', 'HeLa-S3', 'CD14-positive monocyte'])

download_dir = 'files/'
segway_dir = 'segway_runs/'
res_dir = 'reprod_results/'

if os.path.exists(res_dir) == False:
    os.mkdir(res_dir)

print('list of target celltypes', CellType_list)
existing_data = np.array(check_if_data_exists(CellType_list, download_dir))
CellType_list = [CellType_list[i] for i in range(len(CellType_list)) if existing_data[i]==False]

if len(CellType_list) != 0:
    download_encode_files(CellType_list, download_dir, "GRCh38")
else:
    print('No download required!')

CellType_list = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+ct)]

# clean up potential space characters in directory names to prevent later issues
for ct in CellType_list:
    if " " in ct:
        os.system("mv {} {}".format(
            ct.replace(' ', '\ '), ct.replace(" ", "_")
        ))

CellType_list = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+ct)]

if os.path.exists(res_dir+"segway/") == False:
    os.mkdir(res_dir+"segway/")
    
RUN_ALL_REPROD_ANALYSIS(
    segway_dir,  CellType_list, res_dir+"segway/", type="segway", n_processors=8)