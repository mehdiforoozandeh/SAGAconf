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
create_trackname_assay_file(download_dir)

assays = {}
for ct in CellType_list:
    assays[ct] = read_list_of_assays(download_dir+ct)

print(assays)

for c in CellType_list:
    gather_replicates(celltype_dir=download_dir+c)

if os.path.exists(download_dir+"hg38.chrom.sizes") == False:
    sizes_url = 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes'
    sizes_file_dl_response = requests.get(sizes_url, allow_redirects=True)
    open(download_dir+"hg38.chrom.sizes", 'wb').write(sizes_file_dl_response.content)
    print('downloaded the hg38.chrom.sizes file')

# download chromhmm
chromhmm_url = "http://compbio.mit.edu/ChromHMM/ChromHMM.zip"
if os.path.exists("ChromHMM/") == False:
    print("downloading ChromHMM Software")
    chmm_dl_response = requests.get(chromhmm_url, allow_redirects=True)
    open("ChromHMM.zip", 'wb').write(chmm_dl_response.content)
    os.system("unzip ChromHMM.zip")
    
chmm_output_dir = "chromhmm_runs/"
chmm_input_dir = "chmmfiles/"
print("Preparing ChromHMM inputfiles")
prepare_chmm_inputdata("files/K562/", assertion=False)
print("running chmm rep")
ChromHMM_replicate_runs("chmmfiles/K562/", chmm_output_dir, n_thread='0')
print("running chmm concat")
ChromHMM_concat_runs("chmmfiles/K562/", chmm_output_dir, n_thread='0')
print("running chmm param-init")
ChromHMM_paraminit_runs("chmmfiles/K562/", chmm_output_dir, [5, 7],  n_thread='0')