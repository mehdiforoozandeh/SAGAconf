from run import *

CellType_list = np.array(
    ['K562', 'MCF-7', 'GM12878', 'HeLa-S3', 'CD14-positive monocyte'])

download_dir = 'protect_files_/'
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

# convert all bigwigs to bedgraphs (for segway)
for k, v in assays.items():
    for t in v:
        Convert_all_BW2BG(download_dir+k+'/'+t)

# metadata = read_metadata(download_dir)
for c in CellType_list:
    gather_replicates(celltype_dir=download_dir+c)


# download chromosome sizes file for hg38
if os.path.exists(download_dir+"hg38.chrom.sizes") == False:
    sizes_url = 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes'
    sizes_file_dl_response = requests.get(sizes_url, allow_redirects=True)
    open(download_dir+"hg38.chrom.sizes", 'wb').write(sizes_file_dl_response.content)
    print('downloaded the hg38.chrom.sizes file')

# check for existence of genomedata files
gd_exists = []
for ct in CellType_list:
    if os.path.exists(download_dir+ct+'/rep1.genomedata') == False or \
        os.path.exists(download_dir+ct+'/rep2.genomedata') == False:
        gd_exists.append(False)
    
    else:
        gd_exists.append(True)

gd_to_create = [CellType_list[i] for i in range(len(CellType_list)) if gd_exists[i]==False]

if len(gd_to_create) != 0:
    p_obj = mp.Pool(len(gd_to_create))
    p_obj.map(partial(
        create_genomedata, sequence_file=download_dir+"hg38.chrom.sizes"), 
        [download_dir + ct for ct in gd_to_create])

if os.path.exists(segway_dir) == False:
    os.mkdir(segway_dir)

partial_runs_ii = partial(
    RunParse_segway_param_init, 
    replicate_number = 'rep1', output_dir=segway_dir, random_seeds=[5, 7])

p_obj = mp.Pool(len(CellType_list))
p_obj.map(partial_runs_ii, [download_dir+ct for ct in CellType_list])

partial_runs_iii = partial(
    RunParse_segway_param_init, 
    replicate_number = 'rep2', output_dir=segway_dir, random_seeds=[5, 7])
    
p_obj = mp.Pool(len(CellType_list))
p_obj.map(partial_runs_iii, [download_dir+ct for ct in CellType_list])

print('Checking for unparsed posteriors...')
list_of_seg_runs = [
    d for d in os.listdir(segway_dir) if os.path.isdir(segway_dir+'/'+d)]
print(list_of_seg_runs)

for d in list_of_seg_runs:
    print('-Checking for {}  ...'.format(segway_dir+'/'+d+'/parsed_posterior.csv'))

    if os.path.exists(segway_dir+'/'+d+'/parsed_posterior.csv') == False:
        parse_posterior_results(segway_dir+'/'+d, 100, mp=False)

    else:
        print('-Exists!')

print('All parsed!')