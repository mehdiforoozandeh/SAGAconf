from run import *

download_dir = '_protect_files_/'
segway_dir = 'segway_runs/'
res_dir = 'reprod_results/'

CellType_list = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+ct)]

for ct in CellType_list:
    make_pseudo_replicates(download_dir+ct, m_p=True)

# gd_exists = []
# for ct in CellType_list:
#     if os.path.exists(download_dir+ct+'/rep1_psdrep1.genomedata') == False or \
#         os.path.exists(download_dir+ct+'/rep1_psdrep2.genomedata') == False or \
#             os.path.exists(download_dir+ct+'/rep2_psdrep1.genomedata') == False or \
#                 os.path.exists(download_dir+ct+'/rep2_psdrep2.genomedata') == False:

#         gd_exists.append(False)
    
#     else:
#         gd_exists.append(True)

# gd_to_create = [CellType_list[i] for i in range(len(CellType_list)) if gd_exists[i]==False]

# if len(gd_to_create) != 0:
#     p_obj = mp.Pool(len(gd_to_create))
#     p_obj.map(partial(
#         create_psdrep_genomedata, sequence_file=download_dir+"hg38.chrom.sizes"), 
#         [download_dir + ct for ct in gd_to_create])

# if os.path.exists(segway_dir) == False:
#     os.mkdir(segway_dir)


# partial_runs_ii = partial(
#     RunParse_segway_psdreps, 
#     output_dir=segway_dir)

# p_obj = mp.Pool(len(CellType_list))
# p_obj.map(partial_runs_ii, [download_dir+ct for ct in CellType_list])

# print('Checking for unparsed posteriors...')
# list_of_seg_runs = [
#     d for d in os.listdir(segway_dir) if os.path.isdir(segway_dir+'/'+d)]
# print(list_of_seg_runs)

# for d in list_of_seg_runs:
#     print('-Checking for {}  ...'.format(segway_dir+'/'+d+'/parsed_posterior.csv'))

#     if os.path.exists(segway_dir+'/'+d+'/parsed_posterior.csv') == False:
#         parse_posterior_results(segway_dir+'/'+d, 100, mp=False)

#     else:
#         print('-Exists!')

# print('All parsed!')