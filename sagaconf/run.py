import numpy as np
import pandas as pd


def intersect_parsed_posteriors(parsed_df_dir_1, parsed_df_dir_2):
    is_chmm_concat = bool(
        ("chromhmm_runs" in parsed_df_dir_1) and ("chromhmm_runs" in parsed_df_dir_2) and 
        ("concat" in parsed_df_dir_1) and ("concat" in parsed_df_dir_2))

    if is_chmm_concat:
        parsed_df_dir_1 = parsed_df_dir_1.replace("parsed_posterior.csv", "parsed_posterior_rep1.csv")
        parsed_df_dir_2 = parsed_df_dir_2.replace("parsed_posterior.csv", "parsed_posterior_rep2.csv")
    
    if ".bed" in parsed_df_dir_1.lower():
        df1 = pd.read_csv(parsed_df_dir_1, sep="\t", on_bad_lines="skip", encoding_errors="ignore")
    elif ".csv" in parsed_df_dir_1.lower():
        df1 = pd.read_csv(parsed_df_dir_1, on_bad_lines="skip", encoding_errors="ignore").drop("Unnamed: 0", axis=1)

    if ".bed" in parsed_df_dir_2.lower():
        df2 = pd.read_csv(parsed_df_dir_2, sep="\t", on_bad_lines="skip", encoding_errors="ignore")
    elif ".csv" in parsed_df_dir_2.lower():
        df2 = pd.read_csv(parsed_df_dir_2, on_bad_lines="skip", encoding_errors="ignore").drop("Unnamed: 0", axis=1)


    df1.iloc[:, 3:] = df1.iloc[:, 3:].astype("float16")
    df2.iloc[:, 3:] = df2.iloc[:, 3:].astype("float16")

    if df1.columns[3] == "posterior0" and df2.columns[3] == "posterior1":
        df2.columns =  list(df2.columns[:3]) + ["posterior"+str(i)for i in range(len(df2.columns)-3)]
    elif df1.columns[3] == "posterior1" and df2.columns[3] == "posterior0":
        df1.columns =  list(df1.columns[:3]) + ["posterior"+str(i)for i in range(len(df1.columns)-3)]

        
    # to handle concat indexing
    if "_1" in df1.iloc[0, 0] or "_2" in df1.iloc[0, 0]:
        chrdf1 = list(df1.chr)
        for i in range(len(chrdf1)):
            if "_1" in chrdf1[i]:
                chrdf1[i] = chrdf1[i].replace("_1", "")
            elif "_2" in chrdf1[i]:
                chrdf1[i] = chrdf1[i].replace("_2", "")
        df1.chr = np.array(chrdf1)

    if "_1" in df2.iloc[0, 0] or "_2" in df2.iloc[0, 0]:
        chrdf2 = list(df2.chr)
        for i in range(len(chrdf2)):
            if "_2" in chrdf2[i]:
                chrdf2[i] = chrdf2[i].replace("_2", "")
            elif "_1" in chrdf2[i]:
                chrdf2[i] = chrdf2[i].replace("_1", "")
        df2.chr = np.array(chrdf2)

    df1 = df1.rename(columns={'chrom': 'chr'})
    df2 = df2.rename(columns={'chrom': 'chr'})

    intersect = pd.merge(
        df1, 
        df2, 
        how='inner', on=['chr', 'start', 'end'])

    df1 = [intersect.chr, intersect.start, intersect.end]
    df2 = [intersect.chr, intersect.start, intersect.end]

    for c in intersect.columns:
        if c[-1] == 'x':
            df1.append(intersect[c])
        elif c[-1] == 'y':
            df2.append(intersect[c])

    df1 = pd.concat(df1, axis=1)
    df1.columns = [c.replace("_x", "") for c in df1.columns]
    
    df2 = pd.concat(df2, axis=1)
    df2.columns = [c.replace("_y", "") for c in df2.columns]

    return df1, df2


def read_mnemonics(mnemon_file):
    df = pd.read_csv(mnemon_file, sep="\t")
    mnemon = []
    for i in range(len(df)):
        if int(df["old"][0]) == 0:
            mnemon.append(str(df["old"][i])+"_"+df["new"][i])
        elif int(df["old"][0]) == 1:
            mnemon.append(str(int(df["old"][i])-1)+"_"+df["new"][i])
    return mnemon


if __name__=="__main__":
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

    # Run segway replicates     MP
    partial_runs_i = partial(
        RunParse_segway_replicates, output_dir=segway_dir, random_seed=73)
    p_obj = mp.Pool(int(len(CellType_list)/2))
    p_obj.map(partial_runs_i, [download_dir+ct for ct in CellType_list])

    # parse_posteriors 
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
