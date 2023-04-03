import pandas as pd
import numpy as np
import requests, os, itertools, ast


def search_encode(cell, download_dir, target_assembly="GRCh38", check_availability=False):
    # Force return from the server in JSON format
    headers = {'accept': 'application/json'}

    base_search_url = '''https://www.encodeproject.org/search/?type=Experiment&assay_title=Histone+ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&replication_type=isogenic&biosample_ontology.term_name=''' + str(cell)
    
    # GET the object
    response = requests.get(base_search_url, headers=headers)

    # Extract the JSON response as a Python dictionary
    search_results = response.json()
    
    # - iterate over different search results
    tracks_navigation = []
    for sr in search_results['@graph']: 
    
        sr_instance = []
        sr_instance.append(sr['accession'])
        sr_instance.append(sr['target']['label'])

        # - record biosample ID for both replicates
        replicates = []
        for r in sr['replicates']:
            replicates.append(r['library']['biosample']['accession'])
        
        sr_instance.append(tuple(set(replicates)))
        tracks_navigation.append(sr_instance)

    tracks_navigation = pd.DataFrame(tracks_navigation, columns=['accession', 'assay', 'replicates'])
    print(tracks_navigation)

    # - choose the most frequent pair of biosamples. 
    #  (should create all possible combinations if an experiment has more than 2 replicates)

    replicate_pairs = {}
    for i in range(len(tracks_navigation)):
        pairs = [str(sorted(x)) for x in itertools.combinations(tracks_navigation['replicates'][i], 2)]

        for pair in pairs:
            if pair in replicate_pairs.keys():
                replicate_pairs[pair] += 1
            else:
                replicate_pairs[pair] = 1
    
    print(replicate_pairs)
    
    frequent_pair = ast.literal_eval(max(replicate_pairs, key=replicate_pairs.get))
    pair_num_tracks = max(list(replicate_pairs.values()))

    print('most frequent pair: {}, with {} tracks'.format(frequent_pair, pair_num_tracks))

    if check_availability:
        return pair_num_tracks

    # here, we establist rep1 and rep2 based on the order by which they appear in list frequent_pair.
    replicate_glossary = {'rep1':frequent_pair[0], 'rep2':frequent_pair[1]}

    to_drop = []
    for i in range(len(tracks_navigation)):
        statement = frequent_pair[0] in tracks_navigation['replicates'][i] and frequent_pair[1] in tracks_navigation['replicates'][i]
        if statement == False:
            to_drop.append(i)
    tracks_navigation = tracks_navigation.drop(to_drop, axis=0)
    tracks_navigation = tracks_navigation.reset_index(drop=True)
    print(tracks_navigation)

    if not os.path.exists(download_dir):
        os.mkdir(download_dir)
        os.mkdir(download_dir + cell)
    
    if not os.path.exists(download_dir + cell):
        os.mkdir(download_dir + cell)

    with open(download_dir + cell + '/replicate_number.txt' ,'w') as rep_file:
            rep_file.write("rep1:{}\nrep2:{}".format(
                replicate_glossary['rep1'], replicate_glossary['rep2']))

    # for each experiment, get info for all files in that experiment
    for e in range(len(tracks_navigation)):
        print('downloading celltype {} - assay {}'.format(cell, tracks_navigation['assay'][e]))
        exp_url = "https://www.encodeproject.org/experiments/{}".format(tracks_navigation['accession'][e])
        
        exp_respond = requests.get(exp_url, headers=headers)
        exp_results = exp_respond.json()
        
        e_fileslist = list(exp_results['original_files'])
        e_files_navigation = []

        for ef in e_fileslist:
            efile_respond = requests.get("https://www.encodeproject.org{}".format(ef), headers=headers)
            efile_results = efile_respond.json()

            if efile_results['file_format'] == 'bigWig' or efile_results['file_format'] == 'bam':
                try: #ignore files without sufficient info or metadata

                    if ',' not in str(efile_results['origin_batches']):
                        if efile_results['status'] == "released": 
                            #ignore old and depricated versions

                            e_file_biosample = str(efile_results['origin_batches'])
                            e_file_biosample = e_file_biosample.replace('/', '')
                            e_file_biosample = e_file_biosample.replace('biosamples','')
                            # ignore files that contain both replicates 

                            e_files_navigation.append(
                                [
                                    tracks_navigation['assay'][e],
                                    efile_results['accession'], efile_results['file_format'], 
                                    efile_results['output_type'], efile_results['dataset'], 
                                    efile_results['biological_replicates'], e_file_biosample[2:-2],
                                    efile_results['file_size'], efile_results['assembly'], 
                                    "https://www.encodeproject.org{}".format(efile_results['href']), 
                                    efile_results['date_created'], efile_results['status'] 
                                ]
                            )
                except:
                    pass

        e_files_navigation = pd.DataFrame(e_files_navigation, columns=[
            'assay', 'accession', 'file_format', 'output_type', 'experiment', 
            'bio_replicate_number', 'biosample', 'file_size', 'assembly', 
            'download_url', 'date_created', 'status'])

        # define files to be downloaded
        to_download_list = {} #{"rep1_alig":None, "rep2_alig":None, "rep1_fcoc":None, "rep2_fcoc":None}

        if not os.path.exists(download_dir + cell  +"/{}".format(tracks_navigation['assay'][e])):
            os.mkdir(download_dir + cell +"/{}".format(tracks_navigation['assay'][e]))
        
        for f in range(len(e_files_navigation)):
            if e_files_navigation['assembly'][f] == target_assembly:

                if e_files_navigation['output_type'][f] == "fold change over control":
                    if e_files_navigation['biosample'][f] == replicate_glossary['rep1']:

                        if 'rep1_fcoc' not in to_download_list.keys():
                            to_download_list['rep1_fcoc'] = e_files_navigation.iloc[f, :]

                        elif int(str(
                            e_files_navigation['date_created'][f])[:7].replace('-', '')) > int(str(
                                to_download_list['rep1_fcoc']['date_created'])[:7].replace('-', '')):
                            # substitute by a more recent version
                            to_download_list['rep1_fcoc'] = e_files_navigation.iloc[f, :]

                    elif e_files_navigation['biosample'][f] == replicate_glossary['rep2']:

                        if 'rep2_fcoc' not in to_download_list.keys():
                            to_download_list['rep2_fcoc'] = e_files_navigation.iloc[f, :]

                        elif int(str(
                            e_files_navigation['date_created'][f])[:7].replace('-', '')) > int(str(
                                to_download_list['rep2_fcoc']['date_created'])[:7].replace('-', '')):
                            # substitute by a more recent version
                            to_download_list['rep2_fcoc'] = e_files_navigation.iloc[f, :]

                elif e_files_navigation['output_type'][f] == "alignments":
                    if e_files_navigation['biosample'][f] == replicate_glossary['rep1']:
                        
                        if 'rep1_alig' not in to_download_list.keys():
                            to_download_list['rep1_alig'] = e_files_navigation.iloc[f, :]

                        elif int(str(
                            e_files_navigation['date_created'][f])[:7].replace('-', '')) > int(str(
                                to_download_list['rep1_alig']['date_created'])[:7].replace('-', '')):
                            # substitute by a more recent version
                            to_download_list['rep1_alig'] = e_files_navigation.iloc[f, :]
                    
                    elif e_files_navigation['biosample'][f] == replicate_glossary['rep2']:
                        
                        if 'rep2_alig' not in to_download_list.keys():
                            to_download_list['rep2_alig'] = e_files_navigation.iloc[f, :]

                        elif int(str(
                            e_files_navigation['date_created'][f])[:7].replace('-', '')) > int(str(
                                to_download_list['rep2_alig']['date_created'])[:7].replace('-', '')):
                            # substitute by a more recent version
                            to_download_list['rep2_alig'] = e_files_navigation.iloc[f, :]
                
                elif e_files_navigation['output_type'][f] == "signal p-value":
                    if e_files_navigation['biosample'][f] == replicate_glossary['rep1']:
                        
                        if 'rep1_spv' not in to_download_list.keys():
                            to_download_list['rep1_spv'] = e_files_navigation.iloc[f, :]

                        elif int(str(
                            e_files_navigation['date_created'][f])[:7].replace('-', '')) > int(str(
                                to_download_list['rep1_spv']['date_created'])[:7].replace('-', '')):
                            # substitute by a more recent version
                            to_download_list['rep1_spv'] = e_files_navigation.iloc[f, :]
                    
                    elif e_files_navigation['biosample'][f] == replicate_glossary['rep2']:
                        
                        if 'rep2_spv' not in to_download_list.keys():
                            to_download_list['rep2_spv'] = e_files_navigation.iloc[f, :]

                        elif int(str(
                            e_files_navigation['date_created'][f])[:7].replace('-', '')) > int(str(
                                to_download_list['rep2_spv']['date_created'])[:7].replace('-', '')):
                            # substitute by a more recent version
                            to_download_list['rep2_spv'] = e_files_navigation.iloc[f, :]

                
        to_download_list = pd.DataFrame(to_download_list)
        
        # - Navigate all tracks with that pair as replicates
            # - rep1 -> bigwig -> alignments (chromhmm)
            # - rep1 -> bigwig -> fold change over control (segway)
            # - rep2 -> bigwig -> alignments (chromhmm)
            # - rep2 -> bigwig -> fold change over control (segway)

        # << RECORD METADATA FOR DOWNLOADED FILES >> and save in a csv
        # for each file, record:
            # biosample id , download_url, download_file_accession, 
            # signal p-value or fold change over control, genome assembly, date_added

        to_download_list.to_csv(download_dir + cell + "/{}".format(tracks_navigation['assay'][e]+'/track_files_metadata.csv'))

        # download navigated files
        for c in to_download_list.columns:
            if "fcoc" in c:
                save_dir_name = download_dir + cell +"/{}".format(tracks_navigation['assay'][e])+ '/' + to_download_list.loc['accession', c] + ".bigWig"
            elif "alig" in c:
                save_dir_name = download_dir + cell +"/{}".format(tracks_navigation['assay'][e])+ '/' + to_download_list.loc['accession', c] + ".bam"
            elif "spv" in c:
                save_dir_name = download_dir + cell +"/{}".format(tracks_navigation['assay'][e])+ '/' + to_download_list.loc['accession', c] + ".bigWig"
            
            download_link = to_download_list.loc['download_url', c]
            download_response = requests.get(download_link, allow_redirects=True)
            open(save_dir_name, 'wb').write(download_response.content)

def get_data_from_csv(csvfile ="data_summary.csv", downloaddir="files/", bw2bg=True):
    summary = pd.read_csv(csvfile).drop("Unnamed: 0", axis=1)
    summary.Celltype = summary.Celltype.astype(str)

    if os.path.exists(downloaddir) == False:
        os.mkdir(downloaddir)

    listofct = list(set(summary.Celltype.astype(str)))
    for i in range(summary.shape[0]):
        ct = summary["Celltype"][i]
        assay = summary["experiment"][i]
        fileslist = summary["files"][i].split(", ")
        print(ct, assay)
        
        if os.path.exists("""{}/{}/""".format(downloaddir, ct)) == False:
            os.mkdir("""{}/{}/""".format(downloaddir, ct))

        if os.path.exists("""{}/{}/{}/""".format(downloaddir, ct, assay)) == False:
            os.mkdir("""{}/{}/{}/""".format(downloaddir, ct, assay))

        #bam##########################################################################################
        save_dir_name = """{}/{}/{}/{}""".format(downloaddir, ct, assay, fileslist[0]+".bam")
        download_link = """https://www.encodeproject.org/files/{}/@@download/{}.{}""".format(
            fileslist[0], fileslist[0], ".bam"
        )
        # download_link = """https://www.encodeproject.org/files/{}""".format(
        #     fileslist[0]
        # )
        download_response = requests.get(download_link, allow_redirects=True)
        open(save_dir_name, 'wb').write(download_response.content)

        #bigWig#######################################################################################
        save_dir_name = """{}/{}/{}/{}""".format(downloaddir, ct, assay, fileslist[1]+".bigWig")
        download_link = """https://www.encodeproject.org/files/{}/@@download/{}.{}""".format(
            fileslist[1], fileslist[1], ".bigWig"
        )
        # download_link = """https://www.encodeproject.org/files/{}""".format(
        #     fileslist[1]
        # )
        download_response = requests.get(download_link, allow_redirects=True)
        open(save_dir_name, 'wb').write(download_response.content)
        # if 


def create_trackname_assay_file(download_dir):
    tracknames = []
    CellType_list = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+ct)]

    for ct in CellType_list:
        tracks_list = [tr for tr in os.listdir(download_dir+ct) if os.path.isdir(download_dir+ct+'/'+tr)]
        for tr in tracks_list:
            tfmd = pd.read_csv(download_dir+ct+'/'+tr+'/track_files_metadata.csv')
            tfmd.index = list(tfmd['Unnamed: 0'])
            tfmd = tfmd.drop('Unnamed: 0', axis=1)
            try:
                for c in range(len(tfmd.columns)):
                    tracknames.append([tr, tfmd.loc['accession', tfmd.columns[c]]])
            except:
                pass

    with open(download_dir + '/trackname_assay.txt', 'w') as tna:
        for i in range(len(tracknames)):
            tna.write('{}\t{}\n'.format(tracknames[i][1], tracknames[i][0]))

if __name__ == "__main__":
    get_data_from_csv(csvfile ="data_summary.csv")

    exit()
    CellType_list = np.array(
        ['K562', 'MCF-7', 'GM12878', 'HeLa-S3', 'CD14-positive monocyte'])
    for ct in CellType_list:
        search_encode(ct, "files/", target_assembly="GRCh38", check_availability=False)

    # clean up potential space characters in directory names to prevent later issues
    for ct in CellType_list:
        if " " in ct:
            os.system("mv {} {}".format(
                ct.replace(' ', '\ '), ct.replace(" ", "_")
            ))
 
    exit()


    available_data_count = {}
    list_of_ENCODE_celltypes = ['H1', 'IMR-90', 'H9', 'spleen', 'transverse colon', 'trophoblast cell', 'K562', 'MCF-7', 'SK-N-SH', 'motor neuron', 'HCT116', 'neuronal stem cell', 'endothelial cell of umbilical vein', 'keratinocyte', 'GM12878', 'HeLa-S3', 'HepG2', 'gastrocnemius medialis', 'mesendoderm', 'MM.1S', 'OCI-LY1', 'PC-9', 'mammary epithelial cell', 'DND-41', 'DOHH2', 'GM23248', 'Karpas-422', 'Loucy', 'PC-3', 'SU-DHL-6', 'hepatocyte', 'neural cell', 'neural progenitor cell', 'skeletal muscle myoblast', 'smooth muscle cell', 'CD14-positive monocyte', 'GM23338', 'NCI-H929', 'OCI-LY3', 'astrocyte', 'mesenchymal stem cell', 'myotube', 'KMS-11', 'OCI-LY7', 'adrenal gland', 'fibroblast of dermis', 'prostate gland', 'thyroid gland', 'A549', 'HUES64', 'Panc1', 'SK-N-MC', 'body of pancreas', 'cardiac muscle cell', 'endodermal cell', 'skeletal muscle satellite cell', 'A673', 'HAP-1', 'HUES48', 'HUES6', 'NT2/D1', 'RWPE2', 'SJCRH30', 'SJSA1', 'WERI-Rb-1', 'iPS DF 19.11', 'iPS-20b', 'left ventricle myocardium inferior', 'thoracic aorta', 'B cell', 'BE2C', 'Caco-2', 'HEK293', 'UCSF-4', 'iPS DF 6.9', 'ACC112', 'ES-I3', 'KOPT-K1', 'MG63', 'brain microvascular endothelial cell', 'ectodermal cell', 'esophagus muscularis mucosa', 'mesodermal cell', 'mid-neurogenesis radial glial cells', 'mononuclear cell', 'neuroepithelial stem cell', 'radial glial cell', 'BJ', 'GM06990', 'H7', 'HL-60', 'bronchial epithelial cell', 'foreskin fibroblast', 'gastroesophageal sphincter', 'kidney epithelial cell', 'tibial artery', 'AG04450', 'esophagus squamous epithelium', 'fibroblast of lung', 'large intestine', 'right atrium auricular region', 'right lobe of liver', 'uterus', 'vagina', '22Rv1', 'AG04449', 'AG09309', 'AG09319', 'AG10803', 'C4-2B', 'GM08714', 'GM12864', 'GM12865', 'HFF-Myc', 'Jurkat, Clone E6-1', 'LNCaP clone FGC', "Peyer's patch", 'RWPE1', 'VCaP', 'WI38', 'astrocyte of the cerebellum', 'astrocyte of the spinal cord', 'breast epithelium', 'cardiac fibroblast', 'choroid plexus epithelial cell', 'common myeloid progenitor, CD34-positive', 'coronary artery', 'epithelial cell of esophagus', 'epithelial cell of prostate', 'epithelial cell of proximal tubule', 'fibroblast of mammary gland', 'fibroblast of pulmonary artery', 'fibroblast of the aortic adventitia', 'fibroblast of villous mesenchyme', 'foreskin melanocyte', 'heart left ventricle', 'iPS-18c', 'lower leg skin', 'muscle of leg', 'neuron', 'osteoblast', 'pancreas', 'retinal pigment epithelial cell', 'skeletal muscle cell', 'stomach', 'substantia nigra']

    for c in list_of_ENCODE_celltypes:
        try:
            available_data_count[c] = search_encode(c, 'files/', check_availability=True)
        except:
            pass

    available_data_count = {k: v for k, v in sorted(available_data_count.items(), key=lambda item: item[1])}
    print(available_data_count)

    # SELECTED_CELLTYPES = ['K562', 'MCF-7', 'GM12878', 'HeLa-S3', 'CD14-positive monocyte']

    # selection filters
    # 1- tier 1, 2, 2.5 (based on https://www.genome.gov/encode-project-common-cell-types)
    # 2- available replicated tracks => 10