from cgitb import reset
from traceback import print_tb
import pandas as pd
import requests, json, itertools, ast


def search_encode(cell, download_dir, target_assembly="GRCh38"):
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

    print(frequent_pair, pair_num_tracks)

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

    # for each experiment, get info for all files in that experiment
    for e in range(len(tracks_navigation)):
        exp_url = "https://www.encodeproject.org/experiments/{}".format(tracks_navigation['accession'][e])
        
        exp_respond = requests.get(exp_url, headers=headers)
        exp_results = exp_respond.json()
        
        e_fileslist = list(exp_results['original_files'])
        e_files_navigation = []

        for ef in e_fileslist:
            efile_respond = requests.get("https://www.encodeproject.org{}".format(ef), headers=headers)
            efile_results = efile_respond.json()

            if efile_results['file_format'] == 'bigWig':
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
            'accession', 'file_format', 'output_type', 'experiment', 
            'bio_replicate_number', 'biosample', 'file_size', 'assembly', 
            'download_url', 'date_created', 'status'])
        
        # print(e_files_navigation)

        # define files to be downloaded
        
        to_download_list = {} #{"rep1_spv":None, "rep2_spv":None, "rep1_fcoc":None, "rep2_fcoc":None}

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

                
        # - download all tracks with that pair as replicates
            # - rep1 -> bigwig -> signal p-value (chromhmm)
            # - rep1 -> bigwig -> fold change over control (segway)
            # - rep2 -> bigwig -> signal p-value (chromhmm)
            # - rep2 -> bigwig -> fold change over control (segway)

        print(pd.DataFrame(to_download_list))

    

    # << RECORD METADATA FOR DOWNLOADED FILES >> and save in a csv
        # for each file, record:
            # biosample id , download_url, download_file_accession, 
            # signal p-value or fold change over control, genome assembly, date_added
    # create trackname_assays.txt

if __name__ == "__main__":
    search_encode('GM12878', 'files/')