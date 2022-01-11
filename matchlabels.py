from math import exp
from os import link
from matplotlib.pyplot import plot
from numpy import NaN
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import seaborn as sns
from _cluster_matching import *
from _pipeline import *
from _visualize import visualize

def distance_matrix_heatmap(distance_matrix):
        sns.heatmap(
            distance_matrix, annot=True, fmt=".2f",
            linewidths=0.01, center=0, 
            cbar=True, vmin=0, vmax=2)
            
        # plt.colorbar()
        plt.title('Centroid Distance Heatmap')
        plt.xlabel('Replicate 1 Centroids')
        plt.ylabel("Replicate 2 Centroids")
        plt.show()

def coocurence_matrix_heatmap(cooc_mat):
    matmax = cooc_mat.max(axis=1).max(axis=0)
    matmin = cooc_mat.min(axis=1).min(axis=0)
    sns.heatmap(
        cooc_mat.astype('int'), annot=True, fmt="d",
        linewidths=0.01, center=0, cbar=False,
        vmin=matmin/2,
        vmax=matmax/2)
        
    plt.title('Cooccurrence_matrix Heatmap')
    plt.xlabel('Replicate 1 clusters')
    plt.ylabel("Replicate 2 clusters")
    plt.show()

def plain_seg_matching(rep_dir_1, rep_dir_2, matching_strategy='conf_matrix'):

    if matching_strategy == 'conf_matrix':
        parsed_posterior_1 = pd.read_csv(
            '{}/parsed_posterior.csv'.format(rep_dir_1)).drop('Unnamed: 0', axis=1)
        parsed_posterior_2 = pd.read_csv(
            '{}/parsed_posterior.csv'.format(rep_dir_2)).drop('Unnamed: 0', axis=1)
        
        conf_mat = confusion_matrix(parsed_posterior_1, parsed_posterior_2, parsed_posterior_1.shape[1]-3)
        assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')
        match_eval = match_evaluation(conf_mat, assignment_pairs)

        print("label_mapping:\t", assignment_pairs)
        corrected_loci_1, corrected_loci_2 = connect_bipartite(
            parsed_posterior_1, parsed_posterior_2, assignment_pairs)
        print(match_eval)

        vis = visualize(corrected_loci_1, corrected_loci_2, corrected_loci_1.shape[1]-3)
        # vis.confusion_matrix_heatmap()
    
    elif matching_strategy == 'distance_matrix':
        # len1, len2 = read_length_dist_files(
        #     '{}/length_distribution/segment_sizes.tab'.format(rep_dir_1), 
        #     '{}/length_distribution/segment_sizes.tab'.format(rep_dir_2))
        gmtk1, gmtk2 = gmtk_params_files(
            '{}/gmtk_parameters/gmtk_parameters.stats.csv'.format(rep_dir_1), 
            '{}/gmtk_parameters/gmtk_parameters.stats.csv'.format(rep_dir_2))
        
        rep1_data, rep2_data = gmtk1, gmtk2 #curate_dataset(len1, len2, gmtk1, gmtk2) 

        dist_mat = compute_pairwise_centroid_distance(np.array(rep1_data), np.array(rep2_data))
        assignment_pairs = Hungarian_algorithm(dist_mat, conf_or_dis='dist')
        # print(assignment_pairs)
        # distance_matrix_heatmap(dist_mat)

        parsed_posterior_1 = pd.read_csv(
            '{}/parsed_posterior.csv'.format(rep_dir_1)).drop('Unnamed: 0', axis=1)
        parsed_posterior_2 = pd.read_csv(
            '{}/parsed_posterior.csv'.format(rep_dir_2)).drop('Unnamed: 0', axis=1)

        corrected_loci_1, corrected_loci_2 = \
            connect_bipartite(parsed_posterior_1, parsed_posterior_2, assignment_pairs)
        
        cooc_mat = Cooccurrence_matrix(corrected_loci_1, corrected_loci_2)
        eval_results = match_evaluation(cooc_mat, [(i, i) for i in range(cooc_mat.shape[0])])

        # coocurence_matrix_heatmap(cooc_mat)
        print(eval_results)
    
    return corrected_loci_1, corrected_loci_2

def cluster_matching(rep_dir_1, rep_dir_2, n_clust, matching_strategy='distance_matrix', plot_dendogram=False):

    if matching_strategy == 'distance_matrix':
        # len1, len2 = read_length_dist_files(
        #         '{}/length_distribution/segment_sizes.tab'.format(rep_dir_1), 
        #         '{}/length_distribution/segment_sizes.tab'.format(rep_dir_2))

        gmtk1, gmtk2 = gmtk_params_files(
                '{}/gmtk_parameters/gmtk_parameters.stats.csv'.format(rep_dir_1), 
                '{}/gmtk_parameters/gmtk_parameters.stats.csv'.format(rep_dir_2))

        rep1_data, rep2_data = gmtk1, gmtk2 #curate_dataset(len1, len2, gmtk1, gmtk2) 

        clstrer_1 = Clusterer(rep1_data, n_clusters=int(n_clust))
        c1= clstrer_1.fit_hierarchical(metric='euclidean', linkage='ward')
        # tsne_plot(rep1_data, clusters=[str(i) for i in c1.labels_], n_components=2)
        # PCA_plot(rep1_data, PC=2, clusters=[str(i) for i in c1.labels_], px=False)

        clstrer_2 = Clusterer(rep2_data, n_clusters=int(n_clust))
        c2 = clstrer_2.fit_hierarchical(metric='euclidean', linkage='ward')
        # tsne_plot(rep2_data, clusters=[str(i) for i in c2.labels_], n_components=2)
        # PCA_plot(rep2_data, PC=2, clusters=[str(i) for i in c2.labels_], px=False)

        if plot_dendogram:
            clstrer_1.dendrogram(truncate_mode=None, distance_sort='ascending', color_threshold=1)
            clstrer_2.dendrogram(truncate_mode=None, distance_sort='ascending', color_threshold=1)
        
        dist_mat = compute_pairwise_centroid_distance(c1.cluster_centers_, c2.cluster_centers_)
        assignment_pairs = Hungarian_algorithm(dist_mat, conf_or_dis='dist')
        # print(assignment_pairs)
        # distance_matrix_heatmap(dist_mat)
        
        parsed_posterior_1 = pd.read_csv(
            '{}/parsed_posterior.csv'.format(rep_dir_1)).drop('Unnamed: 0', axis=1)
        parsed_posterior_2 = pd.read_csv(
            '{}/parsed_posterior.csv'.format(rep_dir_2)).drop('Unnamed: 0', axis=1)

        clustered_loci1 = update_labels_by_cluster(parsed_posterior_1, c1)
        clustered_loci2 = update_labels_by_cluster(parsed_posterior_2, c2)

        corrected_loci_1, corrected_loci_2 = \
            connect_bipartite(clustered_loci1, clustered_loci2, assignment_pairs)
        
        cooc_mat = Cooccurrence_matrix(corrected_loci_1, corrected_loci_2)
        eval_results = match_evaluation(cooc_mat, [(i, i) for i in range(cooc_mat.shape[0])])

        print(eval_results)
        # coocurence_matrix_heatmap(cooc_mat)

    elif matching_strategy == 'conf_matrix':
        # len1, len2 = read_length_dist_files(
        #         '{}/length_distribution/segment_sizes.tab'.format(rep_dir_1), 
        #         '{}/length_distribution/segment_sizes.tab'.format(rep_dir_2))
        gmtk1, gmtk2 = gmtk_params_files(
                '{}/gmtk_parameters/gmtk_parameters.stats.csv'.format(rep_dir_1), 
                '{}/gmtk_parameters/gmtk_parameters.stats.csv'.format(rep_dir_2))

        rep1_data, rep2_data = gmtk1, gmtk2 #curate_dataset(len1, len2, gmtk1, gmtk2) 

        clstrer_1 = Clusterer(rep1_data, n_clusters=int(n_clust))
        c1= clstrer_1.fit_hierarchical(metric='euclidean', linkage='ward')
        # tsne_plot(rep1_data, clusters=[str(i) for i in c1.labels_], n_components=2)
        # PCA_plot(rep1_data, PC=2, clusters=[str(i) for i in c1.labels_], px=False)

        clstrer_2 = Clusterer(rep2_data, n_clusters=int(n_clust))
        c2 = clstrer_2.fit_hierarchical(metric='euclidean', linkage='ward')
        # tsne_plot(rep2_data, clusters=[str(i) for i in c2.labels_], n_components=2)
        # PCA_plot(rep2_data, PC=2, clusters=[str(i) for i in c2.labels_], px=False)

        if plot_dendogram:
            clstrer_1.dendrogram(truncate_mode=None, distance_sort='ascending', color_threshold=1)
            clstrer_2.dendrogram(truncate_mode=None, distance_sort='ascending', color_threshold=1)
        
        parsed_posterior_1 = pd.read_csv(
            '{}/parsed_posterior.csv'.format(rep_dir_1)).drop('Unnamed: 0', axis=1)
        parsed_posterior_2 = pd.read_csv(
            '{}/parsed_posterior.csv'.format(rep_dir_2)).drop('Unnamed: 0', axis=1)

        clustered_loci1 = update_labels_by_cluster(parsed_posterior_1, c1)
        clustered_loci2 = update_labels_by_cluster(parsed_posterior_2, c2)
        
        conf_mat = confusion_matrix(clustered_loci1, clustered_loci2, clustered_loci1.shape[1]-3)
        assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')
        print("label_mapping:\t", assignment_pairs)

        corrected_loci_1, corrected_loci_2 = \
            connect_bipartite(clustered_loci1, clustered_loci2, assignment_pairs)
        
        cooc_mat = Cooccurrence_matrix(corrected_loci_1, corrected_loci_2)
        eval_results = match_evaluation(cooc_mat, [(i, i) for i in range(cooc_mat.shape[0])])

        print(eval_results)
        # coocurence_matrix_heatmap(cooc_mat)

    return corrected_loci_1, corrected_loci_2


def order_based_clustering(rep_dir_1, rep_dir_2, OE_transform=True, log_transform=True):
    '''
    read posterior
    create cooccurence matrix based on overlap (OE transformed ?)
    match labels using hungarian
    update cooc-matrix (OE transformed)

    use cooc-matrix as similarity matrix (or inverse function to convert it to distance matrix)

    agglomerative clustering
    check different orders of clustering and record matching rate at different levels
    '''

    parsed_posterior_1 = pd.read_csv(
            '{}/parsed_posterior.csv'.format(rep_dir_1)).drop('Unnamed: 0', axis=1)
    parsed_posterior_2 = pd.read_csv(
            '{}/parsed_posterior.csv'.format(rep_dir_2)).drop('Unnamed: 0', axis=1)

    num_labels = parsed_posterior_1.shape[1]-3

    conf_mat = confusion_matrix(
        parsed_posterior_1, parsed_posterior_2, num_labels, 
        OE_transform=OE_transform, log_transform=log_transform)
    
    assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')
    print("label_mapping:\t", assignment_pairs)

    corrected_loci_1, corrected_loci_2 = \
        connect_bipartite(parsed_posterior_1, parsed_posterior_2, assignment_pairs)


    '''
    coocurrence_matrix ->  records the original cooc_mat which is 
    going to be updated in the order of cluster merging'''

    
    coocurrence_matrix = confusion_matrix(
                    corrected_loci_1, corrected_loci_2, num_labels, 
                    OE_transform=False, log_transform=False)
                    
    print('Match Rate:\n {}'.format(
            match_evaluation(coocurrence_matrix, 
            [(i, i) for i in range(num_labels)])
            ))

    conf_mat = confusion_matrix(
        corrected_loci_1, corrected_loci_2, num_labels, 
        OE_transform=OE_transform, log_transform=log_transform)    

    max_similarity = conf_mat.max(axis=1).max(axis=0)
    min_similarity = conf_mat.min(axis=1).min(axis=0)

    # min-max scale the similarity matrix
    conf_mat = (conf_mat - min_similarity) / (max_similarity - min_similarity)

    # calculate pairwise correlation, so that the matrix is symmetrical
    conf_mat = conf_mat.corr()

    # convert similarity (corr) to distance matrix
    distance_matrix = 1 - conf_mat

    # generate linkage matrix to be used for clustering and creating dendograms
    linkage = hc.linkage(sp.distance.squareform(distance_matrix), method='average')

    # use sns to visualize heatmap and dendogram
    c_grid = sns.clustermap(distance_matrix, row_linkage=linkage, col_linkage=linkage, annot=True)

    print(linkage)
    plt.show()

    # merge clusters recursively based on the order of branching 
    orders_of_branching = np.unique(linkage[:,3])

    # record all the intermediate steps and return for further analysis
    record = {
        'new_clust_address':[], 'linkage':linkage,
        'order_of_branching':[], 'num_subclusters':[],
        'loci_1_history':[], 'loci_2_history':[], 
        'overlap_history':[], 'match_rate':[]}

    for o in orders_of_branching:
        for m in range(len(linkage)):
            if linkage[m, 3] == o:
                to_be_merged = [
                    "posterior{}".format(int(linkage[m, 0])), 
                    "posterior{}".format(int(linkage[m, 1]))]

                # print('posterior{}'.format(num_labels + m), to_be_merged, int(np.where(orders_of_branching == o)[0])+1 )

                corrected_loci_1['posterior{}'.format(num_labels + m)] = \
                    corrected_loci_1[to_be_merged[0]] + corrected_loci_1[to_be_merged[1]]

                corrected_loci_1 = corrected_loci_1.drop(to_be_merged, axis=1)

                corrected_loci_2['posterior{}'.format(num_labels + m)] = \
                    corrected_loci_2[to_be_merged[0]] + corrected_loci_2[to_be_merged[1]]

                corrected_loci_2 = corrected_loci_2.drop(to_be_merged, axis=1)

                max_1posteri = find_max_posteri(corrected_loci_1)
                max_2posteri = find_max_posteri(corrected_loci_2)

                observed_overlap = pd.DataFrame(
                    np.zeros((int(corrected_loci_2.shape[1]-3), int(corrected_loci_2.shape[1]-3))),
                    columns=list(corrected_loci_1.columns)[3:],
                    index=list(corrected_loci_1.columns)[3:])

                for i in range(len(corrected_loci_1)):
                    try:
                        observed_overlap.loc[max_1posteri[i], max_2posteri[i]] += 1
                    except:
                        pass
                
                MR = match_evaluation(observed_overlap, [(i, i) for i in range(corrected_loci_1.shape[1]-3)])['all']

                record['new_clust_address'].append(['posterior{}'.format(num_labels + m), to_be_merged])
                record['order_of_branching'].append(int(np.where(orders_of_branching == o)[0])+1)
                record['num_subclusters'].append(o)
                record['loci_1_history'].append(corrected_loci_1)
                record['loci_2_history'].append(corrected_loci_2)
                record['overlap_history'].append(observed_overlap)
                record['match_rate'].append(MR)

                print('Match Rate:\n {}'.format(MR))

    return record    

if __name__=="__main__":
    order_based_clustering('tests/rep1', 'tests/rep2')
    exit()

    pilot_regions_file = 'segway_runs/encodePilotRegions.hg19.bed'

    for n_labels in range(8,17):
        replicate_1_dir = 'segway_runs/num_labels_analysis/rep1_{}labels_res1000'.format(n_labels)
        replicate_2_dir = 'segway_runs/num_labels_analysis/rep2_{}labels_res1000'.format(n_labels)

        print('plain seg => conf_mat // k = {}'.format(n_labels))
        plain_seg_matching(\
            replicate_1_dir, replicate_2_dir, matching_strategy='conf_matrix')

        print('plain seg => dist_mat // k = {}'.format(n_labels))
        plain_seg_matching(\
            replicate_1_dir, replicate_2_dir, matching_strategy='distance_matrix')
        
        replicate_1_dir = 'segway_runs/num_labels_analysis/rep1_{}labels_res1000'.format(16)
        replicate_2_dir = 'segway_runs/num_labels_analysis/rep2_{}labels_res1000'.format(16)
        
        print('post clst => conf_mat // k = {}'.format(n_labels))
        cluster_matching(\
            replicate_1_dir, replicate_2_dir, n_clust=n_labels ,matching_strategy='conf_matrix')
        
        print('post clst => dist_mat // k = {}'.format(n_labels))
        cluster_matching(\
            replicate_1_dir, replicate_2_dir, n_clust=n_labels ,matching_strategy='distance_matrix')

        print('\n\n\n')

