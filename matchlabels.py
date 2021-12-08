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
    sns.heatmap(
        cooc_mat.astype('int'), annot=True, fmt="d",
        linewidths=0.01, center=0, cbar=False,
        vmin=0,
        vmax=250)
        
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

def cluster_matching(rep_dir_1, rep_dir_2, n_clust, matching_strategy='distance_matrix'):

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


if __name__=="__main__":

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

