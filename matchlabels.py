from _cluster_matching import *
from _pipeline import *

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

def main():
    len1, len2 = read_length_dist_files(
            'tests/rep1/length_distribution/segment_sizes.tab', 
            'tests/rep2/length_distribution/segment_sizes.tab')
    gmtk1, gmtk2 = gmtk_params_files(
            'tests/rep1/gmtk_parameters/gmtk_parameters.stats.csv', 
            'tests/rep2/gmtk_parameters/gmtk_parameters.stats.csv')

    rep1_data, rep2_data = create_clustering_data(
        len1, len2, gmtk1, gmtk2)

    c1 = cluster(rep1_data, 10)
    # tsne_plot(rep1_data, clusters=[str(i) for i in c1.labels_], n_components=2)
    # PCA_plot(rep1_data, PC=4, clusters=[str(i) for i in c1.labels_])

    c2 = cluster(rep2_data, 10)
    # tsne_plot(rep2_data, clusters=[str(i) for i in c2.labels_], n_components=2)
    # PCA_plot(rep2_data, PC=4, clusters=[str(i) for i in c2.labels_])

    dist_mat = compute_pairwise_centroid_distance(c1, c2)
    assignment_pairs = Hungarian_algorithm(dist_mat, conf_or_dis='dist')
    # print(assignment_pairs)
    # distance_matrix_heatmap(dist_mat)
    

    parsed_posterior_1 = pd.read_csv(
        'tests/rep1/parsed_posterior.csv').drop('Unnamed: 0', axis=1)
    parsed_posterior_2 = pd.read_csv(
        'tests/rep2/parsed_posterior.csv').drop('Unnamed: 0', axis=1)

    clustered_loci1 = update_labels_by_cluster(parsed_posterior_1, c1)
    clustered_loci2 = update_labels_by_cluster(parsed_posterior_2, c2)

    corrected_loci_1, corrected_loci_2 = \
        connect_bipartite(clustered_loci1, clustered_loci2, assignment_pairs, conf_or_dis='dist')
    
    cooc_mat = Cooccurrence_matrix(corrected_loci_1, corrected_loci_2)
    # coocurence_matrix_heatmap(cooc_mat)
    
    

if __name__=="__main__":
    main()