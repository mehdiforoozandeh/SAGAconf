from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage.measurements import center_of_mass
from scipy.optimize import linear_sum_assignment
from scipy.spatial.kdtree import distance_matrix
from sklearn.cluster import KMeans
import plotly.express as px
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn import preprocessing
from sklearn.cluster import KMeans
from scipy.spatial.distance import pdist


def find_max_posteri(loci, num_labels=10):
    loci_posteri =  loci.loc[:, ['posterior'+str(i) for i in range(int(num_labels))]]
    for lpc in loci_posteri.columns:
        loci_posteri[lpc] = loci_posteri[lpc].astype('float')

    max_posteri = loci_posteri.idxmax(axis=1)
    return max_posteri

def confusion_matrix(loci_1, loci_2, num_labels):
    num_labels = int(num_labels)
    max_1posteri = find_max_posteri(loci_1, num_labels=num_labels)
    max_2posteri = find_max_posteri(loci_2, num_labels=num_labels)

    Confus_Mat = pd.DataFrame(
        np.zeros((int(num_labels), int(num_labels))),
        columns=['posterior'+str(i)for i in range(num_labels)],
        index=['posterior'+str(i)for i in range(num_labels)]
    )
    for i in range(len(loci_1)):
        Confus_Mat.loc[max_1posteri[i], max_2posteri[i]] += 1
    return Confus_Mat

def Hungarian_algorithm(matrix, conf_or_dis='conf'):

    if conf_or_dis == 'conf':
        confusion_matrix = np.array(matrix)
        best_assignments = linear_sum_assignment(confusion_matrix, maximize=True)

        print('Sum of optimal assignment sets / Sum of confusion matrix = {}/{}'.format(
            str(confusion_matrix[best_assignments[0],best_assignments[1]].sum()), 
            str(confusion_matrix.sum())))

        assignment_pairs = [(i, best_assignments[1][i]) for i in range(len(best_assignments[0]))]
        return assignment_pairs

    elif conf_or_dis == 'dist':
        distance_matrix = np.array(matrix)
        best_assignments = linear_sum_assignment(distance_matrix, maximize=False)

        assignment_pairs = [(i, best_assignments[1][i]) for i in range(len(best_assignments[0]))]
        return assignment_pairs

def connect_bipartite(loci_1, loci_2, assignment_matching):
    corrected_loci_2 = []
    corrected_loci_2.append(loci_2.chr)
    corrected_loci_2.append(loci_2.window_start)
    corrected_loci_2.append(loci_2.window_end)

    for i in range(len(assignment_matching)):
        corrected_loci_2.append(loci_2['posterior'+str(assignment_matching[i][1])])

    corrected_loci_2 = pd.concat(corrected_loci_2, axis=1)

    corrected_loci_2.columns = loci_1.columns
    return loci_1, corrected_loci_2
    

def read_length_dist_files(len_file_1, len_file_2):
    length_dist_1 = pd.read_csv(len_file_1, sep='\t')
    length_dist_1 = length_dist_1.loc[1:, ('mean.len', 'stdev.len')]
    length_dist_1 = length_dist_1.reset_index(drop=True)

    length_dist_2 = pd.read_csv(len_file_2, sep='\t')
    length_dist_2 = length_dist_2.loc[1:, ('mean.len', 'stdev.len')]
    length_dist_2 = length_dist_2.reset_index(drop=True)

    return length_dist_1, length_dist_2

def gmtk_params_files(gmtk_file_1, gmtk_file_2):
    gmtk_params_1 = pd.read_csv(gmtk_file_1)
    gmtk_params_1 = gmtk_params_1.drop("Unnamed: 0", axis=1)
    gmtk_params_2 = pd.read_csv(gmtk_file_2)
    gmtk_params_2 = gmtk_params_2.drop("Unnamed: 0", axis=1)

    return gmtk_params_1, gmtk_params_2

def create_clustering_data(length_dist_1, length_dist_2, gmtk_params_1, gmtk_params_2, normalize_len=True):
    if normalize_len:
        min_max_scaler = preprocessing.MinMaxScaler()

        tmp_len1 = length_dist_1.values #returns a numpy array
        x_scaled1 = min_max_scaler.fit_transform(tmp_len1)
        length_dist_1 = pd.DataFrame(x_scaled1, columns=length_dist_1.columns)

        tmp_len2 = length_dist_2.values #returns a numpy array
        x_scaled2 = min_max_scaler.fit_transform(tmp_len2)
        length_dist_2 = pd.DataFrame(x_scaled2, columns=length_dist_2.columns)

    rep1_data = pd.concat([gmtk_params_1, length_dist_1], axis=1)
    rep2_data = pd.concat([gmtk_params_2, length_dist_2], axis=1)

    return rep1_data, rep2_data

def PCA_plot(X, clusters,  PC=5, px=True):
    pca = PCA(n_components=PC)
    components = pca.fit_transform(X)
    labels = {
        str(i): f"PC {i+1} ({var:.1f}%)"
        for i, var in enumerate(pca.explained_variance_ratio_ * 100)
        }

    fig = px.scatter_matrix(
        components,
        labels=labels,
        color=clusters,
        dimensions=range(PC),
        )

    fig.update_traces(diagonal_visible=False)
    fig.update_traces(marker_size=12)
    fig.show()

def tsne_plot(X, clusters, n_components):
    tsne = TSNE(
        n_components=n_components, perplexity=10, random_state=0)

    z = tsne.fit_transform(X)
    sns.scatterplot(
        x=z[:,0], y=z[:,1], hue=clusters, s=100)

    plt.xlabel("TSNE_1")
    plt.ylabel("TSNE_2")
    
    plt.show()

def cluster(clustering_data, k):
    kmeans = KMeans(n_clusters=k, random_state=0)
    kmeans.fit(clustering_data)
    return kmeans

def update_labels_by_cluster(unclustered_loci, clustering_obj): 
    '''
    merges clusters and their corresponding posterior value in 
    each bin
    '''
    new_loci = unclustered_loci.loc[:,('chr', 'window_start', 'window_end')]
    cluster_labels = clustering_obj.labels_

    for i in range(len(cluster_labels)): 
        if "posterior_cluster_{}".format(cluster_labels[i]) not in new_loci.columns:
            new_loci["posterior_cluster_{}".format(cluster_labels[i])] = unclustered_loci["posterior{}".format(i)]
        else:
            new_loci["posterior_cluster_{}".format(cluster_labels[i])] = \
                new_loci["posterior_cluster_{}".format(cluster_labels[i])] + unclustered_loci["posterior{}".format(i)]
    
    return new_loci

def compute_pairwise_centroid_distance(clustering_obj_1, clustering_obj_2):
    '''
    connect centroids using min eucleadian distance.
    '''
    centroid_1 = np.array(clustering_obj_1.cluster_centers_)
    centroid_2 = np.array(clustering_obj_2.cluster_centers_)
    
    distance_matrix = np.zeros(
        (centroid_1.shape[0], centroid_2.shape[0]))

    for i in range(centroid_1.shape[0]):
        for j in range(centroid_2.shape[0]):
            
            dist = pdist(
                np.array([centroid_1[i,:], centroid_2[j,:]]))
            
            distance_matrix[i,j] = dist

    distance_matrix = pd.DataFrame(distance_matrix)
    return distance_matrix
