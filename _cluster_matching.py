from math import log
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from plotly import plot
from scipy.optimize import linear_sum_assignment
from scipy.sparse import coo
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import KMeans, AgglomerativeClustering
import plotly.express as px
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn import preprocessing
from sklearn.cluster import KMeans
from scipy.spatial.distance import pdist


def find_max_posteri(loci):
    loci_posteri =  loci.iloc[:, 3:]
    for lpc in loci_posteri.columns:
        loci_posteri[lpc] = loci_posteri[lpc].astype('float')

    max_posteri = loci_posteri.idxmax(axis=1)
    return max_posteri

def confusion_matrix(loci_1, loci_2, num_labels, OE_transform=True, log_transform=True):
    num_labels = int(num_labels)
    max_1posteri = find_max_posteri(loci_1)
    max_2posteri = find_max_posteri(loci_2)

    observed_overlap = pd.DataFrame(
        np.zeros((int(num_labels), int(num_labels))),
        columns=['posterior{}'.format(i) for i in range(num_labels)],
        index=['posterior{}'.format(i) for i in range(num_labels)]
    )

    for i in range(len(loci_1)):
        try:
            observed_overlap.loc[max_1posteri[i], max_2posteri[i]] += 1
        except:
            pass
    
    if OE_transform:
        epsilon = 1e-3
        expected_overlap = np.zeros((num_labels, num_labels))
        r1_label_bin_count = {}
        r2_label_bin_count = {}
        
        '''
        expected overlap of label x = fraction of genome with label x 
        in rep1 * fraction of genome with label x in rep2
        '''

        for i in range(len(loci_1)):
            if max_1posteri[i] in r1_label_bin_count.keys():
                r1_label_bin_count[max_1posteri[i]] += 1

            else:
                r1_label_bin_count[max_1posteri[i]] = 1


        for i in range(len(loci_2)):
            if max_2posteri[i] in r2_label_bin_count.keys():
                r2_label_bin_count[max_2posteri[i]] += 1

            else:
                r2_label_bin_count[max_2posteri[i]] = 1

        for i in range(num_labels):
            for j in range(num_labels):

                # expected_overlap[i, j] = (r1_label_bin_count['posterior'+str(i)]/corrected_loci_1.shape[0])\
                #      * (r2_label_bin_count['posterior'+str(j)]/corrected_loci_2.shape[0]) * corrected_loci_1.shape[0]

                expected_overlap[i, j] = (r1_label_bin_count['posterior'+str(i)] * r2_label_bin_count['posterior'+str(j)])\
                    / loci_1.shape[0]

        expected_overlap = pd.DataFrame(
            expected_overlap, 
            columns=['posterior'+str(i)for i in range(num_labels)], 
            index=['posterior'+str(i)for i in range(num_labels)])
            
        oe_overlap = (observed_overlap + epsilon) / (expected_overlap + epsilon)
        if log_transform:
            oe_overlap = np.log(oe_overlap)

        return oe_overlap

    else:
        return observed_overlap

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
    corrected_loci_1 = [loci_1.chr, loci_1.start, loci_1.end]
    corrected_loci_2 = [loci_2.chr, loci_2.start, loci_2.end]

    for i in range(len(assignment_matching)):
        corrected_loci_1.append(loci_1['posterior'+str(assignment_matching[i][0])])
        corrected_loci_2.append(loci_2['posterior'+str(assignment_matching[i][1])])

    corrected_loci_1 = pd.concat(corrected_loci_1, axis=1)
    corrected_loci_2 = pd.concat(corrected_loci_2, axis=1)
    corrected_loci_2.columns = corrected_loci_1.columns
    
    return corrected_loci_1, corrected_loci_2
            

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

def curate_dataset(length_dist_1, length_dist_2, gmtk_params_1, gmtk_params_2, normalize_len=True):
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
    if px:
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
    else:
        sns.scatterplot(
        x=components[:,0], y=components[:,1], hue=clusters, s=100)

        plt.title("PCA plot")
        plt.xlabel("PC_1")
        plt.ylabel("PC_2")

        plt.show()

def tsne_plot(X, clusters, n_components):
    tsne = TSNE(
        n_components=n_components, perplexity=10, random_state=0)

    z = tsne.fit_transform(X)
    sns.scatterplot(
        x=z[:,0], y=z[:,1], hue=clusters, s=100)

    plt.title("TSNE plot")
    plt.xlabel("TSNE_1")
    plt.ylabel("TSNE_2")
    
    plt.show()

class Clusterer(object):
    def __init__(self, clustering_data, n_clusters, strategy='hier'):
        self.X = clustering_data
        self.k = n_clusters

        if strategy == 'hier':
            self.model = AgglomerativeClustering
        else:
            self.model = KMeans
            
    def fit_hierarchical(self, metric='euclidean', linkage='ward'):
        self.model = self.model(
                n_clusters=self.k, affinity=metric,
                linkage=linkage,
                compute_distances=True)
        
        self.predicted_labels =self.model.fit_predict(self.X)

        centroids = []
        for i in range(self.k):
            points_with_cluster_i_labels = []

            for j in range(len(self.predicted_labels)):
                if self.predicted_labels[j] == i:
                    points_with_cluster_i_labels.append(self.X.iloc[j, :])

            centroids.append(np.array(
                pd.DataFrame(points_with_cluster_i_labels).mean(axis=0)))

        self.model.cluster_centers_ = np.array(centroids)
        return self.model

    def dendrogram(self, plot=True, **kwargs):
        counts = np.zeros(self.model.children_.shape[0])
        print(counts)
        n_samples = len(self.model.labels_)

        for i, merge in enumerate(self.model.children_):
            current_count = 0
            for child_idx in merge:
                if child_idx < n_samples:
                    current_count += 1  
                else:
                    current_count += counts[child_idx - n_samples]
            counts[i] = current_count

        linkage_matrix = np.column_stack(
            [self.model.children_, self.model.distances_, counts]
        ).astype(float)

        print(linkage_matrix)

        # Plot the corresponding dendrogram
        if plot:
            dendrogram(linkage_matrix, **kwargs)
            plt.xlabel("Label")
            plt.ylabel("Distance")
            plt.show()

        return linkage_matrix

    def distance_based_agglomerative_clustering(self, distance_matrix):
        self.model = AgglomerativeClustering(
            affinity='precomputed', n_clusters=distance_matrix.shape[0], 
            linkage='complete', compute_distances=True).fit(distance_matrix)
        print(self.model.labels_)
        # print(self.model.distances_)

        self.dendrogram()

    def merge_clusters_by_order(self, unclustered_loci, clustering_obj, order=1): 
        '''
        merges clusters and their corresponding posterior value in 
        each bin. merging is done based on the order of hierarchy in 
        the dendogram.
        '''

        num_labels = unclustered_loci.shape[1]-3
        linkage_matrix = self.dendrogram(plot=False)

        new_loci = unclustered_loci.loc[:,('chr', 'start', 'end')]
        last_cluster_name = num_labels - 1


        for _ in range(order):
            merge_criteria = min(np.array(linkage_matrix[:, 3]))
            to_merge = {}
            for i in range(linkage_matrix.shape[0]):
                if linkage_matrix[i, 3] == merge_criteria:
                    last_cluster_name += 1
                    to_merge[last_cluster_name] = (linkage_matrix[i, 0], linkage_matrix[i, 1])


        cluster_labels = clustering_obj.labels_

        for i in range(len(cluster_labels)): 
            if "posterior{}".format(cluster_labels[i]) not in new_loci.columns:
                new_loci["posterior{}".format(cluster_labels[i])] = unclustered_loci["posterior{}".format(i)]
            else:
                new_loci["posterior{}".format(cluster_labels[i])] = \
                    new_loci["posterior{}".format(cluster_labels[i])] + unclustered_loci["posterior{}".format(i)]
        
        return new_loci
    
    def fit_kmeans(self):
        self.model = self.model(
            n_clusters=self.k, random_state=0)
        
        self.model.fit(self.X)
        self.predicted_labels = self.model.predict(self.X)
        return self.model


def update_labels_by_cluster(unclustered_loci, clustering_obj): 
    '''
    merges clusters and their corresponding posterior value in 
    each bin.  merging is done based on the number of clusters.
    '''

    new_loci = unclustered_loci.loc[:,('chr', 'start', 'end')]
    cluster_labels = clustering_obj.labels_

    for i in range(len(cluster_labels)): 
        if "posterior{}".format(cluster_labels[i]) not in new_loci.columns:
            new_loci["posterior{}".format(cluster_labels[i])] = unclustered_loci["posterior{}".format(i)]
        else:
            new_loci["posterior{}".format(cluster_labels[i])] = \
                new_loci["posterior{}".format(cluster_labels[i])] + unclustered_loci["posterior{}".format(i)]
    
    return new_loci

def compute_pairwise_centroid_distance(centroid_1, centroid_2):
    '''
    connect centroids using min eucleadian distance. can also be used to 
    compute the pairwise distance between any two sets of points
    '''
    
    distance_matrix = np.zeros(
        (centroid_1.shape[0], centroid_2.shape[0]))

    for i in range(centroid_1.shape[0]):
        for j in range(centroid_2.shape[0]):
            
            dist = pdist(
                np.array([centroid_1[i,:], centroid_2[j,:]]))
            
            distance_matrix[i,j] = dist

    distance_matrix = pd.DataFrame(distance_matrix)
    return distance_matrix

def Cooccurrence_matrix(loci_1, loci_2):
    cooc_mat = pd.DataFrame(
        np.zeros((int(loci_1.shape[1])-3, int(loci_2.shape[1])-3)), 
        index=loci_1.columns[3:],
        columns=loci_2.columns[3:])
    

    max_1posteri = find_max_posteri(loci_1)
    max_2posteri = find_max_posteri(loci_2)

    for i in range(len(loci_1)):
        try: # to handle some nan values in argmax vectors
            cooc_mat.loc[max_1posteri[i], max_2posteri[i]] += 1

        except:
            pass

    cooc_mat.index = [i.replace('posterior_cluster_','') for i in cooc_mat.index]
    cooc_mat.columns = [i.replace('posterior_cluster_','') for i in cooc_mat.columns]
    return cooc_mat

def match_evaluation(matrix, assignment_pairs):
    '''
    as the matching is performed using either:
    1. confusion matrix and hungarian algorithm
    2. clustering
    
    the approach should be specified to the function.'''
    
    probability_array = {}
    matched_sum = 0
    for i in range(matrix.shape[0]):
        probability_array[matrix.columns[i]] = float(matrix.iloc[i, assignment_pairs[i][1]]) /  (matrix.iloc[i,:].sum() + 1e-3)
        matched_sum += float(matrix.iloc[i, assignment_pairs[i][1]])    

    probability_array['all'] = matched_sum / np.array(matrix).sum()
    return pd.Series(probability_array)

