from matplotlib import pyplot as plt
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import euclidean
from scipy.stats import pearsonr
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import KMeans, AgglomerativeClustering
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn import preprocessing
from sklearn.cluster import KMeans
from scipy.spatial.distance import pdist
from sklearn.metrics import cohen_kappa_score
import math, os
from matplotlib.colors import LinearSegmentedColormap

def overall_overlap_ratio(loci_1, loci_2, w=0):
    joint = joint_overlap_prob(loci_1, loci_2, w=w, symmetric=False) * len(loci_1)
    sum_matched = np.sum([np.max(joint.loc[i, :]) for i in joint.index])
    return sum_matched / len(loci_1)

def overlap_matrix(loci_1, loci_2, type="IoU"):
    if type == "IoU":
        return IoU_overlap(loci_1, loci_2)
    elif type == "ck":
        return Cohens_Kappa_matrix(loci_1, loci_2)
    elif type == "enr":
        return enrichment_of_overlap_matrix(loci_1, loci_2, OE_transform=True)
    elif type == "conditional":
        return enrichment_of_overlap_matrix(loci_1, loci_2, OE_transform=False)
    elif type == "hard_joint":
        return joint_overlap_prob(loci_1, loci_2, w=0)
    elif type == "soft_joint":
        return soft_joint_prob(loci_1, loci_2)
    
def soft_joint_prob(loci_1, loci_2):
    """
    this kind of overlap, takes the posterior prob into account. 
    the posterior prob should not be in logit form.
        if 0<p<1:
            pass
        else:
            p = sigmoid(p)
    """

    num_labels = len(loci_1.columns) -3
    joint = pd.DataFrame(np.zeros((num_labels, num_labels)), columns=loci_2.columns[3:], index=loci_1.columns[3:])

    for k in loci_1.columns[3:]:
        for j in loci_2.columns[3:]:
            soft_overlap = sum(np.array(loci_1.loc[:, k]) * np.array(loci_2.loc[:, j])) / len(loci_1)
            joint.loc[k, j] = soft_overlap
    
    return joint

def soft_coverage(loci_1, loci_2):
    coverage1 = {k:sum(loci_1.loc[:,k])/len(loci_1) for k in loci_1.columns[3:]}
    coverage2 = {k:sum(loci_2.loc[:,k])/len(loci_2) for k in loci_2.columns[3:]}
    return coverage1, coverage2

def joint_overlap_prob(loci_1, loci_2, w=0, symmetric=True):
    num_labels1 = len(loci_1.columns) -3
    num_labels2 = len(loci_2.columns) -3
    resolution = loci_1["end"][0] - loci_1["start"][0] 
    w = math.ceil(w/resolution)

    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

    observed_overlap = {} 

    for k in loci_1.columns[3:]:
        for j in loci_2.columns[3:]:
            observed_overlap[str(k + "|" + j)] = 0

    oo_mat = pd.DataFrame(np.zeros((num_labels1, num_labels2)), columns=loci_2.columns[3:], index=loci_1.columns[3:])

    MAP1 = list(MAP1)
    MAP2 = list(MAP2)
                
    if w == 0:
        for i in range(len(MAP1)):
            k = MAP1[i]
            j = MAP2[i]
            observed_overlap[str(k + "|" + j)] += 1

    elif symmetric==False:
    ############################### NON-SYMMETRIC ##################################
        for i in range(len(MAP1)):
            k = MAP1[i]
            i_neighbors = MAP2[max(0, i-w) : min(i+w, len(MAP2)-1)]
            for j in set(i_neighbors):
                observed_overlap[str(k + "|" + j)] += 1
        
    else:
    ################################# SYMMETRIC ####################################
        for i in range(len(MAP1)):
            R1_window =  MAP1[max(0, i-w) : min(i+w, len(MAP1)-1)]
            R2_window =  MAP2[max(0, i-w) : min(i+w, len(MAP2)-1)]
            for k in R1_window:
                for j in R2_window:
                    observed_overlap[str(k + "|" + j)] += float(1/(len(R1_window)*len(R2_window)))

    for p in observed_overlap.keys():
        oo_mat.loc[p.split("|")[0], p.split("|")[1]] = observed_overlap[p]
    
    return oo_mat / len(loci_1)

def IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False, overlap_coeff=False):
    num_labels1 = len(loci_1.columns) - 3
    num_labels2 = len(loci_2.columns) - 3
    
    IoU = pd.DataFrame(np.zeros((num_labels1, num_labels2)), columns=loci_2.columns[3:], index=loci_1.columns[3:])

    if soft and w == 0:
        joint = soft_joint_prob(loci_1, loci_2)
        coverage1, coverage2 = soft_coverage(loci_1, loci_2)

    else:
        joint = joint_overlap_prob(loci_1, loci_2, w=w, symmetric=symmetric)

        MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
        MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

        coverage1 = {k:len(MAP1.loc[MAP1 == k])/len(loci_1) for k in loci_1.columns[3:]}
        coverage2 = {k:len(MAP2.loc[MAP2 == k])/len(loci_2) for k in loci_2.columns[3:]}

        # coverage1 = {k: joint.loc[k,:].sum() for k in joint.index}
        # coverage2 = {k: joint.loc[:,k].sum() for k in joint.columns}

    for A in loci_1.columns[3:]:
        for B in loci_2.columns[3:]:
            if overlap_coeff:
                IoU.loc[A, B] = (joint.loc[A, B])/np.min([coverage1[A], coverage2[B]])

            else:
                IoU.loc[A, B] = (joint.loc[A, B])/(coverage1[A] + coverage2[B] - ((joint.loc[A, B])))

    return IoU

def __F1_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False):
    num_labels = len(loci_1.columns) - 3
    
    F1 = pd.DataFrame(np.zeros((num_labels, num_labels)), columns=loci_2.columns[3:], index=loci_1.columns[3:])

    if soft and w == 0:
        joint = soft_joint_prob(loci_1, loci_2)
        coverage1, coverage2 = soft_coverage(loci_1, loci_2)

    else:
        joint = joint_overlap_prob(loci_1, loci_2, w=w, symmetric=symmetric)

        MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
        MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

        coverage1 = {k:len(MAP1.loc[MAP1 == k])/len(loci_1) for k in loci_1.columns[3:]}
        coverage2 = {k:len(MAP2.loc[MAP2 == k])/len(loci_2) for k in loci_2.columns[3:]}

        # coverage1 = {k: joint.loc[k,:].sum() for k in joint.index}
        # coverage2 = {k: joint.loc[:,k].sum() for k in joint.columns}

    for A in loci_1.columns[3:]:
        for B in loci_2.columns[3:]:
            F1.loc[A, B] = (2 * joint.loc[A, B]) / (coverage1[A] + coverage2[B])

    return F1

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################

def find_max_posteri(loci):
    # loci_posteri =  loci.iloc[:, 3:]
    loci.iloc[:, 3:] = loci.iloc[:, 3:].astype('float')

    max_posteri =  loci.iloc[:, 3:].idxmax(axis=1)
    return max_posteri

def make_symmetric_matrix(matrix0):
    '''
    this function gets an asymmetric confusion/distance/similarity matrix
    and makes a symmetric one based on the hypothetical overlap of merging 
    all possible pairs of labels.'''
    
    num_labels = matrix0.shape[0]
    matrix1 = pd.DataFrame(
        np.zeros((num_labels,num_labels)), 
        columns=matrix0.columns, index=matrix0.index)
    
    for i in matrix0.columns:
        for j in matrix0.index:

            if i == j:
                matrix1.loc[i, j] = matrix0.loc[i, j]

            else:
                #overlap corresponds to the sum of two values from matrix0 
                overlap = matrix0.loc[i, j] + matrix0.loc[j, i]

                # matrix1[i, j] == matrix1[j, i] == overlap 
                matrix1.loc[i, j] = overlap
                matrix1.loc[j, i] = overlap
    
    return matrix1

def Cohens_Kappa_matrix(loci_1, loci_2):
    num_labels = len(loci_1.columns) -3
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

    coverage1 = {k:len(MAP1.loc[MAP1 == k])/len(loci_1) for k in loci_1.columns[3:]}
    coverage2 = {k:len(MAP2.loc[MAP2 == k])/len(loci_2) for k in loci_2.columns[3:]}

    ck_mat = pd.DataFrame(
        np.zeros((num_labels, num_labels)), 
        columns=loci_2.columns[3:], index=loci_1.columns[3:])

    for k in loci_1.columns[3:]:
        for l in loci_2.columns[3:]:
            o = 0
            for i in range(len(MAP1)):
                if MAP1[i] == k and MAP2[i] == l:
                    o+=1
            
            ck_mat.loc[k, l] = cohen_kappa_score( 
                list(MAP1 == k),
                list(MAP2 == l)
            )

    return ck_mat

def enrichment_of_overlap_matrix(loci_1, loci_2, OE_transform=True):
    num_labels = len(loci_1.columns) -3
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

    coverage1 = {k: len(MAP1.loc[MAP1 == k]) / len(MAP1) for k in loci_1.columns[3:]}
    coverage2 = {k: len(MAP2.loc[MAP2 == k]) / len(MAP2) for k in loci_2.columns[3:]}

    observed_overlap = {} # np.zeros((num_labels, num_labels))
    expected_overlap = {} # np.zeros((num_labels, num_labels))

    for k in loci_1.columns[3:]:
        for j in loci_2.columns[3:]:
            observed_overlap[str(k + "|" + j)] = 0
            expected_overlap[str(k + "|" + j)] = coverage1[k] * coverage2[j] * len(MAP1)

    for i in range(MAP1.shape[0]):
        k = MAP1[i]
        j = MAP2[i]

        observed_overlap[str(k + "|" + j)] += 1

    oo_mat = pd.DataFrame(np.zeros((num_labels, num_labels)), columns=loci_2.columns[3:], index=loci_1.columns[3:])
    eo_mat = pd.DataFrame(np.zeros((num_labels, num_labels)), columns=loci_2.columns[3:], index=loci_1.columns[3:])

    for p in observed_overlap.keys():
        oo_mat.loc[p.split("|")[0], p.split("|")[1]] = observed_overlap[p]
        eo_mat.loc[p.split("|")[0], p.split("|")[1]] = expected_overlap[p]
    
    if OE_transform:
        epsilon = 1e-3

        return np.log(
            (oo_mat + epsilon) / (eo_mat + epsilon)
        )
    
    else:
        for k in oo_mat.index:
            oo_mat.loc[k, :] = oo_mat.loc[k, :] / (coverage1[k]*len(loci_1))
        return oo_mat

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

def connect_bipartite(loci_1, loci_2, assignment_matching, mnemon=True):
    corrected_loci_1 = loci_1.iloc[:, :3]
    corrected_loci_2 = loci_2.iloc[:, :3]

    for i in range(len(assignment_matching)):
        new_name = str(assignment_matching[i][0]) + "|" + str(assignment_matching[i][1])

        if mnemon:
            corrected_loci_1[new_name] = loci_1['posterior'+str(assignment_matching[i][0].split('_')[0])]
            corrected_loci_2[new_name] = loci_2['posterior'+str(assignment_matching[i][1].split('_')[0])]
        else:
            corrected_loci_1[new_name] = loci_1['posterior'+str(assignment_matching[i][0])]
            corrected_loci_2[new_name] = loci_2['posterior'+str(assignment_matching[i][1])]
    
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

def correspondence_based_on_emission(rep_dir1, rep_dir2, outdir, saga="chmm", metric="cosine"):
    if saga == "segway":
        emis_1 = pd.read_csv(rep_dir1 + "/gmtk_parameters/gmtk_parameters.stats.csv").drop("Unnamed: 0", axis=1)
        emis_2 = pd.read_csv(rep_dir2 + "/gmtk_parameters/gmtk_parameters.stats.csv").drop("Unnamed: 0", axis=1)

    elif saga == "chmm":
        for l in os.listdir(rep_dir1):
            if "emissions" in l and ".txt" in l:
                emis_1 = pd.read_csv(f"{rep_dir1}/{l}", sep="\t").drop("State (Emission order)" ,axis=1)

        for l in os.listdir(rep_dir2):
            if "emissions" in l and ".txt" in l:
                emis_2 = pd.read_csv(f"{rep_dir2}/{l}", sep="\t").drop("State (Emission order)" ,axis=1)  

    emis_1 = emis_1.sort_index(axis=1)
    emis_2 = emis_2.sort_index(axis=1)

    boundaries = [x for x in list(np.linspace(0, 1, 20))] + [1] # custom boundaries
    hex_colors = sns.light_palette('navy', n_colors=len(boundaries) * 2 + 2, as_cmap=False).as_hex()
    hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]
    colors=list(zip(boundaries, hex_colors))
    custom_color_map = LinearSegmentedColormap.from_list(
        name='custom_navy',
        colors=colors)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    sns.heatmap(
        emis_1.astype(float), annot=False, fmt=".1f",
        linewidths=0.01,  cbar=True, annot_kws={"size": 6}, 
        cmap=custom_color_map, ax=ax1)
    
    ax1.set_title("Base")
    
    sns.heatmap(
        emis_2.astype(float), annot=False, fmt=".1f",
        linewidths=0.01,  cbar=True, annot_kws={"size": 6}, 
        cmap=custom_color_map, ax=ax2)

    ax2.set_title("Verification")

    plt.tight_layout()
    plt.savefig('{}/emissions.pdf'.format(outdir), format='pdf')
    plt.savefig('{}/emissions.svg'.format(outdir), format='svg')

    sim_mat = np.zeros((emis_1.shape[0], emis_2.shape[0]))

    for i in range(emis_1.shape[0]):
        for j in range(emis_2.shape[0]):

            if metric == "cosine":
                sim_mat[i,j] = cosine_similarity(np.array(emis_1.iloc[i,:]).reshape(1,-1), np.array(emis_2.iloc[j,:]).reshape(1,-1))

            elif metric == "eucl":
                sim_mat[i,j] = euclidean(np.array(emis_1.iloc[i,:]), np.array(emis_2.iloc[j,:]))

            elif metric == "correl":
                corr, _ = pearsonr(np.array(emis_1.iloc[i,:]), np.array(emis_2.iloc[j,:]))

                if math.isnan(corr):
                    sim_mat[i,j] = 0
                else:
                    sim_mat[i,j] = corr

    if metric == "eucl":
        # similarity = exp(-gamma * distance^2)

        sim_mat = np.exp(-1 * sim_mat**2)

    return sim_mat
