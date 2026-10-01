import math
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from ._cluster_matching import IoU_overlap, joint_overlap_prob
from ._reproducibility import overlap_heatmap


def perlabel_is_reproduced(reprod_report, loci_1):
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    labels = MAP1.unique()

    per_label_results = {}
    for k in labels:
        per_label_results[k] = [0, 0] # (True, False)
    
    for t in range(len(MAP1)):
        if reprod_report["is_repr"][t]:
            per_label_results[MAP1[t]][0] +=1

        else:
            per_label_results[MAP1[t]][1] +=1

    return per_label_results


def get_match(loci1, loci2, w):
    joint_overlap = joint_overlap_prob(loci1, loci2, w=w, symmetric=True)
    MAP1 = loci1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci2.iloc[:,3:].idxmax(axis=1)
    coverage1 = {k:len(MAP1.loc[MAP1 == k])/len(loci1) for k in loci1.columns[3:]}
    coverage2 = {k:len(MAP2.loc[MAP2 == k])/len(loci2) for k in loci2.columns[3:]}
    IoU = IoU_overlap(loci1, loci2, w=w, symmetric=False, soft=False)

    F1_rec = {}
    label_rec = {}
    coverage_rec = {}
    for k in IoU.index:
        F1_rec[k] = []
        label_rec[k] = []
        coverage_rec[k] = []

        sorted_k_vectors = IoU.loc[k,:].sort_values(ascending=False)
        i = 0
        while i < len(IoU.columns)-1:
            if i == 0:
                F1_rec[k].append(IoU.loc[k, sorted_k_vectors.index[i]])
                label_rec[k].append(sorted_k_vectors.index[i])
                coverage_rec[k].append(coverage2[sorted_k_vectors.index[i]])

            else:
                irange = list(range(i+1))
                merged_labels = "+".join([sorted_k_vectors.index[j] for j in irange])
                merged_joint = sum([joint_overlap.loc[k, sorted_k_vectors.index[j]] for j in irange])
                merged_coverage = sum([coverage2[sorted_k_vectors.index[j]] for j in irange])
                merged_IoU = (merged_joint) / (coverage1[k] + merged_coverage - merged_joint)

                F1_rec[k].append(merged_IoU)
                label_rec[k].append(merged_labels)
                coverage_rec[k].append(merged_coverage)

            i+=1

    corresp = {}
    for k in F1_rec.keys():
        corresp[k] = label_rec[k][F1_rec[k].index(max(F1_rec[k]))]

        if "+" in corresp[k]:
            corresp[k] = corresp[k].split("+")
        else:
            corresp[k] = [corresp[k]]
    
    return corresp


def dynamic_matches(loci1, loci2, w):
    r1vr2 = get_match(loci1, loci2, w)
    r2vr1 = get_match(loci2, loci1, w)
    
    num_labels = len(loci1.columns) - 3
    bidir_match = pd.DataFrame(np.zeros((num_labels, num_labels)), columns=loci2.columns[3:], index=loci1.columns[3:])

    for k in r1vr2.keys():
        for v in r1vr2[k]:
            bidir_match.loc[k, v] += 1

    for k in r2vr1.keys():
        for v in r2vr1[k]:
            bidir_match.loc[v, k] += 1

    overlap_heatmap(bidir_match)

    corresp = {}
    for k in r1vr2.keys():
        if k not in corresp.keys():
            corresp[k] = []

        for j in r2vr1.keys():

            if bidir_match.loc[k,j] == 2:
                corresp[k].append(j)
        
        if len(corresp[k]) == 0:
            for j in r2vr1.keys():

                if bidir_match.loc[k,j] == 1:
                    corresp[k].append(j)
    
    return corresp


def is_repr_MAP_with_prior(prior, loci_1, loci_2, ovr_threshold, window_bp, matching="static"):
    resolution = loci_1["end"][0] - loci_1["start"][0]
    window_bin = math.ceil(window_bp/resolution)

    enr_ovr = IoU_overlap(loci_1, loci_2, w=window_bp, symmetric=False, soft=False)
    
    if matching == "static":
        best_match = {i:enr_ovr.loc[i, :].idxmax() for i in enr_ovr.index} 
        above_threshold_match = {i: list(enr_ovr.loc[i, enr_ovr.loc[i, :]>ovr_threshold].index) for i in enr_ovr.index}

        for k in above_threshold_match.keys():
            if len(above_threshold_match[k]) == 0:
                above_threshold_match[k] = [best_match[k]]

    elif matching == "dynamic":
        above_threshold_match = dynamic_matches(loci_1, loci_2, w=window_bp)

    MAP1 = list(loci_1.iloc[:,3:].idxmax(axis=1))
    MAP2 = list(loci_2.iloc[:,3:].idxmax(axis=1))

    # reprod_report = loci_1.loc[:, ["chr", "start", "end"]]
    rep_rec = []
    for i in range(len(MAP1)):
        if prior[i] == True:
            rep_rec.append(True)

        else:
            if window_bin > 0:
                neighbor_i = MAP2[max(0, (i-window_bin)) : min((i+window_bin), (len(MAP2)-1))]
            else:
                neighbor_i = [MAP2[i]]

            if (set(above_threshold_match[MAP1[i]]) & set(neighbor_i)):
                rep_rec.append(True)

            else:
                rep_rec.append(False)
    
    return rep_rec


def calibrate(loci_1, loci_2, ovr_threshold, window_bp, numbins=500, matching="static", always_include_best_match=True):
    # based on the resolution determine the number of bins for w
    resolution = loci_1["end"][0] - loci_1["start"][0] 
    window_bin = math.ceil(window_bp/resolution)

    # get the overlap
    if matching == "static":
        enr_ovr = IoU_overlap(loci_1, loci_2, w=window_bp, symmetric=False, soft=False)

    elif matching == "dynamic":
        dyn_matching_map = dynamic_matches(loci_1, loci_2, w=window_bp)
        
    MAP2 = list(loci_2.iloc[:,3:].idxmax(axis=1))

    # determine the number of positions at each posterior bin
    strat_size = int(len(loci_1)/numbins)

    ################# GETTING CALIBRATIONS #################
    calibrations = {}

    for k in loci_1.columns[3:]:

        ####################### DEFINE MATCHES #######################
        if matching == "static":

            best_match = [enr_ovr.loc[k, :].idxmax()]
            above_threshold_match = list(enr_ovr.loc[k, enr_ovr.loc[k, :] > ovr_threshold].index)

            if len(above_threshold_match) == 0:
                if always_include_best_match:
                    above_threshold_match = [best_match[0]]

        elif matching == "dynamic":
            above_threshold_match = dyn_matching_map[k]

        ##############################################################

        vector_pair = pd.concat([loci_1[k], pd.Series(MAP2)], axis=1)
        vector_pair.columns=["pv1", "map2"]
        vector_pair = vector_pair.sort_values("pv1").reset_index()
        
        bins = []
        for b in range(0, vector_pair.shape[0], strat_size):
            # [bin_start, bin_end, num_values_in_bin, ratio_matched]

            subset_vector_pair = vector_pair.iloc[b:b+strat_size, :]
            subset_vector_pair = subset_vector_pair.reset_index(drop=True)
            
            matched = 0
            ####################### CHECK WINDOW #######################
            for i in range(len(subset_vector_pair)):

                if window_bin > 0 :
                    leftmost = max(0, (subset_vector_pair["index"][i] - window_bin))
                    rightmost = min((len(MAP2)-1), (subset_vector_pair["index"][i] + window_bin))

                    neighbors_i = MAP2[leftmost:rightmost]

                else:
                    neighbors_i = [MAP2[subset_vector_pair["index"][i]]]

                if (set(above_threshold_match) & set(neighbors_i)):
                    matched += 1
            #############################################################
            bins.append([
                float(subset_vector_pair["pv1"][subset_vector_pair.index[0]]), 
                float(subset_vector_pair["pv1"][subset_vector_pair.index[-1]]),
                len(subset_vector_pair),
                float(matched/len(subset_vector_pair))])
        
        bins = np.array(bins)
        polyreg = IsotonicRegression(
            y_min=float(np.array(bins[:, 3]).min()), y_max=float(np.array(bins[:, 3]).max()), 
            out_of_bounds="clip", increasing=True)
        
        polyreg.fit(
            np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)), 
            bins[:, 3])
         
        calibrations[k] = polyreg
    ########################################################
    # print("got the calibrations")
    calibrated_loci1 = [loci_1.loc[:, "chr"] ,loci_1.loc[:, "start"] ,loci_1.loc[:, "end"]]

    # here I'm recalibrating the posterior values based on the functions that we learned earlier
    for k in loci_1.columns[3:]:
        calibrated_loci1.append(pd.Series(calibrations[k].predict(np.reshape(np.array(loci_1.loc[:, k]), (-1,1)))))

    calibrated_loci1 = pd.concat(calibrated_loci1, axis=1)
    calibrated_loci1.columns = loci_1.columns
    return calibrated_loci1


def is_repr_posterior(
    loci_1, loci_2, ovr_threshold, window_bp, reproducibility_threshold=0.9, 
    matching="static", always_include_best_match=True, return_r=False, return_mean_r=False):
    MAP1 = list(loci_1.iloc[:,3:].idxmax(axis=1)) 

    calibrated_loci1 = calibrate(loci_1, loci_2, ovr_threshold, window_bp, matching=matching, always_include_best_match=always_include_best_match)
    # here I'm just looking at the calibrated score that is assigned to the MAP1 of each position (~>.9)
    calibrated_reproducibility = []
    for i in range(len(MAP1)):
        calibrated_reproducibility.append(calibrated_loci1.loc[i, MAP1[i]])

    if return_r:
        r = pd.concat([loci_1.chr, loci_1.start, loci_1.end, pd.Series(MAP1), pd.Series(calibrated_reproducibility)], axis=1)
        r.columns = ["chr", "start", "end", "MAP", "r_value"]
        return r

    else:
        # then based on the initial threshold that we had, we say a position is reproduced if it has a calibrated score of >threshold
        binary_isrep = []
        for t in calibrated_reproducibility:
            if t>=reproducibility_threshold:
                binary_isrep.append(True)
            else:
                binary_isrep.append(False)
        if return_mean_r:
            avg_r = np.mean(np.array(calibrated_reproducibility))
            return binary_isrep, avg_r
        else:
            return binary_isrep


def single_point_repr(
    loci_1, loci_2, ovr_threshold, window_bp, posterior=False, reproducibility_threshold=0.9, return_mean_r=False):
    if posterior:
        """
        for all labels in R1:
            creat a calibration curve according to t, w
        
        according to the calibrated p, 
        """
        return is_repr_posterior(
            loci_1, loci_2, ovr_threshold, window_bp, reproducibility_threshold=reproducibility_threshold, return_mean_r=return_mean_r)

    else:
        return is_repr_MAP_with_prior(
            [False for _ in range(len(loci_1))], 
            loci_1, loci_2, ovr_threshold, window_bp)


def keep_reproducible_annotations(loci, bool_reprod_report):
    bool_reprod_report = pd.concat(
        [loci["chr"], loci["start"], loci["end"], pd.Series(bool_reprod_report)], axis=1)
    bool_reprod_report.columns = ["chr", "start", "end", "is_repr"]

    return loci.loc[(bool_reprod_report["is_repr"] == True), :].reset_index(drop=True)


def condense_segments(coordMAP):
    coordMAP = coordMAP.values.tolist()
    i=0
    while i+1 < len(coordMAP):
        if coordMAP[i][0] == coordMAP[i+1][0] and coordMAP[i][3] == coordMAP[i+1][3]:
            coordMAP[i+1][1] = coordMAP[i][1]
            del coordMAP[i]

        else:
            i += 1
    
    return pd.DataFrame(coordMAP, columns=["chr", "start", "end", "MAP"])


def write_MAPloci_in_BED(loci, savedir):
    MAP = loci.iloc[:,3:].idxmax(axis=1)
    coordMAP = pd.concat([loci.iloc[:, :3], MAP], axis=1)
    coordMAP.columns = ["chr", "start", "end", "MAP"]
    denseMAP = condense_segments(coordMAP)

    coordMAP.to_csv(savedir + "/confident_segments.bed", sep='\t', header=False, index=False)
    denseMAP.to_csv(savedir + "/confident_segments_dense.bed", sep='\t', header=False, index=False)
