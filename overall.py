from _cluster_matching import *
import math, scipy
import os
from _reproducibility import *
from datetime import datetime

def is_reproduced(loci_1, loci_2, enr_threshold=3, window_bp=500, raw_overlap_axis=False):
    """
    get enr_ovr matrix
    get MAPs

    for t in map1:
        statement1 = bool(position t in map2 is the best matching label)
        statement2 = bool(position t in map2 has a label with enr_ovr > enr_threshold)
        statement3 = [statement1 or statement2 for (t - 500bp) & (t + 500bp)]

        if statement1 | statement2 | statement3:
            t = "reproduced"
    """
    resolution = loci_1["end"][0] - loci_1["start"][0]
    window_bin = math.ceil(window_bp/resolution)

    if raw_overlap_axis:
        enr_ovr = IoU_overlap(loci_1, loci_2, w=window_bp, symmetric=False, soft=False)
    else:
        enr_ovr = IoU_overlap(loci_1, loci_2, w=window_bp, symmetric=False, soft=False)

    best_match = {i:enr_ovr.loc[i, :].idxmax() for i in enr_ovr.index} 
    above_threshold_match = {i: list(enr_ovr.loc[i, enr_ovr.loc[i, :]>enr_threshold].index) for i in enr_ovr.index}

    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

    reprod_report = loci_1.loc[:, ["chr", "start", "end"]]
    reprod_report["is_repr"] = np.zeros((len(MAP1)))

    reprod_report = reprod_report.values.tolist()

    for t in range(len(MAP1)):
        statement1 = bool(MAP2[t] == best_match[MAP1[t]])
        statement2 = bool(MAP2[t] in above_threshold_match[MAP1[t]])
        statement3 = False
        
        leftmost = max(0, (t - window_bin))
        rightmost = min(len(MAP1), (t + window_bin))

        window_around_t = MAP2[leftmost:rightmost]

        if best_match[MAP1[t]] in list(window_around_t):
            reprod_report[t][3] = True
            
        elif len(list(set(list(window_around_t)) & set(above_threshold_match[MAP1[t]]))) > 0:
            statement3 = True
        
        else:
            statement3 = False

        if statement1 or statement2 or statement3:
            reprod_report[t][3] = True
        else:
            reprod_report[t][3] = False

    # print("reproduced {}".format(float(r/len(MAP1))))
    reprod_report = pd.DataFrame(reprod_report, columns=["chr", "start", "end", "is_repr"])

    # print(len(reprod_report.loc[reprod_report["is_repr"]==True])/ len(reprod_report))
    return reprod_report
        
def is_reproduced_posterior(loci_1, loci_2, enr_threshold=2, window_bp=500, raw_overlap_axis=False, numbins=100, allow_w=True):
    """
    get enr_ovr matrix
    get MAPs
    
    let m be the number of labels
    obtain an m*m posterior calibration functions

    for t in map1:
        statement1 = bool(position t in map2 is the best matching label)
        statement2 = bool(position t in map2 has a label with enr_ovr > enr_threshold)
        statement3 = [statement1 or statement2 for (t - 500bp) & (t + 500bp)]

        if statement1 | statement2 | statement3:
            t = "reproduced"
    """

    resolution = loci_1["end"][0] - loci_1["start"][0]
    window_bin = math.ceil(window_bp/resolution)

    enr_ovr = IoU_overlap(loci_1, loci_2, w=window_bp, symmetric=False, soft=False)
    MP1 = loci_1.iloc[:,3:].max(axis=1) #posterior values 1
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1) # MAP label 1
    
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
    
    above_threshold_coverage = {}
    calibration_funcs = {}
    for i in range(len(loci_1.columns[3:])):
        if loci_1.columns[3:][i] not in calibration_funcs.keys() and loci_1.columns[3:][i] not in above_threshold_coverage.keys():
            calibration_funcs[loci_1.columns[3:][i]] = {}

            best_match = enr_ovr.loc[loci_1.columns[3:][i], :].idxmax() #{j:enr_ovr.loc[j, :].idxmax() for j in enr_ovr.index} 
            above_threshold_match = list(enr_ovr.loc[loci_1.columns[3:][i], enr_ovr.loc[loci_1.columns[3:][i], :]>enr_threshold].index) #{i: list(enr_ovr.loc[i, enr_ovr.loc[i, :]>enr_threshold].index) for i in enr_ovr.index}
            if len(above_threshold_match) == 0:
                above_threshold_match = [best_match]

            print(loci_1.columns[3:][i], above_threshold_match, window_bin)

            calibration_funcs[loci_1.columns[3:][i]],  above_threshold_coverage[loci_1.columns[3:][i]]= \
                    calibrate(
                        loci_1.loc[:, loci_1.columns[3:][i]], MAP2,
                        k=above_threshold_match, 
                        num_bins=numbins, 
                        allow_w=allow_w, 
                        window_bin=window_bin)

        # for j in range(len(loci_2.columns[3:])):
        #     if loci_2.columns[3:][j] not in calibration_funcs[loci_1.columns[3:][i]].keys():
                

        #         # print(loci_1.columns[3:][i], loci_2.columns[3:][j])

        #         calibration_funcs[loci_1.columns[3:][i]][loci_2.columns[3:][j]] = \
        #             calibrate(
        #                 loci_1.loc[:, loci_1.columns[3:][i]], MAP2,
        #                 k=loci_2.columns[3:][j], 
        #                 num_bins=numbins, 
        #                 allow_w=allow_w, 
        #                 window_bin=window_bin)

    posterior_enrichment = []
    # coverage_1 = {k:float(len(MAP1.loc[MAP1 == k]) / len(MAP2)) for k in MAP1.unique()}
    coverage_2 = above_threshold_coverage #{k:float(len(MAP2.loc[MAP2 == k]) / len(MAP2)) for k in MAP2.unique()}

    if raw_overlap_axis:
        for t in range(len(loci_1)): 
            enr_val =calibration_funcs[MAP1[t]].predict(np.reshape(np.array(MP1[t]), (1,1)))[0]
            batchsize = int(len(loci_1)/numbins)
            e = coverage_2[MAP1[t]] * batchsize

            overlap_ratio = (np.exp(enr_val)*e) / batchsize
            
            posterior_enrichment.append(overlap_ratio)

    else:
        posterior_enrichment = calibration_funcs[MAP1[t]].predict(np.reshape(np.array(MP1), (-1,1)))
            
    posterior_enrichment = np.array(posterior_enrichment)
    # print(posterior_enrichment)

    # plt.hist(posterior_enrichment, bins=100)

    is_repr = []
    for t in range(len(posterior_enrichment)):
        leftmost = max(0, t-window_bin)
        rightmost = min(len(posterior_enrichment), t+window_bin)

        max_score = max(posterior_enrichment[leftmost:rightmost])

        if max_score >= enr_threshold:
            is_repr.append(True)
        else:
            is_repr.append(False)

    reprod_report = pd.concat([loci_1["chr"], loci_1["start"], loci_1["end"], pd.Series(is_repr)], axis=1)
    reprod_report.columns = ["chr", "start", "end", "is_repr"]
    # score = len(reprod_report.loc[reprod_report["is_repr"]==True])/ len(reprod_report)

    return reprod_report

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

def reprod_score(loci_1, loci_2, window_bp=500):
    """
    get enr_ovr matrix
    get MAPs

    for t in map1:
        in a window of size w 
            BMS = what is the best match score (BMS) -> max of enr_ovr[R2] within that window
            d_BMS = what is the distance to BMS

            Score_t = BMS/d_BMS 
    """

    resolution = loci_1["end"][0] - loci_1["start"][0]
    window_bin = math.ceil(window_bp/resolution)

    enr_ovr = IoU_overlap(loci_1, loci_2, w=window_bp, symmetric=False, soft=False)

    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

    reprod_report = loci_1.loc[:, ["chr", "start", "end"]]
    reprod_report["repr_score"] = np.zeros((len(MAP1)))
    reprod_report = reprod_report.values.tolist()

    for t in range(len(MAP1)):
        in_window = {}
        in_window[0] = enr_ovr.loc[MAP1[t], MAP2[t]]

        for w in range(1, window_bin+1):
            if (t - w) >= 0 and (t + w) <= len(MAP1)-1 :
                in_window[-1 * w] = enr_ovr.loc[MAP1[t], MAP2[t-w]]
                in_window[w] = enr_ovr.loc[MAP1[t], MAP2[t+w]]
        
        d_BMS = max(in_window, key=in_window.get)
        BMS = max(in_window.values())

        reprod_report[t][3] = (BMS) / ((abs(d_BMS))+1)
    
    reprod_report = pd.DataFrame(reprod_report, columns=["chr", "start", "end", "repr_score"])
    return reprod_report

def max_posterior_to_repr_score(loci_1, reprod_report):
    MP1 = loci_1.iloc[:,3:].max(axis=1)

    pcorrel = scipy.stats.pearsonr(MP1, reprod_report["repr_score"])[0]
    scorrel = scipy.stats.spearmanr(MP1, reprod_report["repr_score"])[0]

    return pcorrel, scorrel
    
def contour_isrep(loci1, loci2, savedir, posterior=True, raw_overlap_axis=True):
    rec = []
    perlabel_rec = {k:[] for k in loci1.columns[3:]}

    if raw_overlap_axis:
        t_range = [0, 120, 20]
    else:
        t_range = [0, 400, 80]

    for w in range(100, 3100, 500):
        for t in range(t_range[0], t_range[1], t_range[2]):

            thresh = float(t)/100   
            if posterior:
                reprod_report = is_reproduced_posterior(
                    loci1, loci2, enr_threshold=thresh, window_bp=w, 
                    raw_overlap_axis=raw_overlap_axis, allow_w=True)
            else:
                reprod_report = is_reproduced(
                    loci1, loci2, enr_threshold=thresh, window_bp=w, 
                    raw_overlap_axis=raw_overlap_axis)

            score = len(reprod_report.loc[reprod_report["is_repr"]==True])/ len(reprod_report)
            rec.append([w, thresh, score])

            lab_rep = perlabel_is_reproduced(reprod_report, loci1)
            for k, v in lab_rep.items():
                perlabel_rec[k].append([w, thresh, v[0]/ (v[0]+v[1])]) #w, thresh, score

    rec = np.array(rec)
    x, y, z = rec[:,0], rec[:,1], rec[:,2]
    rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])
    list_of_ws = np.unique(rec[:, 0])
    cmap = get_cmap('inferno')
    gradients = np.linspace(0, 1, len(list_of_ws))

    for i in range(len(list_of_ws)):
        w = int(list_of_ws[i])
        subset_runs = rec_df.loc[rec_df["w"] == w, :].sort_values(by="t", ascending=False).reset_index(drop=True)
        
        plt.plot(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), label="w = {} bp".format(w))

    if raw_overlap_axis:
        plt.xlabel("Overlap threshold")
    else:
        plt.xlabel("Enrichment of Overlap threshold")

    plt.ylabel("Reproducibility")
    plt.xticks(rotation=45, fontsize=7)
    plt.title("Reproducibility -- General")
    plt.legend()
    plt.tight_layout()
    
    if os.path.exists(savedir+"/contours/") == False:
        os.mkdir(savedir+"/contours/")

    plt.savefig('{}/reprod_lines_general.pdf'.format(savedir+"/contours"), format='pdf')
    plt.savefig('{}/reprod_lines_general.svg'.format(savedir+"/contours"), format='svg')
    plt.clf()

    X, Y = np.meshgrid(x, y)

    Z = np.zeros(X.shape)
    
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            corresponding_line = rec_df.loc[(rec_df["w"] == X[i,j]) & (rec_df["t"]== Y[i,j]), :]
            Z[i, j] = float(corresponding_line["s"])

    plt.contourf(X, Y, Z, cmap='RdGy_r')
    plt.colorbar()
    if raw_overlap_axis:
        plt.ylabel("Overlap threshold")
    else:
        plt.ylabel("Enrichment of Overlap threshold")

    plt.xlabel("Window size")
    plt.xticks(rotation=45, fontsize=7)
    plt.title("Reproducibility -- {}".format(k))
    plt.tight_layout()
    plt.savefig('{}/reprod_contour_{}.pdf'.format(savedir+"/contours", k), format='pdf')
    plt.savefig('{}/reprod_contour_{}.svg'.format(savedir+"/contours", k), format='svg')
    plt.clf()

    for k in perlabel_rec.keys():
        rec = np.array(perlabel_rec[k])

        x, y, z = rec[:,0], rec[:,1], rec[:,2]
        rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])

        list_of_ws = np.unique(rec[:, 0])
        cmap = get_cmap('inferno')
        gradients = np.linspace(0, 1, len(list_of_ws))

        for i in range(len(list_of_ws)):
            w = int(list_of_ws[i])
            subset_runs = rec_df.loc[rec_df["w"] == w, :].sort_values(by="t", ascending=False).reset_index(drop=True)
            
            plt.plot(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), label="w = {} bp".format(w))

        if raw_overlap_axis:
            plt.xlabel("Overlap threshold")
        else:
            plt.xlabel("Enrichment of Overlap threshold")

        plt.ylabel("Reproducibility")
        plt.xticks(rotation=45, fontsize=7)
        plt.title("Reproducibility -- {}".format(k))
        plt.tight_layout()
        plt.legend()
        plt.savefig('{}/reprod_lines_{}.pdf'.format(savedir+"/contours", k), format='pdf')
        plt.savefig('{}/reprod_lines_{}.svg'.format(savedir+"/contours", k), format='svg')
        plt.clf()

        X, Y = np.meshgrid(x, y)

        Z = np.zeros(X.shape)
        
        for i in range(Z.shape[0]):
            for j in range(Z.shape[1]):
                corresponding_line = rec_df.loc[(rec_df["w"] == X[i,j]) & (rec_df["t"]== Y[i,j]), :]
                Z[i, j] = float(corresponding_line["s"])

        plt.contourf(X, Y, Z, cmap='RdGy_r')
        plt.colorbar()
        if raw_overlap_axis:
            plt.ylabel("Overlap threshold")
        else:
            plt.ylabel("Enrichment of Overlap threshold")

        plt.xlabel("Window size")
        plt.xticks(rotation=45, fontsize=7)
        plt.title("Reproducibility -- {}".format(k))
        plt.tight_layout()
        plt.savefig('{}/reprod_contour_{}.pdf'.format(savedir+"/contours", k), format='pdf')
        plt.savefig('{}/reprod_contour_{}.svg'.format(savedir+"/contours", k), format='svg')
        plt.clf()
 
##########################################################################################################################################
################################################################ MODIFIED ################################################################
##########################################################################################################################################

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

def single_point_delta_NMI_repr(loci_1, loci_2, ovr_threshold, window_bp, repr_threshold, posterior=True, matching="static"):

    if posterior:
        isrep1 = is_repr_posterior(loci_1, loci_2, ovr_threshold, window_bp, reproducibility_threshold=repr_threshold, matching=matching)
        isrep2 = is_repr_posterior(loci_2, loci_1, ovr_threshold, window_bp, reproducibility_threshold=repr_threshold, matching=matching)
    else:
        isrep1 = is_repr_MAP_with_prior(
            [False for _ in range(len(loci_1))], 
            loci_1, loci_2, ovr_threshold, window_bp, matching="static")
        
        isrep2 = is_repr_MAP_with_prior(
            [False for _ in range(len(loci_1))], 
            loci_2, loci_1, ovr_threshold, window_bp, matching="static")

    # at this point we can get the intersection of isrep1 and isrep2 as the positions that are considered 
    is_rep_intersec = []

    for i in range(len(isrep1)):
        if isrep1[i]==isrep2[i]==1:
            is_rep_intersec.append(True)
        else:
            is_rep_intersec.append(False)

    calibrated_loci1 = keep_reproducible_annotations(loci_1, is_rep_intersec) 
    calibrated_loci2 = keep_reproducible_annotations(loci_2, is_rep_intersec) 

    pre_calib_NMI = normalized_mutual_information(loci_1, loci_2, soft=False)
    post_calib_NMI = normalized_mutual_information(calibrated_loci1, calibrated_loci2, soft=False)

    return float(post_calib_NMI - pre_calib_NMI)

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

def rescale_calibrated_loci_to_prob(calibrated_loci):
    calloci_array = np.array(calibrated_loci.iloc[:,3:])

    for i in range(calloci_array.shape[0]):
        calloci_array[i,:] = calloci_array[i,:] / np.sum(calloci_array[i,:])

    calibrated_loci.iloc[:,3:] = calloci_array
    return calibrated_loci

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

def OvrWind_contour(
    loci1, loci2, savedir, w_range=[0, 4000, 200], t_range=[0, 11, 1], posterior=True, repr_threshold=0.9):
    # t is the overlap threshold (used to define matches)
    """
    initialize wrange and trange and reverse trange to start from the most strict to easiest
    initialize a matrix of size #t * #w and put a vector of [False for i in range len(loci1)] as each element

    for t in trange
        for w in wrange
            get_reproducibility of A[t, w] for all positions where rep is still False 
            for all reproducibile positions:
                change that position to True if that cell in A has bigger window or lower threshold
    
    after going through A once, return A[t,w] = len(A[t,w]==True) / len(loci1)
    """

    t_range = [i/10 for i in range(t_range[0], t_range[1], t_range[2])]
    t_range = t_range[::-1]

    w_range = list(range(w_range[0], w_range[1], w_range[2]))

    # initialize A
    A = {}
    for t in t_range:
        for w in w_range:
            if t in A.keys():
                A[t][w] = [False for i in range(len(loci1))]
            else:
                A[t] = {}
                A[t][w] = [False for i in range(len(loci1))]
    
    for t in t_range:
        for w in w_range:
            """
            get a subset of loci1 where A[t,w] == False
            run reproducibility_test on the subset of loci1
            update A[t,w] with new reprod_results
            """
            t0 = datetime.now()

            if posterior:
                updated_repr = is_repr_posterior(
                    loci1, loci2, t, w, reproducibility_threshold=repr_threshold, 
                    matching="static", always_include_best_match=False)

            else:
                updated_repr = is_repr_MAP_with_prior(A[t][w], loci1, loci2, t, w, matching="static")

            print(datetime.now() - t0)
            
            # now update all of A with new info
            for tt in t_range:
                for ww in w_range:
                    if tt <= t and ww >= w:
                        A[tt][ww] = updated_repr
    
    rec = []
    perlabel_rec = {k:[] for k in loci1.columns[3:]}
    for thresh in A.keys():
        for w in A[thresh].keys():
            reprod_report = pd.DataFrame(A[thresh][w], columns=["is_repr"])

            score = len(reprod_report.loc[reprod_report["is_repr"]==True])/ len(reprod_report)
            rec.append([w, thresh, score])

            lab_rep = perlabel_is_reproduced(reprod_report, loci1)
            for k, v in lab_rep.items():
                perlabel_rec[k].append([w, thresh, v[0]/ (v[0]+v[1])])
    
    rec = np.array(rec)
    x, y, z = rec[:,0], rec[:,1], rec[:,2]
    rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])
    list_of_ws = np.unique(rec[:, 0])
    cmap = get_cmap('inferno')
    gradients = np.linspace(0, 1, len(list_of_ws))

    first_line=True
    for i in range(len(list_of_ws)):
        w = int(list_of_ws[i])
        subset_runs = rec_df.loc[rec_df["w"] == w, :].sort_values(by="t", ascending=False).reset_index(drop=True)
        
        plt.plot(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), label="w = {} bp".format(w))

        if first_line:
            plt.fill_between(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), alpha=0.6)
        else:
            plt.fill_between(subset_runs.t, previous_line.s, subset_runs.s, color=cmap(gradients[i]), alpha=0.75)
            
        previous_line = subset_runs
        first_line=False

    plt.xlabel("Overlap threshold")
    plt.ylabel("Reproducibility")
    plt.xticks(rotation=45, fontsize=7)
    plt.title("Reproducibility -- General")
    plt.legend()
    plt.tight_layout()
    
    if os.path.exists(savedir+"/contours_OvrWind/") == False:
        os.mkdir(savedir+"/contours_OvrWind/")

    plt.savefig('{}/reprod_lines_general.pdf'.format(savedir+"/contours_OvrWind"), format='pdf')
    plt.savefig('{}/reprod_lines_general.svg'.format(savedir+"/contours_OvrWind"), format='svg')
    plt.clf()

    X, Y = np.meshgrid(x, y)

    Z = np.zeros(X.shape)
    
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            corresponding_line = rec_df.loc[(rec_df["w"] == X[i,j]) & (rec_df["t"]== Y[i,j]), :]
            Z[i, j] = float(corresponding_line["s"])

    plt.contourf(X, Y, Z, cmap='RdGy_r')
    plt.colorbar()

    plt.ylabel("Overlap threshold")

    plt.xlabel("Window size")
    plt.xticks(rotation=45, fontsize=7)
    plt.title("Reproducibility -- General")
    plt.tight_layout()
    plt.savefig('{}/reprod_contour_general.pdf'.format(savedir+"/contours_OvrWind"), format='pdf')
    plt.savefig('{}/reprod_contour_general.svg'.format(savedir+"/contours_OvrWind"), format='svg')
    plt.clf()

    for k in perlabel_rec.keys():
        rec = np.array(perlabel_rec[k])

        x, y, z = rec[:,0], rec[:,1], rec[:,2]
        rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])

        list_of_ws = np.unique(rec[:, 0])
        cmap = get_cmap('inferno')
        gradients = np.linspace(0, 1, len(list_of_ws))

        first_line = True
        for i in range(len(list_of_ws)):
            w = int(list_of_ws[i])
            subset_runs = rec_df.loc[rec_df["w"] == w, :].sort_values(by="t", ascending=False).reset_index(drop=True)
            
            plt.plot(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), label="w = {} bp".format(w))

            if first_line:
                plt.fill_between(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), alpha=0.6)
            else:
                plt.fill_between(subset_runs.t, previous_line.s, subset_runs.s, color=cmap(gradients[i]), alpha=0.75)
                
            previous_line = subset_runs
            first_line=False

        plt.xlabel("Overlap threshold")

        plt.ylabel("Reproducibility")
        plt.xticks(rotation=45, fontsize=7)
        plt.title("Reproducibility -- {}".format(k))
        plt.tight_layout()
        plt.legend()
        plt.savefig('{}/reprod_lines_{}.pdf'.format(savedir+"/contours_OvrWind", k), format='pdf')
        plt.savefig('{}/reprod_lines_{}.svg'.format(savedir+"/contours_OvrWind", k), format='svg')
        plt.clf()

        X, Y = np.meshgrid(x, y)

        Z = np.zeros(X.shape)
        
        for i in range(Z.shape[0]):
            for j in range(Z.shape[1]):
                corresponding_line = rec_df.loc[(rec_df["w"] == X[i,j]) & (rec_df["t"]== Y[i,j]), :]
                Z[i, j] = float(corresponding_line["s"])

        plt.contourf(X, Y, Z, cmap='RdGy_r')
        plt.colorbar()
        plt.ylabel("Overlap threshold")

        plt.xlabel("Window size")
        plt.xticks(rotation=45, fontsize=7)
        plt.title("Reproducibility -- {}".format(k))
        plt.tight_layout()
        plt.savefig('{}/reprod_contour_{}.pdf'.format(savedir+"/contours_OvrWind", k), format='pdf')
        plt.savefig('{}/reprod_contour_{}.svg'.format(savedir+"/contours_OvrWind", k), format='svg')
        plt.clf()

    ##################################################################################################
    
    num_labels = loci1.shape[1]-3
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[25, 16])
    label_being_plotted = 0
    plrk = list(perlabel_rec.keys())

    for i in range(n_rows):
        for j in range(n_cols):
            k = plrk[label_being_plotted]
            rec = np.array(perlabel_rec[k])

            x, y, z = rec[:,0], rec[:,1], rec[:,2]
            rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])

            list_of_ws = np.unique(rec[:, 0])
            cmap = get_cmap('inferno')
            gradients = np.linspace(0, 1, len(list_of_ws))

            first_line = True
            for c in range(len(list_of_ws)):
                w = int(list_of_ws[c])
                subset_runs = rec_df.loc[rec_df["w"] == w, :].sort_values(by="t", ascending=False).reset_index(drop=True)
                
                axs[i,j].plot(subset_runs.t, subset_runs.s, color=cmap(gradients[c]))#, label="w = {} bp".format(w))
                
                if first_line:
                    axs[i,j].fill_between(subset_runs.t, subset_runs.s, color=cmap(gradients[c]), alpha=0.6)
                else:
                    axs[i,j].fill_between(subset_runs.t, previous_line.s, subset_runs.s, color=cmap(gradients[c]), alpha=0.75)
                    
                previous_line = subset_runs
                first_line=False

            label_being_plotted+=1
            axs[i,j].set_title(k, fontsize=13)

    plt.tight_layout()
    plt.savefig('{}/reprod_contour_lines.pdf'.format(savedir+"/contours_OvrWind"), format='pdf')
    plt.savefig('{}/reprod_contour_lines.svg'.format(savedir+"/contours_OvrWind"), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')

def ReprThresWind_contour(
        loci1, loci2, savedir, w_range=[0, 4000, 200], t_range=[50, 100, 5], posterior=True, 
        matching="static", static_thres = 0.75):
    # t is the reproducibility threshold
    """
    initialize wrange and trange and reverse trange to start from the most strict to easiest
    initialize a matrix of size #t * #w and put a vector of [False for i in range len(loci1)] as each element

    for t in trange
        for w in wrange
            get_reproducibility of A[t, w] for all positions where rep is still False 
            for all reproducibile positions:
                change that position to True if that cell in A has bigger window or lower threshold
    
    after going through A once, return A[t,w] = len(A[t,w]==True) / len(loci1)
    """

    t_range = [i/100 for i in range(t_range[0], t_range[1], t_range[2])]
    t_range = t_range[::-1]

    w_range = list(range(w_range[0], w_range[1], w_range[2]))

    # initialize A
    A = {}
    for t in t_range:
        for w in w_range:
            if t in A.keys():
                A[t][w] = [False for i in range(len(loci1))]
            else:
                A[t] = {}
                A[t][w] = [False for i in range(len(loci1))]
    
    for t in t_range:
        for w in w_range:
            """
            get a subset of loci1 where A[t,w] == False
            run reproducibility_test on the subset of loci1
            update A[t,w] with new reprod_results
            """
            t0 = datetime.now()

            if posterior:
                updated_repr = is_repr_posterior(loci1, loci2, static_thres, w, reproducibility_threshold=t, matching=matching)

            else:
                updated_repr = is_repr_MAP_with_prior(A[t][w], loci1, loci2, static_thres, w, matching=matching)

            print(datetime.now() - t0)
            
            # now update all of A with new info
            for tt in t_range:
                for ww in w_range:
                    if tt <= t and ww >= w:
                        A[tt][ww] = updated_repr
    
    rec = []
    perlabel_rec = {k:[] for k in loci1.columns[3:]}
    for thresh in A.keys():
        for w in A[thresh].keys():
            reprod_report = pd.DataFrame(A[thresh][w], columns=["is_repr"])

            score = len(reprod_report.loc[reprod_report["is_repr"]==True])/ len(reprod_report)
            rec.append([w, thresh, score])

            lab_rep = perlabel_is_reproduced(reprod_report, loci1)
            for k, v in lab_rep.items():
                perlabel_rec[k].append([w, thresh, v[0]/ (v[0]+v[1])])
    
    rec = np.array(rec)
    x, y, z = rec[:,0], rec[:,1], rec[:,2]
    rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])
    list_of_ws = np.unique(rec[:, 0])
    cmap = get_cmap('inferno')
    gradients = np.linspace(0, 1, len(list_of_ws))

    first_line=True
    for i in range(len(list_of_ws)):
        w = int(list_of_ws[i])
        subset_runs = rec_df.loc[rec_df["w"] == w, :].sort_values(by="t", ascending=False).reset_index(drop=True)
        
        plt.plot(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), label="w = {} bp".format(w))

        if first_line:
            plt.fill_between(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), alpha=0.6)
        else:
            plt.fill_between(subset_runs.t, previous_line.s, subset_runs.s, color=cmap(gradients[i]), alpha=0.75)
            
        previous_line = subset_runs
        first_line=False

    plt.xlabel("Reproducibility score threshold")
    plt.ylabel("Reproducibility")
    plt.xticks(rotation=45, fontsize=7)
    plt.title("Reproducibility -- General")
    plt.legend()
    plt.tight_layout()
    
    if os.path.exists(savedir+"/contours_ReprThresWind/") == False:
        os.mkdir(savedir+"/contours_ReprThresWind/")

    plt.savefig('{}/reprod_lines_general.pdf'.format(savedir+"/contours_ReprThresWind"), format='pdf')
    plt.savefig('{}/reprod_lines_general.svg'.format(savedir+"/contours_ReprThresWind"), format='svg')
    plt.clf()

    X, Y = np.meshgrid(x, y)

    Z = np.zeros(X.shape)
    
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            corresponding_line = rec_df.loc[(rec_df["w"] == X[i,j]) & (rec_df["t"]== Y[i,j]), :]
            Z[i, j] = float(corresponding_line["s"])

    plt.contourf(X, Y, Z, cmap='RdGy_r')
    plt.colorbar()

    plt.ylabel("Reproducibility score threshold")

    plt.xlabel("Window size")
    plt.xticks(rotation=45, fontsize=7)
    plt.title("Reproducibility -- General")
    plt.tight_layout()
    plt.savefig('{}/reprod_contour_general.pdf'.format(savedir+"/contours_ReprThresWind"), format='pdf')
    plt.savefig('{}/reprod_contour_general.svg'.format(savedir+"/contours_ReprThresWind"), format='svg')
    plt.clf()

    for k in perlabel_rec.keys():
        rec = np.array(perlabel_rec[k])

        x, y, z = rec[:,0], rec[:,1], rec[:,2]
        rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])

        list_of_ws = np.unique(rec[:, 0])
        cmap = get_cmap('inferno')
        gradients = np.linspace(0, 1, len(list_of_ws))

        first_line = True
        for i in range(len(list_of_ws)):
            w = int(list_of_ws[i])
            subset_runs = rec_df.loc[rec_df["w"] == w, :].sort_values(by="t", ascending=False).reset_index(drop=True)
            
            plt.plot(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), label="w = {} bp".format(w))

            if first_line:
                plt.fill_between(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), alpha=0.6)
            else:
                plt.fill_between(subset_runs.t, previous_line.s, subset_runs.s, color=cmap(gradients[i]), alpha=0.75)
                
            previous_line = subset_runs
            first_line=False

        plt.xlabel("Reproducibility score threshold")

        plt.ylabel("Reproducibility")
        plt.xticks(rotation=45, fontsize=7)
        plt.title("Reproducibility -- {}".format(k))
        plt.tight_layout()
        plt.legend()
        plt.savefig('{}/reprod_lines_{}.pdf'.format(savedir+"/contours_ReprThresWind", k), format='pdf')
        plt.savefig('{}/reprod_lines_{}.svg'.format(savedir+"/contours_ReprThresWind", k), format='svg')
        plt.clf()

        X, Y = np.meshgrid(x, y)

        Z = np.zeros(X.shape)
        
        for i in range(Z.shape[0]):
            for j in range(Z.shape[1]):
                corresponding_line = rec_df.loc[(rec_df["w"] == X[i,j]) & (rec_df["t"]== Y[i,j]), :]
                Z[i, j] = float(corresponding_line["s"])

        plt.contourf(X, Y, Z, cmap='RdGy_r')
        plt.colorbar()
        plt.ylabel("Reproducibility score threshold")

        plt.xlabel("Window size")
        plt.xticks(rotation=45, fontsize=7)
        plt.title("Reproducibility -- {}".format(k))
        plt.tight_layout()
        plt.savefig('{}/reprod_contour_{}.pdf'.format(savedir+"/contours_ReprThresWind", k), format='pdf')
        plt.savefig('{}/reprod_contour_{}.svg'.format(savedir+"/contours_ReprThresWind", k), format='svg')
        plt.clf()
    
    ##################################################################################################
    
    num_labels = loci1.shape[1]-3
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=[25, 16])
    label_being_plotted = 0
    plrk = list(perlabel_rec.keys())

    for i in range(n_rows):
        for j in range(n_cols):
            k = plrk[label_being_plotted]
            rec = np.array(perlabel_rec[k])

            x, y, z = rec[:,0], rec[:,1], rec[:,2]
            rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])

            list_of_ws = np.unique(rec[:, 0])
            cmap = get_cmap('inferno')
            gradients = np.linspace(0, 1, len(list_of_ws))

            first_line = True
            for c in range(len(list_of_ws)):
                w = int(list_of_ws[c])
                subset_runs = rec_df.loc[rec_df["w"] == w, :].sort_values(by="t", ascending=False).reset_index(drop=True)
                
                axs[i,j].plot(subset_runs.t, subset_runs.s, color=cmap(gradients[c]))#, label="w = {} bp".format(w))
                
                if first_line:
                    axs[i,j].fill_between(subset_runs.t, subset_runs.s, color=cmap(gradients[c]), alpha=0.6)
                else:
                    axs[i,j].fill_between(subset_runs.t, previous_line.s, subset_runs.s, color=cmap(gradients[c]), alpha=0.75)
                    
                previous_line = subset_runs
                first_line=False

            label_being_plotted+=1
            axs[i,j].set_title(k, fontsize=13)

    plt.tight_layout()
    plt.savefig('{}/reprod_contour_lines.pdf'.format(savedir+"/contours_ReprThresWind"), format='pdf')
    plt.savefig('{}/reprod_contour_lines.svg'.format(savedir+"/contours_ReprThresWind"), format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')

def __OvrWind_delta_NMI_contour(loci1, loci2, savedir, w_range=[0, 4000, 200], t_range=[0, 11, 1], posterior=True, repr_threshold=0.9):
    # NMI is defined for whole annotations. it cannot be used for individal labels
    t_range = [i/10 for i in range(t_range[0], t_range[1], t_range[2])]
    t_range = t_range[::-1]

    w_range = list(range(w_range[0], w_range[1], w_range[2]))

    # initialize A
    A = {}
    for t in t_range:
        for w in w_range:
            if t in A.keys():
                A[t][w] = 0
            else:
                A[t] = {}
                A[t][w] = 0
    
    for t in t_range:
        for w in w_range:
            """
            get a subset of loci1 where A[t,w] == False
            run reproducibility_test on the subset of loci1
            update A[t,w] with new reprod_results
            """
            t0 = datetime.now()

            if posterior:
                delta_NMI = single_point_delta_NMI_repr(
                    loci1, loci2, ovr_threshold=t, window_bp=w, repr_threshold= repr_threshold, posterior=True, matching="static")

            else:
                delta_NMI = single_point_delta_NMI_repr(
                    loci1, loci2, ovr_threshold=t, window_bp=w, repr_threshold= repr_threshold, posterior=False, matching="static")

            print(datetime.now() - t0)
            A[t][w] = delta_NMI

    rec = []
    for thresh in A.keys():
        for w in A[thresh].keys():

            deltaNMI_score = A[thresh][w]
            rec.append([w, thresh, deltaNMI_score])

    rec = np.array(rec)
    x, y, z = rec[:,0], rec[:,1], rec[:,2]
    rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])
    list_of_ws = np.unique(rec[:, 0])
    list_of_ws = np.sort(list_of_ws)[::-1]
    cmap = get_cmap('inferno')
    gradients = np.linspace(0, 1, len(list_of_ws))

    first_line=True
    for i in range(len(list_of_ws)):
        w = int(list_of_ws[i])
        subset_runs = rec_df.loc[rec_df["w"] == w, :].sort_values(by="t", ascending=False).reset_index(drop=True)
        
        plt.plot(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), label="w = {} bp".format(w))

        if first_line:
            plt.fill_between(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), alpha=0.6)
        else:
            plt.fill_between(subset_runs.t, previous_line.s, subset_runs.s, color=cmap(gradients[i]), alpha=0.75)
            
        previous_line = subset_runs
        first_line=False

    plt.xlabel("Overlap threshold")
    plt.ylabel("Delta_NMI")
    plt.xticks(rotation=45, fontsize=7)
    plt.title("Delta_NMI -- General")
    plt.legend()
    plt.tight_layout()
    
    if os.path.exists(savedir+"/contours_OvrWind_deltaNMI/") == False:
        os.mkdir(savedir+"/contours_OvrWind_deltaNMI/")

    plt.savefig('{}/deltaNMI_lines_general.pdf'.format(savedir+"/contours_OvrWind_deltaNMI"), format='pdf')
    plt.savefig('{}/deltaNMI_lines_general.svg'.format(savedir+"/contours_OvrWind_deltaNMI"), format='svg')
    plt.clf()

    X, Y = np.meshgrid(x, y)

    Z = np.zeros(X.shape)
    
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            corresponding_line = rec_df.loc[(rec_df["w"] == X[i,j]) & (rec_df["t"]== Y[i,j]), :]
            Z[i, j] = float(corresponding_line["s"])

    plt.contourf(X, Y, Z, cmap='RdGy_r')
    plt.colorbar()

    plt.ylabel("Overlap threshold")

    plt.xlabel("Window size")
    plt.xticks(rotation=45, fontsize=7)
    plt.title("Delta_NMI -- General")
    plt.tight_layout()
    plt.savefig('{}/deltaNMI_contour_general.pdf'.format(savedir+"/contours_OvrWind_deltaNMI"), format='pdf')
    plt.savefig('{}/deltaNMI_contour_general.svg'.format(savedir+"/contours_OvrWind_deltaNMI"), format='svg')
    plt.clf()

def __ReprThresWind_delta_NMI_contour(
        loci1, loci2, savedir, w_range=[0, 4000, 200], 
        t_range=[50, 100, 5], posterior=True, matching="static", static_thres=0.75):
    # NMI is defined for whole annotations. it cannot be used for individal labels
    t_range = [i/100 for i in range(t_range[0], t_range[1], t_range[2])]
    t_range = t_range[::-1]

    w_range = list(range(w_range[0], w_range[1], w_range[2]))

    # initialize A
    A = {}
    for t in t_range:
        for w in w_range:
            if t in A.keys():
                A[t][w] = 0
            else:
                A[t] = {}
                A[t][w] = 0
    
    for t in t_range:
        for w in w_range:
            """
            get a subset of loci1 where A[t,w] == False
            run reproducibility_test on the subset of loci1
            update A[t,w] with new reprod_results
            """
            t0 = datetime.now()

            if posterior:
                delta_NMI = single_point_delta_NMI_repr(
                    loci1, loci2, ovr_threshold=static_thres, window_bp=w, repr_threshold=t, posterior=True, matching=matching)

            else:
                delta_NMI = single_point_delta_NMI_repr(
                    loci1, loci2, ovr_threshold=static_thres, window_bp=w, repr_threshold=t, posterior=False, matching=matching)

            print(datetime.now() - t0)
            A[t][w] = delta_NMI

    rec = []
    for thresh in A.keys():
        for w in A[thresh].keys():

            deltaNMI_score = A[thresh][w]
            rec.append([w, thresh, deltaNMI_score])

    rec = np.array(rec)
    x, y, z = rec[:,0], rec[:,1], rec[:,2]
    rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])
    list_of_ws = np.unique(rec[:, 0])
    list_of_ws = np.sort(list_of_ws)[::-1]
    cmap = get_cmap('inferno')
    gradients = np.linspace(0, 1, len(list_of_ws))

    first_line=True
    for i in range(len(list_of_ws)):
        w = int(list_of_ws[i])
        subset_runs = rec_df.loc[rec_df["w"] == w, :].sort_values(by="t", ascending=False).reset_index(drop=True)
        
        plt.plot(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), label="w = {} bp".format(w))

        if first_line:
            plt.fill_between(subset_runs.t, subset_runs.s, color=cmap(gradients[i]), alpha=0.6)
        else:
            plt.fill_between(subset_runs.t, previous_line.s, subset_runs.s, color=cmap(gradients[i]), alpha=0.75)
            
        previous_line = subset_runs
        first_line=False

    plt.xlabel("Reproducibility score threshold")
    plt.ylabel("Delta_NMI")
    plt.xticks(rotation=45, fontsize=7)
    plt.title("Delta_NMI -- General")
    plt.legend()
    plt.tight_layout()
    
    if os.path.exists(savedir+"/contours_ReprThresWind_deltaNMI/") == False:
        os.mkdir(savedir+"/contours_ReprThresWind_deltaNMI/")

    plt.savefig('{}/deltaNMI_lines_general.pdf'.format(savedir+"/contours_ReprThresWind_deltaNMI"), format='pdf')
    plt.savefig('{}/deltaNMI_lines_general.svg'.format(savedir+"/contours_ReprThresWind_deltaNMI"), format='svg')
    plt.clf()

    X, Y = np.meshgrid(x, y)

    Z = np.zeros(X.shape)
    
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            corresponding_line = rec_df.loc[(rec_df["w"] == X[i,j]) & (rec_df["t"]== Y[i,j]), :]
            Z[i, j] = float(corresponding_line["s"])

    plt.contourf(X, Y, Z, cmap='RdGy_r')
    plt.colorbar()

    plt.ylabel("Reproducibility score threshold")

    plt.xlabel("Window size")
    plt.xticks(rotation=45, fontsize=7)
    plt.title("Delta_NMI -- General")
    plt.tight_layout()
    plt.savefig('{}/deltaNMI_contour_general.pdf'.format(savedir+"/contours_ReprThresWind_deltaNMI"), format='pdf')
    plt.savefig('{}/deltaNMI_contour_general.svg'.format(savedir+"/contours_ReprThresWind_deltaNMI"), format='svg')
    plt.clf()
