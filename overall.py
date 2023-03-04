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
        enr_ovr = overlap_matrix(loci_1, loci_2, type="IoU")
    else:
        enr_ovr = overlap_matrix(loci_1, loci_2, type="IoU")

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

    enr_ovr = overlap_matrix(loci_1, loci_2, type="IoU")
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

    enr_ovr = overlap_matrix(loci_1, loci_2, type="IoU")

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

def define_matches(loci1, loci2):
    """
    first we get the enrichment of overlap matrix 

    using the overlap matrix, for each R2 label, sort R1 labels 
    as done for granularity matrix, start merging R1 labels recursively 
    
    at each merging step, plot the ck vs genome coverage
    find max ck in the plot
    
    return list of R1 labels that need to be merged in order to get max ck
    then for each R1 label look at the previously calculated correspondences of R2 and find best ones
    """
    # enr_ovr = enrichment_of_overlap_matrix(loci1, loci2, OE_transform=True)
    enr_ovr = Cohens_Kappa_matrix(loci1, loci2)
    MAP1 = loci1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci2.iloc[:,3:].idxmax(axis=1)

    corresp = {}

    for l in loci2.columns[3:]:
        loci_1 = loci1.copy()

        coverage_record = []
        ck_record = []
        progression = []

        sorted_l_vector = enr_ovr.loc[:, l].sort_values(ascending=False)
        query_label = str(sorted_l_vector.index[0])

        for j in range(0, len(sorted_l_vector)):
            if j>0:
                loci_1[query_label + "+" + str(sorted_l_vector.index[j])] = \
                    loci_1[query_label] + loci_1[str(sorted_l_vector.index[j])]

                loci_1 = loci_1.drop([query_label, str(sorted_l_vector.index[j])], axis=1)

                query_label = query_label + "+" + str(sorted_l_vector.index[j])

            MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
            coverage_query = len(MAP1.loc[MAP1 == query_label])/len(loci_1)
            coverage_record.append(coverage_query)
            progression.append(query_label)

            ck =  cohen_kappa_score(
                list(MAP1 == query_label),
                list(MAP2 == l)
            )
            ck_record.append(ck)

        corresp[l] = [
            progression[ck_record.index(max(ck_record))], 
            ck_record[ck_record.index(max(ck_record))], 
            coverage_record[ck_record.index(max(ck_record))]] # best_set_of_labels, best_CK, best_coverage
        
        if "+" in corresp[l][0]:
            corresp[l][0] = corresp[l][0].split("+")
        else:
            corresp[l][0] = [corresp[l][0]]

        # plt.plot(coverage_record, ck_record)
        # plt.title(l + " | maxCK in : " + str(ck_record.index(max(ck_record))))
        # plt.xlabel("coverage")
        # plt.ylabel("Cohen's Kappa")
        # plt.show()
    print(corresp)
    return corresp   

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

def is_repr_MAP_with_prior(prior, loci_1, loci_2, enr_threshold, window_bp):
    resolution = loci_1["end"][0] - loci_1["start"][0]
    window_bin = math.ceil(window_bp/resolution)

    enr_ovr = overlap_matrix(loci_1, loci_2, type="IoU")

    best_match = {i:enr_ovr.loc[i, :].idxmax() for i in enr_ovr.index} 
    above_threshold_match = {i: list(enr_ovr.loc[i, enr_ovr.loc[i, :]>enr_threshold].index) for i in enr_ovr.index}

    for k in above_threshold_match.keys():
        if len(above_threshold_match[k]) == 0:
            above_threshold_match[k] = [best_match[k]]

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

def is_repr_posterior(loci_1, loci_2, enr_threshold, window_bp, numbins=1000):
    # based on the resolution determine the number of bins for w
    resolution = loci_1["end"][0] - loci_1["start"][0] 
    window_bin = math.ceil(window_bp/resolution)

    # get the overlap ratio
    enr_ovr = overlap_matrix(loci_1, loci_2, type="IoU")
    MAP1 = list(loci_1.iloc[:,3:].idxmax(axis=1)) 
    MAP2 = list(loci_2.iloc[:,3:].idxmax(axis=1))

    # determine the number of positions at each posterior bin
    strat_size = int(len(loci_1)/numbins)

    ################# GETTING CALIBRATIONS #################
    calibrations = {}

    for k in loci_1.columns[3:]:

        ####################### DEFINE MATCHES #######################
        best_match = [enr_ovr.loc[k, :].idxmax()]
        above_threshold_match = list(enr_ovr.loc[k, enr_ovr.loc[k, :]>enr_threshold].index)

        if len(above_threshold_match) == 0:
            above_threshold_match = [best_match[0]]
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

    # here I'm just looking at the calibrated score that is assigned to the MAP1 of each position (~>.9)
    calibrated_reproducibility = []
    for i in range(len(MAP1)):
        calibrated_reproducibility.append(calibrated_loci1.loc[i, MAP1[i]])

    # then based on the initial threshold that we had, we say a position is reproduced if it has a calibrated score of >threshold
    binary_isrep = []
    for t in calibrated_reproducibility:
        if t>=enr_threshold:
            binary_isrep.append(True)
        else:
            binary_isrep.append(False)

    return binary_isrep

def single_point_repr(loci_1, loci_2, enr_threshold, window_bp, posterior=False):
    if posterior:
        """
        for all labels in R1:
            creat a calibration curve according to t, w
        
        according to the calibrated p, 
        """
        return is_repr_posterior(loci_1, loci_2, enr_threshold, window_bp, numbins=1000)

    else:
        return is_repr_MAP_with_prior(
            [False for _ in range(len(loci_1))], 
            loci_1, loci_2, enr_threshold, window_bp)

def fast_contour(loci1, loci2, savedir, w_range=[0, 4000, 200], t_range=[0, 11, 1], posterior=True, raw_overlap_axis=True):
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
                updated_repr = is_repr_posterior(loci1, loci2, t, w, numbins=1000)

            else:
                updated_repr = is_repr_MAP_with_prior(A[t][w], loci1, loci2, t, w, raw_overlap_axis=True)

            print(datetime.now() - t0)
            
            # now update all of A with new info
            for tt in t_range:
                for ww in w_range:
                    if tt <= t and ww >= w:
                        A[tt][ww] = updated_repr
            
            # all_scores = pd.DataFrame(np.zeros((len(t_range), len(w_range))), index=t_range, columns=w_range)
            # for ttt in t_range:
            #     for www in w_range:
            #         all_scores.loc[ttt, www] = sum(A[ttt][www]) / len(A[ttt][www])
            # print(all_scores)
    
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

