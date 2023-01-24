from _cluster_matching import *
import math, scipy
import os
from _reproducibility import *


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
        enr_ovr = enrichment_of_overlap_matrix(loci_1, loci_2, OE_transform=False)

    else:
        enr_ovr = enrichment_of_overlap_matrix(loci_1, loci_2, OE_transform=True)

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

def calibrate(vector_1, MAP2, k, window_bin, strat_size="def", num_bins=500, allow_w=True):
    """
    sort based on v1
    for t in range(len(v1)):
        if map2[t] 
    """

    if strat_size == "def":
        strat_size = int(len(vector_1)/num_bins)

    coverage_2 = len(MAP2.loc[MAP2 == k]) / len(MAP2)
    vector_pair = pd.concat([vector_1, MAP2], axis=1)

    vector_pair.columns=["pv1", "map2"]
    MAP2 = list(MAP2)
    
    if allow_w:
        vector_pair = vector_pair.sort_values("pv1").reset_index()

    else:
        vector_pair = vector_pair.sort_values("pv1").reset_index(drop=True)
    
    bins = []
    for b in range(0, vector_pair.shape[0], strat_size):
        # [bin_start, bin_end, num_values_in_bin, num_agreement, num_mislabeled, ratio_agreement]

        if allow_w:
            subset_vector_pair = vector_pair.iloc[b:b+strat_size, :]
            subset_vector_pair = subset_vector_pair.reset_index(drop=True)
            
            observed = 0
            for i in range(len(subset_vector_pair)):
      
                leftmost = max(0, subset_vector_pair["index"][i] - window_bin)
                rightmost = min(len(MAP2), subset_vector_pair["index"][i] + window_bin)

                neighbors_i = MAP2[leftmost:rightmost]

                if k in neighbors_i:
                    observed += 1

        else:
            subset_vector_pair = vector_pair.iloc[b:b+strat_size, :]
            observed = len(subset_vector_pair.loc[subset_vector_pair["map2"] == k])

        expected = (coverage_2 * len(subset_vector_pair))

        if observed!=0:
            # observed += len(subset_vector_pair)*1e-3
            # expected += len(subset_vector_pair)*1e-3

            oe = np.log((observed)/(expected))
            
            bins.append([
                float(subset_vector_pair["pv1"][subset_vector_pair.index[0]]), 
                float(subset_vector_pair["pv1"][subset_vector_pair.index[-1]]),
                len(subset_vector_pair),
                observed,
                len(subset_vector_pair) - observed,
                oe])

    bins = np.array(bins)

    # plt.scatter(
    #     x=np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), 
    #     y=bins[:,5], c='black', s=10000/len(bins))
    
    polyreg = IsotonicRegression(
            y_min=float(np.array(bins[:, 5]).min()), y_max=float(np.array(bins[:, 5]).max()), 
            out_of_bounds="clip", increasing="auto")
    
    polyreg.fit(
        np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)), 
        bins[:, 5])
    
    # plt.plot(
    #     np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), 
    #     polyreg.predict(np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1))), 
    #     '--', c='r', linewidth=3)

    # plt.show()
    return polyreg
        
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

    # enr_ovr = enrichment_of_overlap_matrix(loci_1, loci_2, OE_transform=True)
    MP1 = loci_1.iloc[:,3:].max(axis=1) #posterior values 1
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1) # MAP label 1
    
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
    

    calibration_funcs = {}
    for i in range(len(loci_1.columns[3:])):
        if loci_1.columns[3:][i] not in calibration_funcs.keys():
            calibration_funcs[loci_1.columns[3:][i]] = {}

        for j in range(len(loci_2.columns[3:])):
            if loci_2.columns[3:][j] not in calibration_funcs[loci_1.columns[3:][i]].keys():

                # print(loci_1.columns[3:][i], loci_2.columns[3:][j])

                calibration_funcs[loci_1.columns[3:][i]][loci_2.columns[3:][j]] = \
                    calibrate(
                        loci_1.loc[:, loci_1.columns[3:][i]], MAP2,
                        k=loci_2.columns[3:][j], num_bins=numbins, allow_w=allow_w, window_bin=window_bin)

    posterior_enrichment = []
    # coverage_1 = {k:float(len(MAP1.loc[MAP1 == k]) / len(MAP2)) for k in MAP1.unique()}
    coverage_2 = {k:float(len(MAP2.loc[MAP2 == k]) / len(MAP2)) for k in MAP2.unique()}


    for t in range(len(loci_1)): 
        enr_val =calibration_funcs[MAP1[t]][MAP2[t]].predict(np.reshape(np.array(MP1[t]), (1,1)))[0]

        if raw_overlap_axis:
            batchsize = int(len(loci_1)/numbins)
            e = coverage_2[MAP2[t]] * batchsize

            overlap_ratio = (np.exp(enr_val)*e) / batchsize
            
            posterior_enrichment.append(overlap_ratio)

        else:
            posterior_enrichment.append(enr_val)
            
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

    enr_ovr = enrichment_of_overlap_matrix(loci_1, loci_2, OE_transform=True)

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

    for w in range(100, 3100, 600):
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
