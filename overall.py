from _cluster_matching import *
import math, scipy

def is_reproduced(loci_1, loci_2, enr_threshold=3, window_bp=500):
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
        
        for w in range(1, window_bin+1):
            if (t - w) >= 0 and (t + w) <= len(MAP1)-1 :

                if MAP2[t-w] == best_match[MAP1[t]] or \
                    MAP2[t-w] in above_threshold_match[MAP1[t]] or \
                        MAP2[t+w] == best_match[MAP1[t]] or \
                            MAP2[t+w] in above_threshold_match[MAP1[t]]:

                                statement3 = True


        if statement1 or statement2 or statement3:
            reprod_report[t][3] = True
        else:
            reprod_report[t][3] = False

    # print("reproduced {}".format(float(r/len(MAP1))))
    reprod_report = pd.DataFrame(reprod_report, columns=["chr", "start", "end", "is_repr"])

    # print(len(reprod_report.loc[reprod_report["is_repr"]==True])/ len(reprod_report))
    return reprod_report
            

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

def contour_isrep(loci1, loci2, savedir):
    rec = []
    for w in range(100, 3000, 100):
        for t in range(0, 50, 5):
            thresh = float(t)/10
            reprod_report = is_reproduced(loci1, loci2, enr_threshold=thresh, window_bp=w)

            score = len(reprod_report.loc[reprod_report["is_repr"]==True])/ len(reprod_report)

            rec.append([w, thresh, score])

    rec = np.array(rec)
    x, y, z = rec[:,0], rec[:,1], rec[:,2]
    X, Y = np.meshgrid(x, y)

    rec_df = pd.DataFrame(rec, columns=["w", "t", "s"])

    Z = np.zeros(X.shape)
    
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            corresponding_line = rec_df.loc[(rec_df["w"] == X[i,j]) & (rec_df["t"]== Y[i,j]), :]
            Z[i, j] = float(corresponding_line["s"])

    plt.contourf(X, Y, Z, cmap='RdGy_r')
    plt.colorbar()

    plt.ylabel("Enrichment of Overlap threshold")
    plt.xlabel("Window")
    plt.xticks(rotation=45, fontsize=7)
    plt.tight_layout()

    plt.savefig('{}/reprod_contour.pdf'.format(savedir), format='pdf')
    plt.savefig('{}/reprod_contour.svg'.format(savedir), format='svg')
    plt.clf()


# if __name__=="__main__":
#     replicate_1_dir = "tests/cedar_runs/segway/GM12878_R1/"
#     replicate_2_dir = "tests/cedar_runs/segway/GM12878_R2/"

#     loci1, loci2 = load_data(
#         replicate_1_dir+"/parsed_posterior.csv",
#         replicate_2_dir+"/parsed_posterior.csv",
#         subset=True, logit_transform=True)

#     loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

#     contour_isrep(loci1, loci2)

#     reprep = reprod_score(loci1, loci2, window_bp=500)
#     max_posterior_to_repr_score(loci1, loci2, reprep)

    # is_reproduced(loci1, loci2, enr_threshold=2, window_bp=500)

