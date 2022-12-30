from reports import *
import math

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
            

def reprod_score_posterior_agnostic(loci_1, loci_2):
    """
    get enr_ovr matrix
    get MAPs

    for t in map1:
        in a window of size w 
            BMS = what is the best match score (BMS) -> max of enr_ovr[R2] within that window
            d_BMS = what is the distance to BMS

            Score_t = BMS/d_BMS 
    """
    # can we do the exact same thing but instead of MAPs, on all of logit(p) values
    # we know that for any label in R1, there is score (enr_ovr) connecting it to any label in R2
    # so for each position i and each clusterID j across R1 and R2, we can get the following:
        # score_ij = p_1_ij * p_2_ij * enr_ovr[R1j, R2j]  
    # this score can act as the BMS
    # but, what about the distance parameter???


# now, what if we form a mapping from logit(p) to reprod_score_posterior_agnostic ?
# mapping = {k:[] for k in labels}
# for t in loci:
#     if max(loci[t]) == k:
#         mapping[k].append(tuple(loci[t], score_t))


if __name__=="__main__":
    replicate_1_dir = "tests/cedar_runs/segway/MCF7_R1/"
    replicate_2_dir = "tests/cedar_runs/segway/MCF7_R2/"

    loci1, loci2 = load_data(
        replicate_1_dir+"/parsed_posterior.csv",
        replicate_2_dir+"/parsed_posterior.csv",
        subset=True, logit_transform=False)

    loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    # is_reproduced(loci1, loci2, enr_threshold=10)
    # is_reproduced(loci1, loci2, enr_threshold=5)
    is_reproduced(loci1, loci2, enr_threshold=2, window_bp=500)
    is_reproduced(loci1, loci2, enr_threshold=2, window_bp=1000)
    is_reproduced(loci1, loci2, enr_threshold=2, window_bp=2000)
    # is_reproduced(loci1, loci2, enr_threshold=1)
    # is_reproduced(loci1, loci2, enr_threshold=0)

