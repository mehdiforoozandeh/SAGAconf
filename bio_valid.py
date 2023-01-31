from _reproducibility import *
# from reports import *
from sklearn.linear_model import LinearRegression 
import scipy
import os
import pickle

def load_gene_coords(file, drop_negative_strand=True, drop_overlapping=True):
    gene_coords = pd.read_csv(file)
    gene_coords = gene_coords.drop(["Unnamed: 0"], axis=1)

    gene_coords["start"] = gene_coords["start"].astype("int")
    gene_coords["end"] = gene_coords["end"].astype("int")

    if drop_negative_strand:
        gene_coords = gene_coords.loc[gene_coords["strand"]=="+", :].reset_index(drop=True)
    
    if drop_overlapping:
        todrop = []
        for i in range(len(gene_coords)-1):
            if get_overlap((gene_coords["start"][i], gene_coords["end"][i]),(gene_coords["start"][i+1], gene_coords["end"][i+1])) >0:
                if (gene_coords["end"][i] - gene_coords["start"][i]) <= gene_coords["end"][i+1] - gene_coords["start"][i+1]:
                    todrop.append(i)
                else:
                    todrop.append(i+1)
        gene_coords = gene_coords.drop(todrop).reset_index(drop=True)

    return gene_coords

def load_TSS(TSSdir):
    tss = pd.read_csv(TSSdir, sep="\t", header=None)
    tss.columns = ["chr", "coord", "strand"]
    return tss

def get_coverage(loci):
    MAP = loci.iloc[:,3:].idxmax(axis=1)
    coverage = dict(zip(list(loci.columns[3:]), [0 for _ in range(len(loci.columns[3:]))]))

    for c in list(loci.columns[3:]):
        coverage[c] = len(MAP.loc[MAP == c]) / len(MAP)

    return coverage

def load_transcription_data(file, gene_coord, csv=True):
    if "csv" in file:
        csv==True
        
    elif "pkl" in file:
        csv==False

    if csv:
        trn_data = pd.read_csv(file, sep="\t")

        for j in range(len(trn_data)):
            trn_data.at[j, "gene_id"] = trn_data["gene_id"][j].split(".")[0]
        
        for i in range(len(gene_coord)):
            gene_coord.at[i, "gene_id"] = gene_coord["gene_id"][i].split(".")[0]

        mapped_trn_data = []
        for i in range(len(gene_coord)):
            geneID = gene_coord["gene_id"][i]
            subset = trn_data.loc[trn_data["gene_id"] == geneID, :].reset_index(drop=True)

            if len(subset) > 0:
                mapped_trn_data.append([
                    gene_coord["chr"][i], gene_coord["start"][i], gene_coord["end"][i], geneID, subset["length"][0], subset["TPM"][0], subset["FPKM"][0]
                ])

        mapped_trn_data = pd.DataFrame(mapped_trn_data, columns=["chr", "start", "end", "geneID", "length", "TPM", "FPKM"])
        return mapped_trn_data

    else:
        with open(file, "rb") as f:
            trn_data = pickle.load(f)

        trn_data = pd.DataFrame([[k, v] for k,v in trn_data.items()], columns=["gene_id", "TPM"])

        for j in range(len(trn_data)):
            trn_data.at[j, "gene_id"] = trn_data["gene_id"][j].split(".")[0]

        for i in range(len(gene_coord)):
            gene_coord.at[i, "gene_id"] = gene_coord["gene_id"][i].split(".")[0]
        
        mapped_trn_data = []
        for i in range(len(gene_coord)):
            geneID = gene_coord["gene_id"][i]
            subset = trn_data.loc[trn_data["gene_id"] == geneID, :].reset_index(drop=True)

            if len(subset) > 0:
                mapped_trn_data.append([
                    gene_coord["chr"][i], gene_coord["start"][i], gene_coord["end"][i], geneID, "-", 10**(subset["TPM"][0]), "-"])
        
        mapped_trn_data = pd.DataFrame(mapped_trn_data, columns=["chr", "start", "end", "geneID", "length", "TPM", "FPKM"])
        return mapped_trn_data

def intersect(loci, gene_coord, add_expression=True):
    cols = list(loci.columns)
    
    loci = loci.values.tolist()

    intersect = []
    for i in range(len(loci)):

        gene_coord_chr = gene_coord.loc[gene_coord["chr"] == loci[i][0], :]
        gene_coord_chr = gene_coord_chr.values.tolist()
        o = False
        tpm = 0
        for j in range(len(gene_coord_chr)):
            statement = bool(
                    bool(gene_coord_chr[j][1] <= loci[i][1] <= gene_coord_chr[j][2]) or
                    bool(gene_coord_chr[j][1] <= loci[i][2] <= gene_coord_chr[j][2]))

            if statement:
                o = True
                if add_expression:
                    tpm = float(gene_coord_chr[j][5])
                continue
        
        if o:
            if add_expression:
                intersect.append([tpm] + loci[i])
            else:
                intersect.append(loci[i])
    
    if add_expression:
        intersect = pd.DataFrame(intersect, columns=["TPM"]+cols)
    else:
        pd.DataFrame(intersect, columns=cols)
    return intersect

def condense_MAP(coordMAP):
    coordMAP = coordMAP.values.tolist()
    i=0
    while i+1 < len(coordMAP):
        if coordMAP[i][2] == coordMAP[i+1][1]:
            coordMAP[i][2] = coordMAP[i+1][2]
            del coordMAP[i+1]

        else:
            i += 1
    
    return pd.DataFrame(coordMAP, columns=["chr", "start", "end", "MAP"])

def condense_posterior(coord_posterior):
    cols = coord_posterior.columns
    coord_posterior = coord_posterior.values.tolist()
    i=0
    while i+1 < len(coord_posterior):
        statement = bool(
            coord_posterior[i][0]==coord_posterior[i+1][0] and\
                 coord_posterior[i][2] == coord_posterior[i+1][1] and\
                     coord_posterior[i][3]==coord_posterior[i+1][3])
        if statement:
            coord_posterior[i][2] = coord_posterior[i+1][2]
            del coord_posterior[i+1]

        else:
            i += 1
    
    return pd.DataFrame(coord_posterior, columns=cols)

def get_overlap(tup1, tup2):

    x = range(tup1[0], tup1[1])
    y = range(tup2[0], tup2[1])

    return len( range(max(x[0], y[0]), min(x[-1], y[-1])+1))        

def general_transcribed_enrich(loci, k, gene_coords):
    
    """
    label k
        e = what percentage of k overlaps with known gene-body regions (+- 5kb)
        c = coverage of label k
        return e/c
    """

    if "gene_type" not in gene_coords.columns:
        gene_coords["gene_type"] = np.array(["protein_coding" for _ in range(len(gene_coords))])

    MAP = pd.concat(
        [loci["chr"], loci["start"], loci["end"], 
        loci.iloc[:,3:].idxmax(axis=1)],
        axis=1)
    MAP.columns = ["chr", "start", "end", "MAP"]

    MAP_k = MAP.loc[MAP["MAP"]==k].reset_index(drop=True)

    denseMAP_k = condense_MAP(MAP_k)

    report = []
    total_segment_len = sum([denseMAP_k["end"][i] - denseMAP_k["start"][i] for i in range(denseMAP_k.shape[0])])

    for i in range(denseMAP_k.shape[0]):
        gene_coords_subset_chr = gene_coords.loc[
            gene_coords["chr"] == denseMAP_k["chr"][i]].reset_index(drop=True)

        for j in range(len(gene_coords_subset_chr)):
            o = get_overlap(
                (gene_coords_subset_chr["start"][j], gene_coords_subset_chr["end"][j]), 
                (denseMAP_k["start"][i], denseMAP_k["end"][i])
            )

            if o>0:
                report.append([o, gene_coords_subset_chr["gene_type"][j]])
    
    report = pd.DataFrame(report, columns=["overlap", "geneType"])
    genetypes = list(report["geneType"].unique())

    new_report = {g:0 for g in genetypes}
    for i in range(len(report)):
        new_report[report["geneType"][i]] +=  float(report["overlap"][i])/total_segment_len

    return new_report

def plot_general_transc(loci, gene_coords, savedir, exp=True):
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)

    transc_reports = {}
    for k in loci.columns[3:]:
        trans_enr_k = general_transcribed_enrich(loci, k, gene_coords)
        transc_reports[k] = trans_enr_k
        
    transc_reports = pd.DataFrame(transc_reports).T
    transc_reports = transc_reports.reindex(sorted(transc_reports.columns), axis=1)

    ax = transc_reports.plot.bar()
    plt.xlabel("labels")
    plt.xticks(rotation=45)
    plt.ylabel("ratio overlap with gene-body")

    # Put a legend to the right of the current axis
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25),
    #     ncol=4, fancybox=True, shadow=True, fontsize=9)

    plt.tight_layout()

    if exp:
        plt.savefig(savedir+"/general_transc_expressed.pdf", format='pdf')
        plt.savefig(savedir+"/general_transc_expressed.svg", format='svg')
    else:
        plt.savefig(savedir+"/general_transc_not_expressed.pdf", format='pdf')
        plt.savefig(savedir+"/general_transc_not_expressed.svg", format='svg')
    
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')

def posterior_transcribed_enrichment(loci, gene_coords, savedir, num_bins="def"):
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)
    """
    define bins of size B for posterior

    for posteriors in range B
        for k in labels

            what percentage of positions in that range overlap with protein-coding genes
                t = number of intersected regions with posterior in range B
                e = number of all regions with posterior in range B
    """
    cv = get_coverage(loci)

    # print(loci)
    # print(gene_coords)
    print("intersecting")
    intersection = intersect(loci, gene_coords)

    reports = {}
    for k in loci.columns[3:]:
        if num_bins=="dynamic":
            strat_size = int(len(loci)/(cv[k]*1000))
        
        else:
            strat_size = int(len(loci)/40)
        
        
        reports[k] = []
        posterior_vector = loci.loc[:, ["chr","start","end",k]]
        posterior_vector = posterior_vector.sort_values(by=k).reset_index(drop=True)

        for b in range(0, posterior_vector.shape[0], strat_size):
            subset_vector_pair = posterior_vector.iloc[b:b+strat_size, :].reset_index(drop=True)
            bin_range = [
                float(subset_vector_pair.iloc[0, 3]), 
                float(subset_vector_pair.iloc[-1, 3])
                ]
            
            t = len(intersection.loc[
                (bin_range[0]<=intersection[k])&
                (intersection[k]<=bin_range[1])  ,:])

            e = len(posterior_vector.loc[
                (bin_range[0]<=posterior_vector[k])&
                (posterior_vector[k]<=bin_range[1])  ,:])

            epsilon = 1e-3

            if t != 0:
                epsilon = 0
                obs = float(t) / len(intersection)
                exp = float(e) / len(posterior_vector)
                reports[k].append([bin_range, np.log(float((obs+epsilon)/(exp+epsilon)))])
    
    num_labels = len(loci.columns[3:])
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(10,8))
    colors = [i for i in get_cmap('tab20').colors]
    label_being_plotted = 0
    
    for i in range(n_rows):
        for j in range(n_cols):
            k = loci.columns[3:][label_being_plotted]

            x = [(kk[0][0]+kk[0][1])/2 for kk in reports[k]]
            y = [kk[1] for kk in reports[k]]

            axs[i,j].bar(
                x, y,  
                width=(np.max(loci[k].values)-np.min(loci[k].values))/len(reports[k]),
                color=colors[label_being_plotted],
                alpha=0.5)
            
            try:
                # if len(x)<=3:
                #     f = UnivariateSpline(x= np.array(x), y=y, k=1)
                # else:
                f = UnivariateSpline(x= np.array(x), y=y, k=1)

                axs[i,j].plot(
                    x, 
                    f(x),
                    label=k, c=colors[label_being_plotted], 
                    linewidth=1)
            except:
                pass

            axs[i,j].set_title(k, fontsize=7)
            label_being_plotted += 1

    plt.tight_layout()

    plt.savefig(savedir+"/transc_enrich.pdf", format='pdf')
    plt.savefig(savedir+"/transc_enrich.svg", format='svg')
    
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')

def posterior_transcription_correlation(loci, trans_data, savedir):
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)
    """
    get the intersection of trans_data and loci

    for k in labels:
        for i in intersected_loci
            for j in trans_data

                what posterior corresponds to what TPM?
    """
    intersection = intersect(loci, trans_data)

    report = {}
    for k in loci.columns[3:]:
        report[k] = []
        for i in range(len(intersection)):
            subset_trans_chr = trans_data.loc[(trans_data["chr"] == loci["chr"][i]),:].reset_index(drop=True)

            for j in range(len(subset_trans_chr)):
                o = get_overlap(
                    (subset_trans_chr["start"][j], subset_trans_chr["end"][j]), 
                    (intersection["start"][i], intersection["end"][i])
                )
                if o > 0:
                    report[k].append(
                        [intersection[k][i], subset_trans_chr["TPM"][j]]
                    )
    
    for k in loci.columns[3:]:
        report[k] = pd.DataFrame(report[k], columns=["posterior","TPM"]).sort_values(by="posterior").reset_index(drop=True)

    num_labels = len(loci.columns[3:])
    P_correlations = {}
    S_correlations = {}
    for k in loci.columns[3:]:
        x = np.array(report[k]["posterior"])
        y = np.array(report[k]["TPM"]) 
        P_correlations[k] = scipy.stats.pearsonr(x, y)[0]
        S_correlations[k] = scipy.stats.spearmanr(x, y)[0]

    ##########################################################################################

    plt.bar(P_correlations.keys(), P_correlations.values(), color="black", alpha=0.5)
    plt.ylabel("Pearson's Correlation")
    plt.xticks(rotation=45, fontsize=8)
    plt.tight_layout()

    plt.savefig(savedir+"/poster_trans_pearson_correl.pdf", format='pdf')
    plt.savefig(savedir+"/poster_trans_pearson_correl.svg", format='svg')
    
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')

    ##########################################################################################

    plt.bar(S_correlations.keys(), S_correlations.values(), color="black", alpha=0.5)
    plt.ylabel("Spearman's Correlation")
    plt.xticks(rotation=45, fontsize=8)
    plt.tight_layout()

    plt.savefig(savedir+"/poster_trans_spearman_correl.pdf", format='pdf')
    plt.savefig(savedir+"/poster_trans_spearman_correl.svg", format='svg')
    
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')

def posterior_transcription_enrichment(loci, trans_data, savedir, num_bins=30):
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)
    """
    get the intersection of trans_data and loci

    for k in labels:
        sort intersect based on posteriors of k
        break it down into chunks
        for each bin in chnks
            find overlap with transdata
    """
    MAPestimate = list(loci.iloc[:,3:].idxmax(axis=1))
    intersection = intersect(loci, trans_data, add_expression=True)

    labels = loci.columns[3:]

    enr_dict = {}

    strat_size = int(len(intersection) / num_bins)
    for l in labels:
        
        enr_l = []
        posterior_vector = intersection[l].astype("float").sort_values().reset_index(drop=True)
        for b in range(0, posterior_vector.shape[0], strat_size):

            subset_vector_pair = posterior_vector.iloc[b:b+strat_size].reset_index(drop=True)
            bin_range = [
                float(subset_vector_pair.iloc[subset_vector_pair.index[0]]), 
                float(subset_vector_pair.iloc[subset_vector_pair.index[-1]])
            ]

            c = loci.loc[
                (bin_range[0] <= loci[l].astype("float"))&
                (loci[l].astype("float") <= bin_range[1]), ["chr", "start", "end", l]]

            t = intersection.loc[
                (bin_range[0] <= intersection[l].astype("float"))&
                (intersection[l].astype("float") <= bin_range[1]), :]
            
            MAP_frac = len([m for m in c.index if MAPestimate[m] == l]) / len(c)

            if len(t) > 0:
                enr_l.append([bin_range[0], bin_range[1], np.mean(t["TPM"]), MAP_frac])
        
        enr_dict[l] = pd.DataFrame(enr_l, columns=["left", "right", "TPM", "map_frac"])

    list_binsdict = list(enr_dict.keys())
    num_labels = len(list_binsdict)
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(15,10))
    label_being_plotted = 0

    for i in range(n_rows):
        for j in range(n_cols):

            l = list_binsdict[label_being_plotted]
            enr = enr_dict[l]
            # print(enr)

            enr_to_plot_MAP = enr.loc[
                enr["map_frac"]>0, :] # removing bins where MAP != l
            
            enr_to_plot_notMAP = enr.loc[(enr["map_frac"]==0), :]

            axs[i,j].bar(
                [i for i in enr_to_plot_notMAP.index], enr_to_plot_notMAP.TPM, 
                align='center', width=1, edgecolor='red', color='red', alpha=0.5)

            axs[i,j].bar(
                [i for i in enr_to_plot_MAP.index], enr_to_plot_MAP.TPM, 
                align='center', width=1, edgecolor='green', color='green', alpha=0.5)

            f = UnivariateSpline(x= np.array([i for i in range(len(enr))]), y=enr.TPM, k=3, s=1000)

            axs[i,j].plot(
                [i for i in range(len(enr))], 
                f([i for i in range(len(enr))]),
                c="black", linewidth=1.5)

            axs[i,j].plot(
                [i for i in range(len(enr))], [0 for i in range(len(enr))], 
                color='black', linewidth=0.5)
                
            axs[i,j].set_xticks([])

            axs[i,j].set_title(l, fontsize=6)
            label_being_plotted +=1
    
    plt.tight_layout()
    plt.savefig(savedir+"/posterior_transcription_enrichment.pdf", format='pdf')
    plt.savefig(savedir+"/posterior_transcription_enrichment.svg", format='svg')
    sns.reset_orig
    plt.close("all")
    plt.style.use('default')    

                    
def overal_TSS_enrichment(loci, pltsavedir):
    if os.path.exists(pltsavedir) == False:
        os.mkdir(pltsavedir)
    TSS_obj = TSS_enrichment(loci, TSSdir="biovalidation/RefSeqTSS.hg38.txt", savedir=pltsavedir)
    TSS_obj.tss_enrich(m_p=False)
    # TSS_obj.tss_enrich_vs_repr()
    TSS_obj.tss_enr_vs_posterior_rank()

# def get_all_bioval(replicate_1_dir, replicate_2_dir, genecode_dir, rnaseq=None):
#     loci1, loci2 = load_data(
#         replicate_1_dir+"/parsed_posterior.csv",
#         replicate_2_dir+"/parsed_posterior.csv",
#         subset=True, logit_transform=False)

#     loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

#     gene_coords = load_gene_coords(genecode_dir)
#     if rnaseq != None:
#         trans_data = load_transcription_data(rnaseq, gene_coords)

#         trans_data = trans_data.drop(trans_data[trans_data.TPM==0].index).reset_index(drop=True)

#         trans_data_exp = trans_data[np.log10(trans_data.TPM) > 2]
#         trans_data_notexp = trans_data[np.log10(trans_data.TPM) < 0.5]

#         trans_data_exp.TPM = np.log10(trans_data_exp.TPM)
#         trans_data_notexp.TPM = np.log10(trans_data_notexp.TPM)
#         trans_data.TPM = np.log10(trans_data.TPM)

#         plot_general_transc(loci1, trans_data_exp, savedir=replicate_1_dir+"/general_transc_exp", exp=True)
#         plot_general_transc(loci1, trans_data_notexp, savedir=replicate_1_dir+"/general_transc_notexp", exp=False)

#         plot_general_transc(loci2, trans_data_exp, savedir=replicate_2_dir+"/general_transc_exp", exp=True)
#         plot_general_transc(loci2, trans_data_notexp, savedir=replicate_2_dir+"/general_transc_notexp", exp=False)

#         posterior_transcription_correlation(loci1, trans_data, savedir=replicate_1_dir+"/trans_post_correl")
#         posterior_transcription_correlation(loci2, trans_data, savedir=replicate_2_dir+"/trans_post_correl")

#     overal_TSS_enrichment(loci1, replicate_1_dir+"/tss_enr")
#     overal_TSS_enrichment(loci2, replicate_2_dir+"/tss_enr")

if __name__=="__main__":
    gene_coords = load_gene_coords("biovalidation/parsed_genecode_data_hg38_release42.csv")
    for f in os.listdir("biovalidation/RNA_seq/MCF-7/"):
        filename = "biovalidation/RNA_seq/MCF-7/" + f
        print(filename)

        if "tsv" in filename:
            data = load_transcription_data(
                filename, 
                gene_coords, csv=True)

        elif "pkl" in filename:
            data = load_transcription_data(
                filename, 
                gene_coords, csv=False)

        # print(data)


#     replicate_1_dir = "tests/cedar_runs/chmm/GM12878_R1/"
#     replicate_2_dir = "tests/cedar_runs/chmm/GM12878_R2/"

    # get_all_bioval(
    #     replicate_1_dir, replicate_2_dir, 
    #     genecode_dir="biovalidation/parsed_genecode_data_hg38_release42.csv", 
    #     rnaseq="biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv")
    
    # loci1, loci2 = load_data(
    #     replicate_1_dir+"/parsed_posterior.csv",
    #     replicate_2_dir+"/parsed_posterior.csv",
    #     subset=True, logit_transform=False)

    # loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    # tss_coords = load_TSS("biovalidation/RefSeqTSS.hg38.txt")
    # gene_coords = load_gene_coords("biovalidation/parsed_genecode_data_hg38_release42.csv")
    # trans_data = load_transcription_data("biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv", gene_coords)

    # trans_data = trans_data.drop(trans_data[trans_data.TPM==0].index).reset_index(drop=True)

    # trans_data_exp = trans_data[np.log10(trans_data.TPM) > 2]#trans_data.TPM.quantile(.75)]
    # trans_data_notexp = trans_data[np.log10(trans_data.TPM) < 0.5]#trans_data.TPM.quantile(.75)]

    # trans_data_exp.TPM = np.log10(trans_data_exp.TPM)
    # trans_data_notexp.TPM = np.log10(trans_data_notexp.TPM)
    # trans_data.TPM = np.log10(trans_data.TPM)

    # plt.hist(np.log10(trans_data_exp.TPM), bins=100)
    # plt.hist(np.log10(trans_data_notexp.TPM), bins=100)

    # plot_general_transc(loci1, gene_coords, "--")
    # plot_general_transc(loci1, trans_data_exp, "--")
    # plot_general_transc(loci1, trans_data_notexp, "--")


    # posterior_transcribed_enrichment(loci1, gene_coords)
    # posterior_transcribed_enrichment(loci1, trans_data_exp)
    # posterior_transcribed_enrichment(loci1, trans_data_notexp)

    # overal_TSS_enrichment(loci1, replicate_1_dir+"/tss_enr")

    # posterior_transcription_correlation(loci1, trans_data_exp)
    # posterior_transcription_correlation(loci1, trans_data_notexp)
    # posterior_transcription_correlation(loci1, trans_data)