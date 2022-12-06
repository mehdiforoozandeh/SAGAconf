from _reproducibility import *
from reports import *
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

def load_transcription_data(file, gene_coord):

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

def intersect(loci, gene_coord):
    cols = loci.columns
    
    loci = loci.values.tolist()

    intersect = []
    for i in range(len(loci)):

        gene_coord_chr = gene_coord.loc[gene_coord["chr"] == loci[i][0], :]
        gene_coord_chr = gene_coord_chr.values.tolist()
        o = False
        for j in range(len(gene_coord_chr)):
            statement = bool(
                    bool(gene_coord_chr[j][1] <= loci[i][1] <= gene_coord_chr[j][2]) or
                    bool(gene_coord_chr[j][1] <= loci[i][2] <= gene_coord_chr[j][2]))

            if statement:
                o = True
                continue
        
        if o:
            intersect.append(loci[i])
    
    intersect = pd.DataFrame(intersect, columns=cols)
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

def plot_general_transc(loci, gene_coords, savedir):
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
    plt.show()

def posterior_transcribed_enrichment(loci, gene_coords, strat_size="def"):
    """
    define bins of size B for posterior

    for posteriors in range B
        for k in labels

            what percentage of positions in that range overlap with protein-coding genes
                t = number of intersected regions with posterior in range B
                e = number of all regions with posterior in range B
    """

    print(loci)
    print(gene_coords)
    print("intersecting")
    intersection = intersect(loci, gene_coords)

    if strat_size=="def":
        strat_size = int(len(loci)/40)

    reports = {}
    for k in loci.columns[3:]:
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

            if e != 0 and t != 0:
                obs = float(t) / len(intersection)
                exp = float(e) / len(posterior_vector)
                reports[k].append([bin_range, np.log(float(obs/exp))])
    
    num_labels = len(loci.columns[3:])
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
    colors = [i for i in get_cmap('tab20').colors]
    label_being_plotted = 0
    
    for i in range(n_rows):
        for j in range(n_cols):
            k = loci.columns[3:][label_being_plotted]

            x = [(kk[0][0]+kk[0][1])/2 for kk in reports[k]]
            y = [kk[1] for kk in reports[k]]

            axs[i,j].bar(
                x, y,  
                width=0.05,
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

    plt.show()

def posterior_transcription_correlation(loci, trans_data, k, stratsize = "def"):
    posterior_k = loci.loc[:, ["chr", "start", "end", k]]
    posterior_k = condense_posterior(posterior_k)

    report = []

    for b in range(len(posterior_k)):
        trans_data_subset = trans_data.loc[(trans_data["chr"] == posterior_k["chr"][b]),:].reset_index(drop=True)

        for t in range(len(trans_data_subset)):
            
            o = get_overlap(
                (posterior_k["start"][b], posterior_k["end"][b]), 
                (trans_data_subset["start"][t], trans_data_subset["end"][t])
            )
    
            if o>0:
                report.append(
                    [posterior_k[k][b], (float(o)/(posterior_k["end"][b] - posterior_k["start"][b]))*trans_data_subset["TPM"][t]]
                )

    report = pd.DataFrame(report, columns=["posterior", "transcription"]).sort_values(by="posterior").reset_index(drop=True)
    
    return report

def plot_poster_transc_correl(loci, trans_data, savedir, strat_size="def"):
    reports = {}
    for k in loci.columns[3:]:
        print(k)
        reports[k] = posterior_transcription_correlation(loci, trans_data, k)

    for k, report in reports.items():
        print(k, len(report))
        if strat_size == "def":
            strat_size = int(len(report)/40)

        binned_data = []
        for i in range(0, len(report), strat_size):
            report_subset = report.iloc[i:i+strat_size, :].reset_index(drop=True)
            bin_range = (report_subset["posterior"][0], list(report_subset["posterior"])[-1])
            mean_p = np.mean(report_subset["posterior"])
            mean_transc = np.mean(report_subset["transcription"])
            std_transc = np.std(report_subset["transcription"])

            binned_data.append([
                bin_range[0], bin_range[1], mean_p, mean_transc, std_transc
                ])
        
        # binned_data = pd.DataFrame(binned_data, columns=["start", "end", "mean_p", "mean_transc", "std_transc"])
        print(binned_data)

        #######################################################################################
        colors = [i for i in get_cmap('tab20').colors]
        list_binsdict = list(binned_data.keys())
        num_labels = len(list_binsdict)
        n_cols = math.floor(math.sqrt(num_labels))
        n_rows = math.ceil(num_labels / n_cols)

        fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
        label_being_plotted = 0
        ci = 0
        for i in range(n_rows):
            for j in range(n_cols):
                k = list_binsdict[label_being_plotted]
                x = [(kk[0]+kk[1])/2 for kk in binned_data[k]]
                y = [kk[3] for kk in binned_data[k]]
                
                suggested_width = float(np.array(x).max() - np.array(x).min()) / len(x)

                axs[i,j].bar(
                    x, y,  
                    width=0.05,
                    color=colors[ci],
                    alpha=0.5)
                
                try:
                    if len(x)<=3:
                        f = UnivariateSpline(x= np.array(x), y=y, k=1)
                    else:
                        f = UnivariateSpline(x= np.array(x), y=y, k=3)

                    axs[i,j].plot(
                        x, 
                        f(x),
                        label=k, c=colors[ci], 
                        linewidth=1)
                except:
                    pass

                axs[i,j].set_title(k, fontsize=7)
                ci += 1
                label_being_plotted += 1
                
        
        # plt.savefig('{}/rep_vs_TSS_enrichment_bars.pdf'.format(savedir), format='pdf')
        # plt.savefig('{}/rep_vs_TSS_enrichment_bars.svg'.format(savedir), format='svg')
        # plt.clf()
        plt.tight_layout()
        plt.show()

def overal_TSS_enrichment(loci, pltsavedir):
    if os.path.exists(pltsavedir) == False:
        os.mkdir(pltsavedir)
    TSS_obj = TSS_enrichment(loci, TSSdir="biovalidation/RefSeqTSS.hg38.txt", savedir=pltsavedir)
    TSS_obj.tss_enrich(m_p=False)
    TSS_obj.tss_enrich_vs_repr()


if __name__=="__main__":
    replicate_1_dir = "tests/cedar_runs/chmm/GM12878_R1/"
    replicate_2_dir = "tests/cedar_runs/chmm/GM12878_R2/"
    
    loci1, loci2 = load_data(
        replicate_1_dir+"/parsed_posterior.csv",
        replicate_2_dir+"/parsed_posterior.csv",
        subset=True)

    loci1, loci2 = process_data(loci1, loci2, replicate_1_dir, replicate_2_dir, mnemons=True, match=False)

    tss_coords = load_TSS("biovalidation/RefSeqTSS.hg38.txt")
    gene_coords = load_gene_coords("biovalidation/parsed_genecode_data_hg38_release42.csv")
    trans_data = load_transcription_data("biovalidation/RNA_seq/GM12878/preferred_default_ENCFF240WBI.tsv", gene_coords)

    # trans_data = trans_data.drop(trans_data[trans_data.TPM==0].index).reset_index(drop=True)

    # print(np.sum(trans_data.TPM))
    # plt.hist(np.log10(trans_data.TPM+1), bins=100)
    # plt.show()

    trans_data_exp = trans_data[np.log10(trans_data.TPM) > 2]#trans_data.TPM.quantile(.75)]
    trans_data_notexp = trans_data[np.log10(trans_data.TPM) < 0.5]#trans_data.TPM.quantile(.75)]

    # plt.hist(np.log10(trans_data_exp.TPM), bins=100)
    # plt.hist(np.log10(trans_data_notexp.TPM), bins=100)

    plot_general_transc(loci1, gene_coords, "--")
    plot_general_transc(loci1, trans_data_exp, "--")
    plot_general_transc(loci1, trans_data_notexp, "--")

    posterior_transcribed_enrichment(loci1, gene_coords)
    posterior_transcribed_enrichment(loci1, trans_data_exp)
    posterior_transcribed_enrichment(loci1, trans_data_notexp)

    # overal_TSS_enrichment(loci1, replicate_1_dir+"/tss_enr")

    # plot_poster_transc_correl(loci1, trans_data, "--")
    # plot_general_transc(loci2, gene_coords, "--")
    # plot_poster_transc_correl(loci2, trans_data, "--")