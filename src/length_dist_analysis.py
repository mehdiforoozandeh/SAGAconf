import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

def read_bed(bedfile):
    bed =  pd.read_csv(bedfile, sep="\t", header=None)
    bed.columns = ["chr", "start", "end", "label", "?", "??", "???", "????", "?????"]
    return bed

def get_length(bed):
    bed["length"] = bed.end - bed.start
    bed = bed.drop(["start","end","?", "??", "???", "????", "?????"], axis=1)
    return bed

def clean_bed(bed):
    bed = bed.loc[bed["chr"] != "chrUn"]
    bed = bed.loc[bed["chr"] != "chrM"]
    bed = bed.loc[bed["chr"] != "chrEBV"]
    bed = bed.loc[bed['chr'].str.contains("random")==False]
    bed = bed.reset_index(drop=True)
    return bed 

def ECDF(bed):
    print(bed)
    q = bed["length"].quantile(0.999)
    bed = bed.loc[bed["length"] < q]

    sns.ecdfplot(data=bed, x="length", hue="label")
    plt.tight_layout()
    plt.show()

def chr_lengthdist(bed):
    unique_chr = list(np.unique(np.array(bed["chr"])))
    unique_labels = list(np.unique(np.array(bed["label"])))

    
    unique_chr = [c for c in unique_chr if "chrUn" not in c and "random" not in c and "chrM" not in c and "chrEBV" not in c]
    print(unique_chr)
    
    num_labels = len(unique_labels)
    n_cols = math.floor(math.sqrt(num_labels))
    n_rows = math.ceil(num_labels / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
    label_being_plotted = 0

    length_dist = {}
    for i in range(n_rows):
        for j in range(n_cols):
            l = unique_labels[label_being_plotted]
            if l not in length_dist.keys():
                length_dist[l] = {}
            for c in unique_chr:
                selection = bed.loc[(bed["chr"]==c) &( bed["label"]==l)]
                length_dist[l][c] = np.array(selection["length"])

            axs[i,j].boxplot([v for v in length_dist[l].values()], patch_artist = True, notch ='True', showfliers=False)
            axs[i,j].set_title(l, fontsize=8)
            axs[i,j].set_xticks(ticks=[i+1 for i in range(len(length_dist[l].keys()))])
            axs[i,j].set_xticklabels([k for k in length_dist[l].keys()], rotation=90, fontsize=7)
            label_being_plotted+=1

    plt.tight_layout()
    plt.show()

if __name__=="__main__":
    bed = read_bed("tests/length_dist_anal/ENCFF194OGV.bed")
    len_bed = get_length(bed)
    len_bed = clean_bed(len_bed)
    print(len_bed.loc[len_bed["length"]==np.max(len_bed.length)])

    # chr_lengthdist(len_bed)
    ECDF(len_bed)

    bed = read_bed("tests/length_dist_anal/ENCFF505SNZ.bed")
    len_bed = get_length(bed)
    len_bed = clean_bed(len_bed)
    print(len_bed.loc[len_bed["length"]==np.max(len_bed.length)])

    # chr_lengthdist(len_bed)
    ECDF(bed)
