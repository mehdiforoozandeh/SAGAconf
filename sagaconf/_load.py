"""Loading and intersecting the parsed posteriors (kept free of the plotting and model imports,
so the CLI can load the data while those are imported)."""
import os

import numpy as np

from ._utils import logit_array
from .run import intersect_parsed_posteriors, read_mnemonics


def load_data(posterior1_dir, posterior2_dir, subset=False, logit_transform=False, force_WG=False):
    print("loading and intersecting")
    loci_1, loci_2 = intersect_parsed_posteriors(
        posterior1_dir, 
        posterior2_dir)

    if subset and force_WG==False:
        loci_1 = loci_1.loc[loci_1["chr"]=="chr21"].reset_index(drop=True)
        loci_2 = loci_2.loc[loci_2["chr"]=="chr21"].reset_index(drop=True)

    print("the shapes of the input matrices are: {}, {}".format(str(loci_1.shape), str(loci_2.shape)))

    if logit_transform:
        loci_1.iloc[:,3:] = logit_array(np.array(loci_1.iloc[:,3:]))
        loci_2.iloc[:,3:] = logit_array(np.array(loci_2.iloc[:,3:]))

    return loci_1, loci_2


def process_data(loci_1, loci_2, replicate_1_dir, replicate_2_dir, mnemons=True, vm="NA", bm="NA", match=False, custom_order=True):
    if match:
        from ._cluster_matching import Hungarian_algorithm, connect_bipartite, overlap_matrix
    # print('generating confmat 1 ...')
    num_labels = loci_1.shape[1]-3

    loci_1.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    loci_2.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    
    if mnemons:
        print("loading mnemonics...")
        if os.path.exists(
            "/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics_rep1.txt") and os.path.exists(
            "/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics_rep2.txt"):

            print("reading concat mnemonics")
            loci_1_mnemon = read_mnemonics("/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics_rep1.txt")
            loci_2_mnemon = read_mnemonics("/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics_rep2.txt")
        else:
            print("reading mnemonics")
            if bm == "NA":
                bm = "/".join(replicate_1_dir.split("/")[:-1])+"/mnemonics.txt"
            if vm == "NA":
                vm = "/".join(replicate_2_dir.split("/")[:-1])+"/mnemonics.txt"

            loci_1_mnemon = read_mnemonics(bm)
            loci_2_mnemon = read_mnemonics(vm)

        mnemon1_dict = {}
        for i in loci_1_mnemon:
            if len(i.split("_")) == 2:
                mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]
            if len(i.split("_")) == 3:
                mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4] + '_' + i.split("_")[2][:3]

        mnemon2_dict = {}
        for i in loci_2_mnemon:
            if len(i.split("_")) == 2:
                mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4]
            if len(i.split("_")) == 3:
                mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:4] + '_' + i.split("_")[2][:3]

        #handle missing mnemonics
        for i in range(num_labels):
            if str(i) not in mnemon1_dict.keys():
                mnemon1_dict[str(i)] = str(f"{i}_Unkn")
            if str(i) not in mnemon2_dict.keys():
                mnemon2_dict[str(i)] = str(f"{i}_Unkn")
        
        if match:
            conf_mat = overlap_matrix(loci_1, loci_2, type="IoU")

            assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')
            for i in range(len(assignment_pairs)):
                assignment_pairs[i] = (mnemon1_dict[str(assignment_pairs[i][0])], mnemon2_dict[str(assignment_pairs[i][1])])
            print(assignment_pairs)



            loci_1, loci_2 = \
                connect_bipartite(loci_1, loci_2, assignment_pairs, mnemon=True)
            
            print('connected barpartite')

        else:
            loci_1.columns = list(loci_1.columns[:3]) + [mnemon1_dict[c.replace("posterior","")] for c in loci_1.columns[3:]]
            loci_2.columns = list(loci_2.columns[:3]) + [mnemon2_dict[c.replace("posterior","")] for c in loci_2.columns[3:]]

    else:
        if match:
            conf_mat = overlap_matrix(loci_1, loci_2, type="IoU")

            assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')
            loci_1, loci_2 = \
                connect_bipartite(loci_1, loci_2, assignment_pairs, mnemon=False)

            print('connected barpartite')
    
    if mnemons and custom_order:
        SORT_ORDER = {"Prom": 0, "Prom_fla":1, "Enha":2, "Enha_low":3, "Biva":4, "Tran":5, "Cons":6, "Facu":7, "K9K3":8, "Quie":9, "Unkn":10}
        try:
            new_columns = []
            for c in loci_1.columns[3:]:
                l = "_".join(c.split("_")[1:])
                new_columns.append(str(SORT_ORDER[l])+"_"+c)
                
            new_columns.sort()
            for i in range(len(new_columns)):
                new_columns[i] = new_columns[i][2:]

            loci_1 = loci_1[["chr", "start", "end"] + new_columns]
        except:
            pass

        ##########################################################################################
        ##########################################################################################
        try:
            new_columns = []
            for c in loci_2.columns[3:]:
                l = "_".join(c.split("_")[1:])
                new_columns.append(str(SORT_ORDER[l])+"_"+c)
                
            new_columns.sort()
            for i in range(len(new_columns)):
                new_columns[i] = new_columns[i][2:]

            loci_2 = loci_2[["chr", "start", "end"] + new_columns]
        except:
            pass

    return loci_1, loci_2
