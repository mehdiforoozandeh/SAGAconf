import numpy as np
import pandas as pd
from ._cluster_matching import IoU_overlap
from ._reproducibility import NMI_from_matrix


def get_coverage(loci, query_label="all"):
    MAP = loci.iloc[:,3:].idxmax(axis=1)
    coverage = dict(zip(list(loci.columns[3:]), [0 for _ in range(len(loci.columns[3:]))]))

    if query_label == "all":
        for c in list(loci.columns[3:]):
            coverage[c] = len(MAP.loc[MAP == c]) / len(MAP)
        return coverage

    else:
        return len(MAP.loc[MAP == query_label]) / len(MAP)


def get_agr_nonsymmetric(MAP1, MAP2, query_label1, query_label2):
    record = [0,0] # [agreement_count, disagreement_count]

    if isinstance(MAP1.index, pd.RangeIndex) and MAP1.index.start == 0 and MAP1.index.step == 1 \
            and len(MAP2) >= len(MAP1) and isinstance(MAP2.index, pd.RangeIndex) \
            and MAP2.index.start == 0 and MAP2.index.step == 1:
        sel = MAP1.to_numpy() == query_label1
        hit = int(np.count_nonzero(MAP2.to_numpy()[:len(MAP1)][sel] == query_label2))
        record[0] = hit
        record[1] = int(np.count_nonzero(sel)) - hit
    else:
        for i in range(len(MAP1)):
            if MAP1[i] == query_label1:

                if MAP2[i] == query_label2:
                    record[0] += 1
                else:
                    record[1] += 1

    try:
        agr = float(record[0] / (record[0] + record[1]))
    except:
        agr = 0
    return agr


def granularity_vs_agreement_nonsymmetric(loci_1, loci_2, k, disregard_posterior=True):
    '''
    input locis must NOT be matched -- optionally with updated mnemonics
    '''
    num_labels = loci_1.shape[1]-3
    MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
    MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
    if disregard_posterior and _onehot_merging_is_exact(loci_1, loci_2, MAP1, MAP2):
        return _granularity_from_counts(loci_1, loci_2, k, MAP1, MAP2)

    if disregard_posterior:
        # to avoid MAP changes during the merging process, disregard probability granularity and convert 
        # to hard zero vs one matrix
        
        for c in loci_1.columns[3:]:
            is_c = MAP1==c
            loci_1.loc[is_c, c] = 1
            loci_1.loc[~is_c, c] = 0
        
        for c in loci_2.columns[3:]:
            is_c = MAP2==c
            loci_2.loc[is_c, c] = 1
            loci_2.loc[~is_c, c] = 0

    confmat = IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False)
    
    sorted_k_vector = confmat.loc[k,:].sort_values(ascending=False)
    # print(k, sorted_k_vector)

    query_label = str(sorted_k_vector.index[0])
    coverage_record = [0]
    agreement_record = [0]

    for j in range(0, len(sorted_k_vector)):
        if j == 0:
            pass
        else:
            loci_2[query_label + "+" + str(sorted_k_vector.index[j])] = \
                loci_2[query_label] + loci_2[str(sorted_k_vector.index[j])]

            loci_2 = loci_2.drop([query_label, str(sorted_k_vector.index[j])], axis=1)

            query_label = query_label + "+" + str(sorted_k_vector.index[j])
    
        if j == 0:  # loci_1 does not change inside this loop
            MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
        MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)

        jth_order_agr = get_agr_nonsymmetric(MAP1, MAP2, query_label1=k, query_label2=query_label)
        # same as get_coverage(loci_2, query_label=query_label), reusing MAP2
        jth_order_cov = len(MAP2.loc[MAP2 == query_label]) / len(MAP2)

        coverage_record.append(jth_order_cov)
        agreement_record.append(jth_order_agr)
    return coverage_record, agreement_record, sorted_k_vector.index


def _onehot_merging_is_exact(loci_1, loci_2, MAP1, MAP2):
    """True if the hard 0/1 matrices keep each row's MAP and merged-column names cannot clash,
    so every MAP2 after merging is the merged label for rows whose MAP2 was merged, else MAP2."""
    labels1, labels2 = loci_1.columns[3:], loci_2.columns[3:]
    return (labels1.is_unique and labels2.is_unique and len(labels2) > 0
            and all(isinstance(c, str) for c in labels1) and all(isinstance(c, str) and "+" not in c for c in labels2)
            and len(MAP1) == len(MAP2) and len(MAP1) > 0
            and not MAP1.isna().any() and not MAP2.isna().any()
            and isinstance(MAP1.index, pd.RangeIndex) and MAP1.index.start == 0 and MAP1.index.step == 1
            and isinstance(MAP2.index, pd.RangeIndex) and MAP2.index.start == 0 and MAP2.index.step == 1)


def _granularity_from_counts(loci_1, loci_2, k, MAP1, MAP2):
    """granularity_vs_agreement_nonsymmetric with disregard_posterior, from integer counts.

    The IoU of the 0/1 matrices depends only on the MAPs, which binarising keeps. Merging the
    first j+1 labels of sorted_k_vector makes MAP2 the merged label exactly where the original
    MAP2 is one of them, so the agreement and coverage counts are cumulative label counts."""
    confmat = IoU_overlap(loci_1, loci_2, w=0, symmetric=True, soft=False, maps=(MAP1, MAP2))
    sorted_k_vector = confmat.loc[k,:].sort_values(ascending=False)

    rank2 = sorted_k_vector.index.get_indexer(pd.Index(MAP2.to_numpy(), dtype=object))
    is_k = MAP1.to_numpy() == k
    n_k = int(np.count_nonzero(is_k))
    cov_counts = np.cumsum(np.bincount(rank2, minlength=len(sorted_k_vector)))
    agr_counts = np.cumsum(np.bincount(rank2[is_k], minlength=len(sorted_k_vector)))

    coverage_record = [0]
    agreement_record = [0]
    for j in range(0, len(sorted_k_vector)):
        hit = int(agr_counts[j])
        try:
            agr = float(hit / (hit + (n_k - hit)))
        except:
            agr = 0
        coverage_record.append(int(cov_counts[j]) / len(MAP2))
        agreement_record.append(agr)
    return coverage_record, agreement_record, sorted_k_vector.index


def merging_gain(joint_overlap):
    coverage1 = {k:sum(joint_overlap.loc[k,:]) for k in joint_overlap.index}
    coverage2 = {k:sum(joint_overlap.loc[:,k]) for k in joint_overlap.columns}

    gain = pd.DataFrame(
        np.zeros([joint_overlap.shape[0], joint_overlap.shape[0]]), 
        index=joint_overlap.index,
        columns=joint_overlap.index)
    
    for i in joint_overlap.index:
        for j in joint_overlap.index:
            joint_overlap_prime = joint_overlap.copy()

            if i == j:
                gain.loc[i, j] = 0

            else:
                MI0 = NMI_from_matrix(joint_overlap_prime, return_MI=True)
                H0 = -1 * np.sum([
                    np.sum(joint_overlap_prime.loc[p, :]) * np.log2(np.sum(joint_overlap_prime.loc[p, :]) + 1e-9) for p in joint_overlap_prime.index
                ])

                joint_prime_ij = joint_overlap.loc[i, :] + joint_overlap.loc[j, :]
                coverage_ij = np.sum(joint_prime_ij)

                joint_overlap_prime.loc[str(i) + " + " + str(j)] = joint_prime_ij
                joint_overlap_prime = joint_overlap_prime.drop([i, j], axis=0)

                MI1 = NMI_from_matrix(joint_overlap_prime, return_MI=True)
                H1 = -1 * np.sum([
                    np.sum(joint_overlap_prime.loc[p, :]) * np.log2(np.sum(joint_overlap_prime.loc[p, :]) + 1e-9) for p in joint_overlap_prime.index
                ])
                
                gain.loc[i, j] = 1 - ((MI1 - MI0) / (H1 - H0))

    return gain


def merge_clusters(joint, loci_1, loci_2, r1=True):
    m = 0

    #########################################################################################################
    if r1:
        num_labels = loci_1.shape[1]-3
        r1_gain = pd.DataFrame(merging_gain(joint), columns=joint.index, index=joint.index)

        max_value = r1_gain.max().max()

        # Find the indices of the maximum value in the similarity matrix
        max_index = np.where(r1_gain == max_value)

        # Get the row and column labels of the maximum value
        row_label = r1_gain.index[max_index[0][0]]
        col_label = r1_gain.columns[max_index[1][0]]

        print(f"The most similar pair is ({row_label}, {col_label}) with a similarity of {max_value}")

        linkage_1 = np.array([[max_index[0][0], max_index[1][0]]])

        merged_label_ID_1 = {}
        labels = loci_1.iloc[:, 3:].columns
        for i in range(num_labels):
            merged_label_ID_1[i] = labels[i]
        
        to_be_merged_1 = [
            merged_label_ID_1[int(linkage_1[m, 0])],
            merged_label_ID_1[int(linkage_1[m, 1])],
        ]

        merged_label_ID_1[num_labels + m] = str(
            merged_label_ID_1[int(linkage_1[m, 0])] + "+" + merged_label_ID_1[int(linkage_1[m, 1])]
        )

        loci_1[merged_label_ID_1[num_labels + m]] = \
            loci_1[to_be_merged_1[0]] + loci_1[to_be_merged_1[1]]
        loci_1 = loci_1.drop(to_be_merged_1, axis=1)

    #########################################################################################################
    else:
        num_labels = loci_2.shape[1]-3
        r2_gain = pd.DataFrame(merging_gain(joint.T), columns=joint.columns, index=joint.columns)

        max_value = r2_gain.max().max()

        # Find the indices of the maximum value in the similarity matrix
        max_index = np.where(r2_gain == max_value)

        # Get the row and column labels of the maximum value
        row_label = r2_gain.index[max_index[0][0]]
        col_label = r2_gain.columns[max_index[1][0]]

        print(f"The most similar pair is ({row_label}, {col_label}) with a similarity of {max_value}")

        linkage_2 = np.array([[max_index[0][0], max_index[1][0]]])

        merged_label_ID_2 = {}
        labels = loci_2.iloc[:, 3:].columns
        for i in range(num_labels):
            merged_label_ID_2[i] = labels[i]

        to_be_merged_2 = [
            merged_label_ID_2[int(linkage_2[m, 0])],
            merged_label_ID_2[int(linkage_2[m, 1])],
        ]

        merged_label_ID_2[num_labels + m] = str(
            merged_label_ID_2[int(linkage_2[m, 0])] + "+" + merged_label_ID_2[int(linkage_2[m, 1])]
        )

        loci_2[merged_label_ID_2[num_labels + m]] = \
            loci_2[to_be_merged_2[0]] + loci_2[to_be_merged_2[1]]
        loci_2 = loci_2.drop(to_be_merged_2, axis=1)

    return loci_1, loci_2


if __name__=="__main__":
    run_on_subset = True
    mnemons = True
    symmetric = False

    replicate_1_dir = "tests/cedar_runs/chmm/MCF7_R1/"
    replicate_2_dir = "tests/cedar_runs/chmm/MCF7_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons, symmetric)
    run(replicate_2_dir, replicate_1_dir, run_on_subset, mnemons, symmetric)

    replicate_1_dir = "tests/cedar_runs/segway/MCF7_R1/"
    replicate_2_dir = "tests/cedar_runs/segway/MCF7_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons, symmetric)
    run(replicate_2_dir, replicate_1_dir, run_on_subset, mnemons, symmetric)

    replicate_1_dir = "tests/cedar_runs/chmm/GM12878_R1/"
    replicate_2_dir = "tests/cedar_runs/chmm/GM12878_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons, symmetric)
    run(replicate_2_dir, replicate_1_dir, run_on_subset, mnemons, symmetric)

    replicate_1_dir = "tests/cedar_runs/segway/GM12878_R1/"
    replicate_2_dir = "tests/cedar_runs/segway/GM12878_R2/"
    run(replicate_1_dir, replicate_2_dir, run_on_subset, mnemons, symmetric)
    run(replicate_2_dir, replicate_1_dir, run_on_subset, mnemons, symmetric)
