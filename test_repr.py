from cProfile import run
from os import lseek
from timeit import timeit
from _reproducibility import *
import pandas as pd
from matchlabels import *
from _visualize import visualize
from run import *

def prepare_data(replicate_1_dir, replicate_2_dir, cluster=False, num_labels=None, plot_dend=False):
    if cluster:
        loci_1, loci_2 = cluster_matching(\
                    replicate_1_dir, replicate_2_dir, 
                    n_clust=num_labels, matching_strategy='conf_matrix', plot_dendogram=plot_dend)
    else:
        loci_1, loci_2 = plain_seg_matching(\
                    replicate_1_dir, replicate_2_dir, 
                    matching_strategy='conf_matrix')
    
    return loci_1, loci_2

def prepare_concatenated_results(rep_dir_1, rep_dir_2):
    parsed_posterior_1 = pd.read_csv(
        '{}/parsed_posterior.csv'.format(rep_dir_1)).drop('Unnamed: 0', axis=1)
    parsed_posterior_2 = pd.read_csv(
        '{}/parsed_posterior.csv'.format(rep_dir_2)).drop('Unnamed: 0', axis=1)

    loci_1 = parsed_posterior_1.iloc[:int(len(parsed_posterior_1)/2),:]
    loci_2 = parsed_posterior_2.iloc[int(len(parsed_posterior_2)/2):,:]

    loci_2 = loci_2.reset_index(drop=True)

    return loci_1, loci_2

def report_reproducibility(rep_1_dir, rep_2_dir, pltsavedir, concat=False, verbose=True):
    if concat:
        loci_1, loci_2 = prepare_concatenated_results(rep_1_dir, rep_2_dir)
    else:
        loci_1, loci_2 = prepare_data(
            rep_1_dir, rep_2_dir, cluster=False, num_labels=9, plot_dend=False)
    
    
    if os.path.exists(pltsavedir)==False:
        os.mkdir(pltsavedir)

    if os.path.exists(pltsavedir+"/agreement")==False:
        os.mkdir(pltsavedir+"/agreement")

    if os.path.exists(pltsavedir+"/cc")==False:
        os.mkdir(pltsavedir+"/cc")

    if os.path.exists(pltsavedir+"/calib")==False:
        os.mkdir(pltsavedir+"/calib")

    if os.path.exists(pltsavedir+"/sankey")==False:
        os.mkdir(pltsavedir+"/sankey")

    agr = Agreement(loci_1, loci_2, pltsavedir+"/agreement")
    cc = correspondence_curve(loci_1, loci_2, pltsavedir+"/cc")
    vis = sankey(loci_1, loci_2, pltsavedir+"/sankey")
    calb = posterior_calibration(loci_1, loci_2, log_transform=False, ignore_overconf=True, filter_nan=True)

    if verbose:
        # print('perlabel agreement:    ', agr.per_label_agreement())
        print('general agreement:    ', agr.general_agreement())
        # print('expected agreement:    ', agr.expected_agreement)
        # print('perlabel o/e ratio:    ', agr.per_label_OE_ratio())
        print('general o/e ratio:    ', agr.general_OE_ratio(log_transform=False))
        # print('perlabel CK:    ', agr.per_label_cohens_kappa())
        print('general CK:    ', agr.general_cohens_kappa())

    agr.plot_agreement()
    agr.plot_CK()
    agr.plot_OE()

    cc.plot_curve(plot_general=True, merge_plots=False)
    vis.sankey_diag()
    
def test_calibration(replicate_1_dir, replicate_2_dir):
    loci_1, loci_2 = intersect_parsed_posteriors(
        replicate_1_dir+"/parsed_posterior.csv", 
        replicate_2_dir+"/parsed_posterior.csv")

    num_labels = loci_1.shape[1]-3
    loci_1.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    loci_2.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]

    print('generating confmat 1')
    
    conf_mat = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=False)

    assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')
    print(assignment_pairs)


    loci_1, loci_2 = \
        connect_bipartite(loci_1, loci_2, assignment_pairs)
    # label_coverage_1 = {}

    # for c in loci_1.columns[3:]: 
    #     label_coverage_1[c] = len(
    #         loci_1.loc[
    #             (loci_1.iloc[:,3:].idxmax(axis=1) == c), :
    #         ]
    #     ) / len(loci_1) 
    # print(label_coverage_1)

    # label_coverage_2 = {}

    # for c in loci_2.columns[3:]: 
    #     label_coverage_2[c] = len(
    #         loci_2.loc[
    #             (loci_2.iloc[:,3:].idxmax(axis=1) == c), :
    #         ]
    #     ) / len(loci_2) 
    # print(label_coverage_2)
    print('connected barpartite')
    
    agr = Agreement(loci_1, loci_2, "agreement")
    print('general agreement:    ', agr.per_label_agreement())
    # print(loci_1)
    # print(loci_2)
    
    calb = posterior_calibration(
        loci_1, loci_2, log_transform=False, ignore_overconf=False, filter_nan=True, oe_transform=True, savedir="tests/reprod_plots")
    calibrated_loci_1 = calb.perlabel_calibration_function(degree=5, num_bins=10, return_caliberated_matrix=True)

    # vext = validate_EXT()
    # vext.read_feat_agg_enrichment(
    #     replicate_1_dir+"/aggregations/feature_aggregation.tab", 
    #     agr.per_label_agreement())
    # vext.TSS_vs_agreement_plot()
    # tss_coords = vext.TSS_from_gtf("gencode.v29.primary_assembly.annotation_UCSC_names.gtf")

    # print(tss_coords)
    # # vext.TSS_vs_calibrated_plot(calibrated_loci_1)
    # vext.TSS_enrich_vs_reproducibility(calibrated_loci_1)

def subsampleloci(replicate_1_dir, samplesize=100):
    with open(replicate_1_dir+"/parsed_posterior_short.csv", 'w') as newfile:
        with open(replicate_1_dir+"/parsed_posterior.csv", 'r') as origfile:
            lines = origfile.readlines()
            for l in range(0, len(lines), samplesize):
                newfile.write(lines[l])
    

def TSS_test(loci_1, loci_2, tss, replicate_1_dir, replicate_2_dir):
    num_labels = loci_1.shape[1]-3
    loci_1.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    loci_2.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]

    print('generating confmat 1')
    
    conf_mat = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=False)

    assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')

    loci_1_mnemon = read_mnemonics(replicate_1_dir.replace("/parsed_posterior.csv","/mnemonics.txt"))
    loci_2_mnemon = read_mnemonics(replicate_2_dir.replace("/parsed_posterior.csv","/mnemonics.txt"))
    
    mnemon1_dict = {}
    for i in loci_1_mnemon:
        mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:5]

    mnemon2_dict = {}
    for i in loci_2_mnemon:
        mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:5]
    
    for i in range(len(assignment_pairs)):
        assignment_pairs[i] = (mnemon1_dict[str(assignment_pairs[i][0])], mnemon2_dict[str(assignment_pairs[i][1])])
    print(assignment_pairs)

    loci_1, loci_2 = \
        connect_bipartite(loci_1, loci_2, assignment_pairs, mnemon=True)

    print('connected barpartite')
    print(loci_1.shape[0])
    # plot_posterior_histogram(loci_1)
    # plot_posterior_track(loci_1)

    # plot_posterior_histogram(loci_2)
    # plot_posterior_track(loci_2)


    vis = sankey(loci_1, loci_2, "tests/clb_1_oe")
    vis.heatmap() #make this without mnemonics
    
    # calb = posterior_calibration(loci_1, loci_2, "tests/clb_1_logoe", log_transform=True, ignore_overconf=False, filter_nan=True)
    # loci_1_c = calb.perlabel_calibration_function(method="isoton_reg", degree=3, num_bins=10, return_caliberated_matrix=True)
    # plot_posterior_distribution(loci_1_c)

    # print(loci_1_c)
    # calb = posterior_calibration(loci_2, loci_1, "tests/clb_2_logoe", log_transform=True, ignore_overconf=True, filter_nan=True)
    # loci_2_c = calb.perlabel_calibration_function(method="isoton_reg", degree=3, num_bins=10, return_caliberated_matrix=True)
    # print(loci_2_c)

    calb = posterior_calibration(loci_1, loci_2, "tests/clb_1_oe", log_transform=False, ignore_overconf=False, filter_nan=True)
    loci_1_c = calb.perlabel_calibration_function(method="isoton_reg", degree=3, num_bins=10, return_caliberated_matrix=True)
    print(loci_1_c)
    
    # calb = posterior_calibration(loci_2, loci_1, "tests/clb_2_oe", log_transform=False, ignore_overconf=False, filter_nan=True)
    # loci_2_c = calb.perlabel_calibration_function(method="isoton_reg", degree=3, num_bins=10, return_caliberated_matrix=True)
    # print(loci_2_c)

    # tssenr1 = tss_enrich(loci_1_c, tss)
    # plt.bar(list(tssenr1.keys()), list(tssenr1.values()))
    # plt.ylabel("log(O/E) TSS enrichment")
    # plt.xticks(rotation=45, fontsize=7)
    # plt.show()

    # tss_enrich_vs_repr(loci_1_c, tss, m_p=True)

    # tssenr2 = tss_enrich(loci_2_c, tss)
    # plt.bar(list(tssenr2.keys()), list(tssenr2.values()))
    # plt.ylabel("log(O/E) TSS enrichment")
    # plt.xticks(rotation=45, fontsize=7)
    # plt.show()

    # tss_enrich_vs_repr(loci_2_c, tss, num_bins=4)

def benchmark(replicate_1_dir, replicate_2_dir, pltsavedir):
    t0 = datetime.now()
    loci_1, loci_2 = intersect_parsed_posteriors(
        replicate_1_dir+"/parsed_posterior_short.csv", 
        replicate_2_dir+"/parsed_posterior_short.csv")
    
    print("took {} to load and intersect the data".format(datetime.now() - t0))

    num_labels = loci_1.shape[1]-3
    loci_1.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]
    loci_2.columns = ["chr", "start", "end"]+["posterior{}".format(i) for i in range(num_labels)]

    print('generating confmat 1')
    
    conf_mat = confusion_matrix(
        loci_1, loci_2, num_labels, 
        OE_transform=True, symmetric=False)

    assignment_pairs = Hungarian_algorithm(conf_mat, conf_or_dis='conf')

    if "concat" in replicate_1_dir or "concat" in replicate_2_dir:
        loci_1_mnemon = read_mnemonics(replicate_1_dir+"/mnemonics_rep1.txt")
        loci_2_mnemon = read_mnemonics(replicate_2_dir+"/mnemonics_rep2.txt")
    else:
        loci_1_mnemon = read_mnemonics(replicate_1_dir+"/mnemonics.txt")
        loci_2_mnemon = read_mnemonics(replicate_2_dir+"/mnemonics.txt")
    
    mnemon1_dict = {}
    for i in loci_1_mnemon:
        mnemon1_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:5]

    mnemon2_dict = {}
    for i in loci_2_mnemon:
        mnemon2_dict[i.split("_")[0]] = i.split("_")[0]+'_'+i.split("_")[1][:5]
    
    for i in range(len(assignment_pairs)):
        assignment_pairs[i] = (mnemon1_dict[str(assignment_pairs[i][0])], mnemon2_dict[str(assignment_pairs[i][1])])
    print(assignment_pairs)

    loci_1, loci_2 = \
        connect_bipartite(loci_1, loci_2, assignment_pairs, mnemon=True)
    
    print('connected barpartite')

    print(loci_1)
    print(loci_2)

    ####

    if os.path.exists(pltsavedir) == False:
        os.mkdir(pltsavedir)

    if os.path.exists(pltsavedir+"/cc") == False:
        os.mkdir(pltsavedir+"/cc")

    if os.path.exists(pltsavedir+"/agr") == False:
        os.mkdir(pltsavedir+"/agr")
    
    if os.path.exists(pltsavedir+"/snk") == False:
        os.mkdir(pltsavedir+"/snk")

    if os.path.exists(pltsavedir+"/clb") == False:
        os.mkdir(pltsavedir+"/clb")

    cc = correspondence_curve(loci_1, loci_2, pltsavedir+"/cc")
    cc.plot_curve(plot_general=False, merge_plots=False)
    del cc

    calb = posterior_calibration(
        loci_1, loci_2, log_transform=False, ignore_overconf=False, filter_nan=True, 
        oe_transform=True, savedir=pltsavedir+"/clb")
    calibrated_loci_1 = calb.perlabel_calibration_function(degree=5, num_bins=25, return_caliberated_matrix=True)
    print(calibrated_loci_1)
    
    agr = Agreement(loci_1, loci_2, pltsavedir+"/agr")
    print('per-label agreement:    ', agr.per_label_agreement())
    print('general agreement:    ', agr.general_agreement())
    print('general log(o/e) agreement:    ', agr.general_OE_ratio())
    print('general Cohens Kappa score:    ', agr.general_cohens_kappa())
    agr.plot_agreement()
    agr.plot_CK()
    agr.plot_OE()
    del agr

    vis = sankey(loci_1, loci_2, pltsavedir+"/snk")
    vis.sankey_diag()
    vis.heatmap()
    del vis
    
    print("the whole process took:  ",datetime.now() - t0)
    
if __name__=="__main__":
    replicate_1_dir = "tests/chmm/GM12878_rep1/parsed_posterior.csv"
    replicate_2_dir = "tests/chmm/GM12878_rep2/parsed_posterior.csv"

    full_reproducibility_report(replicate_1_dir, replicate_2_dir, "tests/chmmgm12878/", run_on_subset=True)

    # replicate_1_dir = "tests/segway/GM12878_rep1/parsed_posterior_short.csv"
    # replicate_2_dir = "tests/segway/GM12878_rep2/parsed_posterior_short.csv"

    # full_reproducibility_report(replicate_1_dir, replicate_2_dir, "tests/pltsavedir_segw/")
    

    exit()
    # replicate_1_dir = 'tests/chromhmm_runs/gm12878_rep1/parsed_posterior_short.csv'
    # replicate_2_dir = 'tests/chromhmm_runs/gm12878_rep2/parsed_posterior_short.csv'
    # test_calibration(replicate_1_dir, replicate_2_dir)

    # replicate_1_dir = 'tests/segway_runs/gm12878_rep1/parsed_posterior_short.csv'
    # replicate_2_dir = 'tests/segway_runs/gm12878_rep2/parsed_posterior_short.csv'
    # test_calibration(replicate_1_dir, replicate_2_dir)
    # benchmark(replicate_1_dir, replicate_2_dir,"tests/reprod_segw")

    # replicate_1_dir = 'tests/chromhmm_runs/GM12878_rep1/'
    # replicate_2_dir = 'tests/chromhmm_runs/GM12878_rep2/'
    # benchmark(replicate_1_dir, replicate_2_dir,  "tests/reprod_chmm")

    # full_reproducibility_report(replicate_1_dir, replicate_2_dir, "TEST_REPR") 
    
    # replicate_1_dir = 'tests/chmm/GM12878_rep1/parsed_posterior.csv'
    # replicate_2_dir = 'tests/chmm/GM12878_rep2/parsed_posterior.csv'
    # loci_1, loci_2 = intersect_parsed_posteriors(replicate_1_dir, replicate_2_dir)

    # loci_1 = loci_1.iloc[:500000,:]
    # loci_2 = loci_2.iloc[:500000,:]

    # tss = pd.read_csv("tests/RefSeqTSS.hg38.txt", sep="\t", header=None)
    # tss.columns = ["chr", "coord", "strand"]
    # tss.insert(2, "end", tss["start"].astype(int)+1)
    # loci_1.start = loci_1.start.astype("int")
    # loci_2.start = loci_2.start.astype("int")
    # loci_1.end = loci_1.end.astype("int")
    # loci_2.end = loci_2.end.astype("int")
    # # tss.start = tss.start.astype("int")
    # # tss.end = tss.end.astype("int")
    # print(loci_1)
    
    # if os.path.exists("tests/clb_1_logoe") == False:
    #     os.mkdir("tests/clb_1_logoe")
    # if os.path.exists("tests/clb_2_logoe") == False:
    #     os.mkdir("tests/clb_2_logoe")

    # if os.path.exists("tests/clb_1_oe") == False:
    #     os.mkdir("tests/clb_1_oe")
    # if os.path.exists("tests/clb_2_oe") == False:
    #     os.mkdir("tests/clb_2_oe")
    
    # TSS_test(loci_1, loci_2, tss, replicate_1_dir, replicate_2_dir)
    

