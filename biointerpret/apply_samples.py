# for each of your samples
import pandas as pd
import util
import sys, os, gzip, pickle

def get_mnemon(source_dir, output_dir):
    # the original 16 features that go to the classifier - just keeping it here
    feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

    mapping_file = "{}/trackname_assay.txt".format(source_dir)

    sample_list = [ct for ct in os.listdir(source_dir) if os.path.isdir(source_dir+ct)]

    #iterate through all folders in source_dir
    
    track_assay_map = {}
    inputTrack_list = []
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[1]] = fields[1]
            inputTrack_list.append(fields[1])

    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    for s in sample_list:
        ann_features, ann_label_bases, ann_feature_names = util.features_from_segtools_dir(
            "{}/{}/feature_aggregation.tab".format(source_dir, s), 
            "{}/{}/signal_distribution.tab".format(source_dir, s), 
            track_assay_map)

        df = pd.DataFrame(ann_features)
        dft = df.T
        dftr = dft[feature_names]

        # load the model
        model_file = "model.pickle.gz"
        with gzip.open(model_file, "r") as f:
            the_model = pickle.load(f)

        labels = the_model.predict(dftr)
        old_labels = dft.index

        if os.path.exists(output_dir+ "/" + s) == False:
            os.mkdir(output_dir+ "/" + s)
            
        with open(output_dir+ "/" + s+ "/mnemonics.txt", 'w') as f:
            f.write("old\tnew\n")
            for i in range(len(old_labels)):
                f.write("{}\t{}\n".format(old_labels[i], labels[i]))


if __name__ == "__main__":
    source_dir = sys.argv[1]
    output_dir = sys.argv[2]
    get_mnemon(source_dir, output_dir)