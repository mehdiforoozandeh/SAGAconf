#!/bin/env python

##################################
# begin header
import sys
import os
import argparse
import subprocess
import gzip
import math
import random
from path import Path
from collections import defaultdict
#import bedtools
import shutil
import os

sys.path.append(".")
import util

parser = argparse.ArgumentParser()
parser.add_argument("workdir", type=Path)
parser.add_argument("--model-path", type=Path)
parser.add_argument("--input-path", type=Path)
args = parser.parse_args()
workdir = args.workdir

if not workdir.exists():
    workdir.makedirs()

shutil.copy(__file__, workdir / "apply_samples.py")
experiment_dir = Path(os.getcwd())
project_dir = experiment_dir / ".." / ".."
os.chdir(workdir)
workdir = Path(".")

import logging
logging.basicConfig(format="%(asctime)s %(levelname)s:%(message)s")
logger = logging.getLogger('log')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler("stdout.txt")
fh.setLevel(logging.DEBUG) # >> this determines the file level
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)# >> this determines the output level
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s')
ch.setFormatter(formatter)
fh.setFormatter(formatter)
# add the handlers to logger
logger.addHandler(ch)
logger.addHandler(fh)

# end header
##################################

from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.ensemble import RandomForestClassifier
import numpy
import pandas
import pickle
from copy import deepcopy

import psutil
import resource
process = psutil.Process(os.getpid())
def log_mem():
    psutil_mem = float(process.memory_info().rss) / float(2 ** 30)
    resource_mem = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / float(2 ** 20)
    logger.info("Memory usage -- resource: %.1f gb; psutil: %.1f gb", resource_mem, psutil_mem )

#######################################################
# Input files
#######################################################

# Trained sklearn model for classification
model_fname = str(args.model_path) if args.model_path is not None else "../model.pickle.gz" # model in the experimentdir
model = pickle.load(gzip.open(model_fname))

segtools_dir = args.input_path or Path('../segwayOutput') # contains a directory for each cell type, where segtool outputs are

trackname_mapping_fname = segtools_dir / 'trackname_assay.txt'

segtools_trackname_mapping = {}
with open(trackname_mapping_fname, "r") as f:
     for line in f:
         line = line.split()
         segtools_trackname_mapping[line[0]] = line[1]

print('1', segtools_trackname_mapping) ###testprint

########################################################
# For Plotting
#######################################################
feature_order_R = """c('H3K9me3', 'H3K27me3', 'H3K36me3', 'H3K4me3', 'H3K4me1', 'H3K27ac', "5' flanking (1000-10000 bp)", "5' flanking (1-1000 bp)", 'initial exon', 'initial intron', 'internal exons', 'internal introns', 'terminal exon', 'terminal intron', "3' flanking (1-1000 bp)", "3' flanking (1000-10000 bp)")"""

bio_label_order_R =  """c('Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence')"""
bio_label_order_ref_R =  """c('Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'Unclassified')"""


###################### Choose one of the following ################
# This one for colour sclae same as the paper with scales predefined
feature_heatmap_disc_vals_R = "c(-Inf, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, Inf)"
# This one for variable colour sclae containing 11 colours
#feature_heatmap_disc_vals_R =  11


#######################################################
# Get target anns
#######################################################
# for each annotation file, the required info for the classifier is the address of the signal_distribution and feature_aggregation files, given in Path objects.
# The other fields are meta information and include: ann_id, celltype, dataset_key, url, assembly, tool (for us it is all segway) - for the test run I only filled the celltype and the ann_id

target_anns = []
next_id = 0
# get a list of cell-types or annotations

celltypedirs = [name for name in os.listdir(segtools_dir) if os.path.isdir(os.path.join(segtools_dir, name))]
print('2', segtools_dir, celltypedirs)###testprint
for celltype in celltypedirs:
    signal_dist_fname = segtools_dir / celltype / "signal_distribution.tab"
    gene_agg_fname = segtools_dir / celltype / "feature_aggregation.tab"
    ann = {"ann_id": next_id, "url": '', "celltype": celltype, "tool": '', "dataset_key": '', 
           "concatenation_key": '', "assembly": '', 
           "signal_dist_fname": signal_dist_fname, "gene_agg_fname": gene_agg_fname}
    next_id += 1
    target_anns.append(ann)
print('3', target_anns)###testprint
#######################################################
# Classify our annotations
#######################################################
logger.info("Starting classifying our annotation labels...")
log_mem()

################## Providing feature order for the classifier - same as training  #####################
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

num_anns = len(target_anns)
model_classes_dict = {label: i for i, label in enumerate(model.classes_)}
print('4', model_classes_dict) ###testprint
classification_dir = Path("classification")
if not classification_dir.exists():
    classification_dir.makedirs()

bio_labels = {} # {ann_id: {int_label: bio_label}}
bio_label_probs = {} # {ann_id: {int_label: {bio_label: prob}}}

for ann_index, ann in enumerate(target_anns):
    print('5', target_anns[ann_index]['celltype'], ann_index) ###testprint
    ann_id = ann["ann_id"]

    ann_features, _, feature_list = util.features_from_segtools_dir(ann["gene_agg_fname"], ann["signal_dist_fname"], segtools_trackname_mapping)

    print('6', ann_features) ###testprint
    print('7', feature_list) ###testprint

    if (len(feature_list) < 16): # if histone tracks are missing
        logger.info("Histone tracks are missing, skipping.".format(**locals()))
        # continue

    target_dir = classification_dir / str(target_anns[ann_index]['celltype'])
    print('8', target_dir) ###testprint
    if not target_dir.exists():
        target_dir.makedirs()
    
    print('9', 'made target dir', target_dir) ###testprint
    mnemonics_outfn = target_dir / "mnemonics.txt"
    classifier_data_outfn = target_dir / "classifier_data.tab"

    print('10', mnemonics_outfn) ###testprint
    print('11', classifier_data_outfn) ###testprint

    if mnemonics_outfn.exists():
        logger.info("Found mnemonics for  {ann_id}, skipping.".format(**locals()))
    else:
        logger.info("Starting ann {ann_index}/{num_anns} {ann_id}...".format(**locals()))
        
        labels = list(ann_features.keys())
        c = 0
        with open(classifier_data_outfn, "w") as f:
            f.write("orig_label\t")
            f.write("\t".join(feature_names))
            f.write("\n")
            for label in ann_features:
                print(c, label)
                f.write(label)
                for feature_name in feature_names:
                    f.write("\t")
                    f.write(str(ann_features[label][feature_name]))
                f.write("\n")
                c += 1
        classifier_data_frame = pandas.read_csv(classifier_data_outfn, sep="\t")
        features_mat = util.features_frame_to_matrix(classifier_data_frame, feature_names)
        if numpy.any(numpy.isnan(features_mat)):
            logger.error("Missing features for ann {ann_id}!".format(**locals()))
            raise Exception("Missing features for ann {ann_id}!".format(**locals()))
            predicted_labels = ["FeaturesMissing" for i in range(len(labels))]
            probs = numpy.ones(shape=(len(labels), len(model.classes_))) * (float(1)/len(model.classes_))
        else:
            predicted_labels = model.predict(features_mat)
            probs = model.predict_proba(features_mat)
        if not target_dir.exists(): target_dir.makedirs()
        bio_labels[ann_id] = {}
        with open(mnemonics_outfn, "w") as mnemonics_out:
            mnemonics_out.write("old\tnew\n")
            for i in range(len(labels)):
                label = labels[i]
                predicted_label = predicted_labels[i]
                if ((predicted_label == "Unclassified")
                    or ((predicted_label != "FeaturesMissing")
                        and (probs[i, model_classes_dict[predicted_label]] <= 0.25))):
                    predicted_label = "LowConfidence"
                bio_labels[ann_id][label] = predicted_label
                mnemonics_out.write("{label}\t{predicted_label}\n".format(**locals()))
        bio_label_probs[ann_id] = {}
        probs_outfn = target_dir / "probs.txt"
        with open(probs_outfn, "w") as f:
            f.write("label\t" + "\t".join(map(str, model.classes_)) + "\n")
            for i, label in enumerate(labels):
                bio_label_probs[ann_id][label] = {}
                for bio_label_index, bio_label in enumerate(model.classes_):
                    bio_label_probs[ann_id][label][bio_label] = probs[i,bio_label_index]
                f.write(str(label))
                f.write("\t")
                f.write("\t".join(map(str, probs[i,:])))
                f.write("\n")

exit()

##################################
##################################
# Diagnostic plots
##################################
##################################

##################################
# Read bio label assignments
##################################

bio_label_dir = classification_dir
print("bio_label_dir", bio_label_dir) ###testprint
bio_labels = {} # {celltype: {int_label: bio_label}}

for celltype in bio_label_dir.listdir():
    print('bio_label_dir.listdir()', bio_label_dir.listdir()) ###testprint
    print('celltype', celltype)###testprint
    celltype = celltype.basename()
    print('celltypebasename', celltype)###testprint
    mnem_fname = bio_label_dir / celltype / "mnemonics.txt"
    print('mnem_fname', mnem_fname)###testprint
    if mnem_fname.exists():
        print('mnem_fname.exists()')###testprint
        celltype_bio_labels = {}
        features_missing = False
        with open(mnem_fname, "r") as f:
            for line in f:
                if "old" in line: continue
                line = line.split()
                int_label = line[0]
                bio_label = line[1]
                if bio_label == "FeaturesMissing":
                    features_missing = True
                celltype_bio_labels[int_label] = bio_label
        if not features_missing:
            bio_labels[str(celltype)] = celltype_bio_labels

#bio_labels = {celltype: bio_labels[celltype] for celltype in bio_labels if bio_labels[celltype]['0'] != "FeaturesMissing"}
bio_labels = {celltype: bio_labels[celltype] for celltype in bio_labels}  

all_bio_labels = set.union(*map(set, map(lambda x: x.values(), bio_labels.values())))

############################
# Read probs assignment
############################
bio_label_probs = {}
for celltype in bio_label_dir.listdir():
    celltype = celltype.basename()
    probs_fname = bio_label_dir / celltype / "probs.txt"
    if probs_fname.exists(): 
        celltype_bio_labels_prob = {}
        with open(probs_fname, "r") as f:
            for line in f:
                header = line
                header = header.split()
                break
            for line in f:
                if "label" in line: continue
                line = line.split()
                int_label = line[0]
                temp = {}
                for i in range(len(header) - 1):
                    temp[header[i+1]] = line[i+1]
                celltype_bio_labels_prob[int_label] = temp
        bio_label_probs[str(celltype)] = celltype_bio_labels_prob


bio_labels_outfn = workdir / "bio_labels.pickle.gz"
with gzip.open("bio_labels.pickle.gz", "w") as f:
    pickle.dump(bio_labels, f)

bio_label_prob_outfn = workdir / "bio_label_probs.pickle.gz"
with gzip.open("bio_label_probs.pickle.gz", "w") as f:
    pickle.dump(bio_label_probs, f)

#######################################################
# Distribution of label probs
#######################################################
log_mem()

if not Path("stats").exists():
    Path("stats").makedirs()

if True:
    key = "stats/model_predicted_prob"

    data_fname = workdir / "{key}.tab".format(**locals())
    with open(data_fname, "w") as f:
        f.write("celltype\tint_label\tprob\n")
        for celltype_index, celltype in enumerate(bio_labels):
            for int_label in bio_label_probs[celltype]:
                prob = max(bio_label_probs[celltype][int_label].values())
                f.write("{celltype}\t{int_label}\t{prob}\n".format(**locals()))
                
    script_fname = workdir / "{key}.R".format(**locals())
    script = \
    """
    require(Cairo)
    require(ggplot2)
    require(reshape2)
    require(RColorBrewer)

    data <- read.delim("{data_fname}", header=TRUE)

    p <- ggplot(data) +
      aes(x=prob) +
      stat_ecdf() +
      scale_x_continuous(limits=c(0,1)) +
      theme_bw()
    ggsave(paste("{key}_ecdf.pdf", sep=""), p, width=10, height=6, units="in")
    """.format(**locals())

    with open(script_fname, "w") as f:
        f.write(script)

    cmd = ["Rscript", script_fname]
    logger.info(" ".join(cmd))
    subprocess.check_call(cmd)

log_mem()
##################################
# Read length distributions
##################################
log_mem()

length_dist_dir = segtools_dir
length_dists = {}
for celltype in length_dist_dir.listdir():
    celltype = celltype.basename()
    if celltype in ["run.py", "stdout.txt", "trackname_assay.txt"]: continue
    length_dist_fname = length_dist_dir / celltype / "segment_sizes.tab"
    if length_dist_fname.exists():
        print(celltype)
        length_dists[celltype] = pandas.read_table(length_dist_fname)


##################################
# common cell types
##################################
common_celltypes = set(length_dists.keys()) & set(bio_labels.keys())
print(common_celltypes)

#######################################################
# Plot indiv features
#######################################################
log_mem()

############### Normalized Features ##################

if not Path("indiv").exists():
    Path("indiv").makedirs()

if not Path("indiv/Normalized").exists():
    Path("indiv/Normalized").makedirs()

if not Path("indiv/Normalized/Bio-labels").exists():
    Path("indiv/Normalized/Bio-labels").makedirs()

if not Path("indiv/Normalized/Clusters").exists():
    Path("indiv/Normalized/Clusters").makedirs()


feature_names = set()
key = "indiv/Normalized/classification_features_normalized"
key_bio = "indiv/Normalized/Bio-labels/classification_features_normalized"
key_cluster = "indiv/Normalized/Clusters/classification_features_normalized"
data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("celltype\tfeature_name\tint_label\tbio_label\tlabel_key\tfeature_val\n")
    for celltype_index, celltype in enumerate(common_celltypes):
        if not celltype in bio_labels: continue
        print(celltype)
        target_dir = classification_dir / celltype
        classifier_data_outfn = target_dir / "classifier_data.tab"
        classifier_data_frame = pandas.read_csv(classifier_data_outfn, sep="\t")
        for col in list(classifier_data_frame.columns):
            if col == "orig_label":
                continue
            classifier_data_frame[col] = (classifier_data_frame[col] - classifier_data_frame[col].mean())/classifier_data_frame[col].std()
        for row_index in range(classifier_data_frame.shape[0]):
            int_label = classifier_data_frame.loc[row_index, "orig_label"]
            int_label = str(int(int_label))
            bio_label = bio_labels[celltype][int_label]
            label_key = "{bio_label}_{int_label}".format(**locals())
            for feature_name in classifier_data_frame.columns:
                if feature_name == "orig_label": continue
                feature_names.add(feature_name)
                feature_val = classifier_data_frame.loc[row_index, feature_name]
                f.write("{celltype}\t{feature_name}\t{int_label}\t{bio_label}\t{label_key}\t{feature_val}\n".format(**locals()))
feature_names = list(feature_names)
#feature_heatmap_disc_vals_R = 11 #Max upto 11 colours in RdBu
script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
require(reshape2)
require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)

data$feature_val_disc <- cut(data$feature_val, breaks={feature_heatmap_disc_vals_R})
palette="RdBu"
colors <- brewer.pal(palette, n=length(levels(data$feature_val_disc)))
color_mapping <- length(levels(data$feature_val_disc))
colors <- rev(colors)
data$label_celltype_key <- as.factor(paste(data$celltype, "_", data$int_label, "___", data$bio_label, sep=""))

p <- ggplot(data) +
  aes(x=feature_val) +
  stat_ecdf() +
  scale_x_continuous(breaks=seq(-2,5,0.5), limits=c(-3,5)) +
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste("{key}_ecdf.pdf", sep=""), p, width=10, height=6, units="in")

common_celltypes = unique(data$celltype)
for (i in common_celltypes) {{
print(i)
p <- ggplot(data[data$celltype==i,], aes(x=label_key, y=feature_name, fill=feature_val_disc)) +
  geom_tile() +
  ggtitle(i) +
  scale_fill_manual(values=colors) +
  scale_y_discrete(name="Feature") +
  scale_x_discrete(name="Label") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste("{key_cluster}_",i,".pdf", sep=""), p)
}}

for (i in common_celltypes) {{
print(i)
p <- ggplot(data[data$celltype==i,], aes(x=bio_label, y=feature_name, fill=feature_val_disc)) +
  geom_tile() +
  ggtitle(i) +
  scale_fill_manual(values=colors) +
  scale_y_discrete(name="Feature") +
  scale_x_discrete(name="Label") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste("{key_bio}_",i,".pdf", sep=""), p)
}}

label_key_order <- c()
cast_data <- acast(data, label_celltype_key ~ feature_name, value.var="feature_val")
ord <- hclust( dist(cast_data, method = "euclidean") )$order
data$label_celltype_key <- factor(data$label_celltype_key, levels=rownames(cast_data)[ord])

p <- ggplot(data, aes(y=label_celltype_key, x=feature_name, fill=feature_val_disc)) +
  geom_tile() +
  scale_fill_manual(values=colors, breaks=color_mapping, drop=FALSE) +
  scale_x_discrete(name="Feature") +
  scale_y_discrete(name="Label") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("{key}_all.pdf", p, width=15, height=350, units="in", limitsize=FALSE)
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)
cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)

log_mem()

############### Non-Normalized Features ##################

if not Path("indiv").exists():
    Path("indiv").makedirs()

if not Path("indiv/Not-Normalized").exists():
    Path("indiv/Not-Normalized").makedirs()

if not Path("indiv/Not-Normalized/Bio-labels").exists():
    Path("indiv/Not-Normalized/Bio-labels").makedirs()

if not Path("indiv/Not-Normalized/Clusters").exists():
    Path("indiv/Not-Normalized/Clusters").makedirs()

feature_names = set()
key = "indiv/Not-Normalized/classification_features"
key_bio = "indiv/Not-Normalized/Bio-labels/classification_features_normalized"
key_cluster = "indiv/Not-Normalized/Clusters/classification_features_normalize"
data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("celltype\tfeature_name\tint_label\tbio_label\tlabel_key\tfeature_val\n")
    for celltype_index, celltype in enumerate(common_celltypes):
        if not celltype in bio_labels: continue
        print(celltype)
        target_dir = classification_dir / celltype
        classifier_data_outfn = target_dir / "classifier_data.tab"
        classifier_data_frame = pandas.read_csv(classifier_data_outfn, sep="\t")
        for row_index in range(classifier_data_frame.shape[0]):
            int_label = classifier_data_frame.loc[row_index, "orig_label"]
            int_label = str(int(int_label))
            bio_label = bio_labels[celltype][int_label]
            label_key = "{bio_label}_{int_label}".format(**locals())
            for feature_name in classifier_data_frame.columns:
                if feature_name == "orig_label": continue
                feature_names.add(feature_name)
                feature_val = classifier_data_frame.loc[row_index, feature_name]
                f.write("{celltype}\t{feature_name}\t{int_label}\t{bio_label}\t{label_key}\t{feature_val}\n".format(**locals()))
feature_names = list(feature_names)
common_celltypes_list = set(common_celltypes)
#feature_heatmap_disc_vals_R = 11 #Max upto 11 colours in RdBu
script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
require(reshape2)
require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)

data$feature_val_disc <- cut(data$feature_val, breaks={feature_heatmap_disc_vals_R})
palette="RdBu"
colors <- brewer.pal(palette, n=length(levels(data$feature_val_disc)))
color_mapping <- length(levels(data$feature_val_disc))
colors <- rev(colors)
data$label_celltype_key <- as.factor(paste(data$celltype, "_", data$int_label, "___", data$bio_label, sep=""))

p <- ggplot(data) +
  aes(x=feature_val) +
  stat_ecdf() +
  scale_x_continuous(breaks=seq(-2,5,0.5), limits=c(-3,5)) +
  scale_y_continuous(breaks=seq(0,1,0.1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste("{key}_ecdf.pdf", sep=""), p, width=10, height=6, units="in")

common_celltypes = unique(data$celltype)
for (i in common_celltypes) {{
print(i)
p <- ggplot(data[data$celltype==i,], aes(x=label_key, y=feature_name, fill=feature_val_disc)) +
  geom_tile() +
  ggtitle(i) +
  scale_fill_manual(values=colors) +
  scale_y_discrete(name="Feature") +
  scale_x_discrete(name="Label") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste("{key_cluster}_",i,".pdf", sep=""), p)
}}

for (i in common_celltypes) {{
print(i)
p <- ggplot(data[data$celltype==i,], aes(x=bio_label, y=feature_name, fill=feature_val_disc)) +
  geom_tile() +
  ggtitle(i) +
  scale_fill_manual(values=colors) +
  scale_y_discrete(name="Feature") +
  scale_x_discrete(name="Label") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste("{key_bio}_",i,".pdf", sep=""), p)
}}

label_key_order <- c()
cast_data <- acast(data, label_celltype_key ~ feature_name, value.var="feature_val")
ord <- hclust( dist(cast_data, method = "euclidean") )$order
data$label_celltype_key <- factor(data$label_celltype_key, levels=rownames(cast_data)[ord])

p <- ggplot(data, aes(y=label_celltype_key, x=feature_name, fill=feature_val_disc)) +
  geom_tile() +
  scale_fill_manual(values=colors, breaks=color_mapping, drop=FALSE) +
  scale_x_discrete(name="Feature") +
  scale_y_discrete(name="Label") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("{key}_all.pdf", p, width=15, height=350, units="in", limitsize=FALSE)
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)
cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)

log_mem()

##################################
# Length Distributions
##################################

if not Path("Length_distribution").exists():
    Path("Length_distribution").makedirs()

################## Cluster Length Distributions ###############
if not Path("Length_distribution/Clusters").exists():
    Path("Length_distribution/Clusters").makedirs()

for key in common_celltypes:
	length_distr_fname = length_dist_dir / key / "length_distribution.tab"
	length_distr_df = pandas.read_table(length_distr_fname)

	clusters = length_distr_df.label.unique()

	new_df = length_distr_df.groupby(by='label').count()
	new_df['>=1e4'] = length_distr_df[length_distr_df['length'] >= 10000].groupby(by='label').count()
	new_df['=100'] = length_distr_df[length_distr_df['length'] == 100].groupby(by='label').count()
	new_df['<=200'] = length_distr_df[length_distr_df['length'] <= 200].groupby(by='label').count()

	print(new_df)

	with open("Length_distribution/Clusters/{key}_outliers_table.tab".format(**locals()), "w") as f:
		f.write(str(new_df))

	script_fname = "Length_distribution/Clusters/Length_dist_clusters.R".format(**locals())
	script = \
	"""
	require(Cairo)
	require(ggplot2)

	data <- read.delim("{length_distr_fname}", header=TRUE)

	p <- ggplot(data[data$length <= 3000 ,], aes(x=length)) +
	  geom_histogram(binwidth=100, colour="black", fill="white") + 
	  #geom_density() +
	  facet_wrap( ~ label ,scales = "free")
	ggsave("Length_distribution/Clusters/{key}_hist.pdf", p)
	""".format(**locals())

	with open(script_fname, "w") as f:
	    f.write(script)

	cmd = ["Rscript", script_fname]
	subprocess.check_call(cmd)

###################### Bio Labels Length Distribution ##############

if not Path("Length_distribution/Bio_labels").exists():
    Path("Length_distribution/Bio_labels").makedirs()

#with gzip.open("bio_labels.pickle.gz", "r") as f:
		#bio_labels = pickle.load(f)

for key in common_celltypes:
	length_distr_fname = length_dist_dir / key / "length_distribution.tab"
	length_distr_df = pandas.read_table(length_distr_fname)

	length_distr_df['label'] = length_distr_df['label'].astype(str)
	length_distr_df['bio_label'] = length_distr_df['label'].map(bio_labels[key])
	length_distr_df.to_csv('Length_distribution/Bio_labels/length_distr_df_temp.csv')

	new_df = length_distr_df.drop(['label'], axis=1, inplace=True)
	new_df = length_distr_df.groupby(by='bio_label').count()
	
	print(new_df)
	
	new_df['>=1e4'] = length_distr_df[length_distr_df['length'] >= 10000].groupby(by='bio_label').count()
	new_df['=100'] = length_distr_df[length_distr_df['length'] == 100].groupby(by='bio_label').count()
	new_df['<=200'] = length_distr_df[length_distr_df['length'] <= 200].groupby(by='bio_label').count()

	print(new_df)

	with open("Length_distribution/Bio_labels/{key}_outliers_table.tab".format(**locals()), "w") as f:
		f.write(str(new_df))

	script_fname = "Length_distribution/Bio_labels/Length_dist_bio-label.R".format(**locals())
	script = \
	"""
	require(Cairo)
	require(ggplot2)

	data <- read.csv("Length_distribution/Bio_labels/length_distr_df_temp.csv", header=TRUE)

	p <- ggplot(data[data$length <= 3000 ,], aes(x=length)) +
	  geom_histogram(binwidth=100, colour="black", fill="white") + 
	  #geom_density() +
	  facet_wrap( ~ bio_label ,scales = "free")
	ggsave("Length_distribution/Bio_labels/{key}_hist.pdf", p)
	""".format(**locals())

	with open(script_fname, "w") as f:
		f.write(script)

	cmd = ["Rscript", script_fname]
	subprocess.check_call(cmd)

##################################
# num bp ~ bio label
##################################
print("Plotting num bp ~ bio label")
log_mem()
if not Path("stats").exists():
    Path("stats").makedirs()

key = "stats/num_bp"
num_bp = {}
for celltype in common_celltypes:
    for row_id in length_dists[celltype].index:
        int_label = length_dists[celltype].loc[row_id, "label"]
        if int_label == "all": continue
        bio_label = bio_labels[celltype][int_label]
        if not bio_label in num_bp:
            num_bp[bio_label] = 0
        num_bp[bio_label] += length_dists[celltype].loc[row_id, "num.bp"]

data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("bio_label\tnum_bp\n")
    for bio_label in num_bp:
        num_bp_val = num_bp[bio_label]
        f.write("{bio_label}\t{num_bp_val}\n".format(**locals()))

#data$bio_label <- ordered(data$bio_label, {bio_label_order_R})
script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
#require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)

data$frac_bp <- data$num_bp / sum(data$num_bp)
print(head(data))
data$bio_label <- ordered(data$bio_label, {bio_label_order_R})
p <- ggplot(data, aes(x=bio_label, y=frac_bp)) +
  geom_bar(stat="identity") +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Fraction of genome") +
  coord_flip() +
  theme_classic()
ggsave("{key}.pdf", p, width=4, height=2.5, units="in")

p <- ggplot(data, aes(x=bio_label, y=frac_bp)) +
  geom_point() +
  scale_y_log10(breaks=c(0.001,0.003,0.01,0.03,0.1,0.3), name="Fraction of bases") +
  scale_x_discrete(name="") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggsave("{key}_log.pdf", p, width=4, height=5, units="in")
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)

print("Plot finished")

##################################
# Mean observed classifier feature values for each label
##################################
print("Plotting Mean observed classifier feature values for each label")
log_mem()

key = "stats/mean_feature_vals"
data_fname = workdir / "{key}.tab".format(**locals())

mean_features = {bio_label: None for bio_label in all_bio_labels}
mean_features_num_bp = {bio_label: 0 for bio_label in all_bio_labels}
for celltype in common_celltypes:
    logger.info("Starting reading features for celltype {celltype}...".format(**locals()))
    target_dir = classification_dir / celltype
    classifier_data_outfn = target_dir / "classifier_data.tab"
    classifier_data_frame = pandas.read_csv(classifier_data_outfn, sep="\t")
    features_mat = util.features_frame_to_matrix(classifier_data_frame, feature_names)
    for example_index in range(classifier_data_frame.shape[0]):
        int_label = str(classifier_data_frame.orig_label[example_index])
        bio_label = bio_labels[celltype][int_label]
        label_num_bp = length_dists[celltype].loc[numpy.where(length_dists[celltype].label == int_label)[0][0], "num.bp"]
        if mean_features_num_bp[bio_label] == 0:
            mean_features[bio_label] = features_mat[example_index,:]
        else:
            mean_features[bio_label] = (mean_features[bio_label]*mean_features_num_bp[bio_label] + features_mat[example_index,:]*label_num_bp) / (mean_features_num_bp[bio_label] + label_num_bp)
        mean_features_num_bp[bio_label] = label_num_bp

with open(data_fname, "w") as f:
    f.write("bio_label\tfeature_name\tfeature_mean\n")
    for bio_label in all_bio_labels:
        if not (mean_features[bio_label] is None):
            for feature_index, feature_name in enumerate(feature_names):
                feature_mean = mean_features[bio_label][feature_index]
                f.write("{bio_label}\t{feature_name}\t{feature_mean}\n".format(**locals()))

#data$bio_label <- ordered(data$bio_label, {bio_label_order_R})
script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)

data$bio_label <- ordered(data$bio_label, {bio_label_order_R})
data$feature_mean_disc <- cut(data$feature_mean, breaks={feature_heatmap_disc_vals_R})
palette="RdBu"
colors <- brewer.pal(palette, n=length(levels(data$feature_mean_disc)))
color_mapping <- levels(data$feature_mean_disc)
colors <- rev(colors)

p <- ggplot(data) +
  aes(x=bio_label, y=feature_name, fill=feature_mean_disc) +
  geom_tile() +
  scale_fill_manual(values=colors, breaks=color_mapping, drop=FALSE, name="Feature\nvalue") +
  scale_y_discrete(name="") +
  scale_x_discrete(name="") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("{key}.pdf", p, width=6, height=4.5, units="in")
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)

print("Plot finished")

##################################
# num_assignments ~ bio label
##################################
print("Plotting num_assignments ~ bio label")
log_mem()
key = "stats/num_assignments"
num_assignments = {}
for celltype in common_celltypes:
    for row_id in length_dists[celltype].index:
        int_label = length_dists[celltype].loc[row_id, "label"]
        if int_label == "all": continue
        bio_label = bio_labels[celltype][int_label]
        if not bio_label in num_assignments:
            num_assignments[bio_label] = 0
        num_assignments[bio_label] += 1

data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("bio_label\tnum_assignments\tfrac_assignments\n")
    for bio_label in num_assignments:
        num_assignments_val = num_assignments[bio_label]
        frac_assignments_val = float(num_assignments[bio_label]) / len(common_celltypes)
        f.write("{bio_label}\t{num_assignments_val}\t{frac_assignments_val}\n".format(**locals()))

#data$bio_label <- ordered(data$bio_label, {bio_label_order_R})   
script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
#require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)
data$bio_label <- ordered(data$bio_label, {bio_label_order_R})   

p <- ggplot(data, aes(x=bio_label, y=frac_assignments)) +
  geom_bar(stat="identity") +
  scale_y_continuous(name="Average number of\nlabels per cell type") +
  scale_x_discrete(name="") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggsave("{key}.pdf", p, width=4, height=4, units="in")
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)
print("Plot finished")
##################################
# number of celltypes with at least one assignment ~ bio label
##################################
print("Plotting number of celltypes with at least one assignment ~ bio label")
log_mem()
key = "stats/num_celltypes"
num_celltypes = {}
for celltype in common_celltypes:
    celltype_bio_labels = set()
    for row_id in length_dists[celltype].index:
        int_label = length_dists[celltype].loc[row_id, "label"]
        if int_label == "all": continue
        bio_label = bio_labels[celltype][int_label]
        celltype_bio_labels.add(bio_label)
    for bio_label in celltype_bio_labels:
        if not bio_label in num_celltypes:
            num_celltypes[bio_label] = 0
        num_celltypes[bio_label] += 1

data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("bio_label\tnum_celltypes\n")
    for bio_label in num_celltypes:
        num_celltypes_val = float(num_celltypes[bio_label]) / len(common_celltypes)
        f.write("{bio_label}\t{num_celltypes_val}\n".format(**locals()))

#data$bio_label <- ordered(data$bio_label, {bio_label_order_R})    
script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
#require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)
data$bio_label <- ordered(data$bio_label, {bio_label_order_R})   

p <- ggplot(data, aes(x=bio_label, y=num_celltypes)) +
  scale_y_continuous(name="Fraction of cell types with\nat least one label", limits=c(0,1)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw()
ggsave("{key}.pdf", p, width=6, height=5, units="in")
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)

print("Plot finished")

##################################
# meaning ~ distribution of number of labels with that label over CTs
##################################
print("Plotting meaning ~ distribution of number of labels with that label over CTs  ")
log_mem()
key = "stats/num_assignments_per_celltype"
data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("celltype\tbio_label\tnum_assignments\n")
    for celltype_index, celltype in enumerate(bio_labels):
        num_assignments = {bio_label: 0 for bio_label in all_bio_labels}
        for int_label in bio_labels[celltype]:
            bio_label = bio_labels[celltype][int_label]
            num_assignments[bio_label] += 1
        for bio_label in all_bio_labels:
            num = num_assignments[bio_label]
            f.write("{celltype}\t{bio_label}\t{num}\n".format(**locals()))

#data$bio_label <- ordered(data$bio_label, {bio_label_order_R}) 
script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)
data$bio_label <- ordered(data$bio_label, {bio_label_order_R})   

p <- ggplot(data, aes(x=num_assignments)) +
  stat_bin() +
  facet_grid(bio_label ~ .) +
  scale_x_continuous(name="Number of assignments in a given cell type") +
  scale_y_continuous(name="Number of cell types") +
  theme_bw()
ggsave("{key}.pdf", p, width=6, height=8, units="in")
""".format(**locals())
with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))    
subprocess.check_call(cmd)

print("Plot finished")
exit()

