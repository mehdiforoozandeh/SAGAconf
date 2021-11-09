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
args = parser.parse_args()
workdir = args.workdir

if not workdir.exists():
    workdir.makedirs()

shutil.copy(__file__, workdir / "apply.py")
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
model_fname = "model.pickle.gz"
model = pickle.load(gzip.open(model_fname))

segtools_dir = Path('/home/ishangoel04062001/interpretation_final/segwayOutput') # contains a directory for each cell type, where segtool outputs are

trackname_mapping_fname = segtools_dir / 'trackname_assay.txt'

segtools_trackname_mapping = {}
with open(trackname_mapping_fname, "r") as f:
     for line in f:
         line = line.split()
         segtools_trackname_mapping[line[0]] = line[1]

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
for celltype in celltypedirs:
    signal_dist_fname = segtools_dir / celltype / "signal_distribution.tab"
    gene_agg_fname = segtools_dir / celltype / "feature_aggregation.tab"
    ann = {"ann_id": next_id, "url": '', "celltype": celltype, "tool": '', "dataset_key": '', 
           "concatenation_key": '', "assembly": '', 
           "signal_dist_fname": signal_dist_fname, "gene_agg_fname": gene_agg_fname}
    next_id += 1
    target_anns.append(ann)

##################################
##################################
# Diagnostic plots
##################################
##################################

##################################
# Read bio label assignments
##################################
classification_dir = Path("classification")
bio_label_dir = classification_dir
bio_labels = {} # {celltype: {int_label: bio_label}}

for celltype in bio_label_dir.listdir():
    celltype = celltype.basename()
    mnem_fname = bio_label_dir / celltype / "mnemonics.txt"
    if mnem_fname.exists():
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
    length_dist_fname = length_dist_dir / celltype / "plots" / "segment_sizes.tab"
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
	length_distr_fname = length_dist_dir / key / "plots" / "length_distribution.tab"
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
	length_distr_fname = length_dist_dir / key / "plots" / "length_distribution.tab"
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
##################################
# Functional importantance ~ bio label
##################################
log_mem()

label_features_fname = experiment_dir / "../35_2016-08-12_stats_genomedata/nobackup/03_2016-10-11_featuresonly_mean-identity_75q_all/label_feature_vals.pickle"
# {feature_name: {celltype: {label: feature_val} } }
label_feature_vals = pickle.load(open(label_features_fname, "r"))
#feature_name = "mean"

#if not path("stats/functionality").exists():
    #path("stats/functionality").makedirs()

feature_name = "75q_all"
key = "stats/functionality".format(**locals())
with open("{key}.tab".format(**locals()), "w") as f:
    f.write("celltype\tint_label\tbio_label\tfunctionality\n")
    for celltype in bio_labels:
        if celltype in label_feature_vals:
            for int_label in bio_labels[celltype]:
                bio_label = bio_labels[celltype][int_label]
                functionality = label_feature_vals[celltype][int_label]
                f.write("{celltype}\t{int_label}\t{bio_label}\t{functionality}\n".format(**locals()))

data_fname = workdir / "{key}.tab".format(**locals())
script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
#require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)
data$bio_label <- ordered(data$bio_label, {bio_label_order_R})

p <- ggplot(data, aes(x=functionality, color=bio_label, y=1-..y..)) +
  stat_ecdf(size=1.3) +
  scale_x_continuous() +
  theme_bw()
ggsave("{key}_ecdf.pdf", p, width=7, height=5, units="in")

p <- ggplot(data, aes(x=functionality)) +
  stat_bin() +
  facet_grid(bio_label ~ ., scales="free") +
  scale_x_continuous(name="Functionality score", limits=c(0.5, 1.1)) +
  scale_y_continuous(name="Density", breaks=NULL) +
  theme_classic() +
  theme(strip.text.y = element_text(angle=0))
ggsave("{key}_histogram.pdf", p, width=4, height=4, units="in")
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)

##################################
# Overlap with reference annotations
##################################


log_mem()
total_counts = 0
ref_total_counts = {ref_bio_label: 0 for ref_bio_label in all_ref_bio_labels}
ann_total_counts = {ann_bio_label: 0 for ann_bio_label in all_bio_labels}
reference_overlap = {ref_bio_label: {ann_bio_label: 0 for ann_bio_label in all_bio_labels} for ref_bio_label in all_ref_bio_labels}
reference_overlap_raw = {}

reference_overlap_dir = project_dir / "results/29_2016-06-03_reference_ann_overlap/nobackup/01_2016-06-03_run"
for ann_index, ann in enumerate(reference_anns):
    celltype = ann["celltype"]
    concat_key = ann["concatenation_key"]
    ann_id = "{ann_index}_{celltype}_{concat_key}".format(**locals())
    reference_overlap_raw[ann_id] = {}
    reference_overlap_path = reference_overlap_dir / ann["celltype"] / ann["concatenation_key"] / "overlap.tab"
    if not reference_overlap_path.exists():
        logger.warning("Missing {reference_overlap_path}".format(**locals()))
        continue
    if not ann["celltype"] in bio_labels:
        celltype = ann["celltype"]
        logger.warning("bio label assignments is missing celltype: {celltype}".format(**locals()))
        continue
    with open(reference_overlap_path, "r") as f:
        ref_orig_labels = []
        for line in f:
            # header line
            if line.startswith("# mode=bases"):
                continue
            # ref orig labels
            elif line.startswith("\t"):
                ref_orig_labels = line.split()[:-2] # remove "none", "total"
                for ref_orig_label in ref_orig_labels:
                    ref_bio_label = label_mappings[ann["concatenation_key"]][ref_orig_label]
                    ref_label_id = "{ref_bio_label}_{ref_orig_label}".format(**locals())
                    reference_overlap_raw[ann_id][ref_label_id] = {}
            else:
                line = line.split()
                ann_int_label = line[0]
                ann_bio_label = bio_labels[ann["celltype"]][ann_int_label]
                ann_label_id = "{ann_bio_label}_{ann_int_label}".format(**locals())
                counts = map(int, line[1:-2]) # remove label, "none", "total"
                assert (len(counts) == len(ref_orig_labels))
                for ref_orig_label_index, ref_orig_label in enumerate(ref_orig_labels):
                    #reference_overlap[ann_index][ref_orig_label][ann_int_label] = counts[ref_orig_label_index]
                    count = counts[ref_orig_label_index]
                    ref_bio_label = label_mappings[ann["concatenation_key"]][ref_orig_label]
                    ref_label_id = "{ref_bio_label}_{ref_orig_label}".format(**locals())
                    reference_overlap_raw[ann_id][ref_label_id][ann_label_id] = count
                    if ref_bio_label in all_ref_bio_labels:
                        total_counts += count
                        ann_total_counts[ann_bio_label] += count
                        ref_total_counts[ref_bio_label] += count
                        reference_overlap[ref_bio_label][ann_bio_label] += count


key = "stats/reference_overlap"
data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("ref_bio_label\tann_bio_label\tcount\texpected\tenr\n")
    for ref_bio_label in all_ref_bio_labels:
        for ann_bio_label in all_bio_labels:
            count = reference_overlap[ref_bio_label][ann_bio_label]
            expected = (float(total_counts)
                        * (float(ref_total_counts[ref_bio_label]) / total_counts)
                        * (float(ann_total_counts[ann_bio_label]) / total_counts))
            if expected > 0:
                enr = numpy.log2(float(count) / expected)
            else:
                print("expected = 0; ", count, ref_total_counts[ref_bio_label], ann_total_counts[ann_bio_label], ref_bio_label, ann_bio_label)
                enr = 0.0
            f.write("{ref_bio_label}\t{ann_bio_label}\t{count}\t{expected}\t{enr}\n".format(**locals()))

script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)
data$enr_disc = cut(data$enr, breaks=c(-Inf, 0, 1.5, 4.0, Inf))
palette="Reds"
colors <- brewer.pal(palette, n=length(levels(data$enr_disc)))
color_mapping <- levels(data$enr_disc)

data$ann_bio_label <- ordered(data$ann_bio_label, {bio_label_order_R})
data$ref_bio_label <- ordered(data$ref_bio_label, {bio_label_order_ref_R})
data$enr_text <- sprintf("%.1f", data$enr)

  #geom_text(size=4) +
p <- ggplot(data, aes(x=ref_bio_label, y=ann_bio_label, fill=enr_disc, label=enr_text)) +
  geom_tile() +
  scale_fill_manual(values=colors, breaks=color_mapping, drop=FALSE, name="Overlap\nenrichment\nlog2(obs/expected)") +
  scale_x_discrete(name="Reference label") +
  scale_y_discrete(name="New label") +
  theme_bw() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("{key}.pdf", p, width=5, height=3.5, units="in")
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)


if not Path("indiv").exists():
    Path("indiv").makedirs()
key = "indiv/reference_overlap"
data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("ann_id\tref_label_id\tann_label_id\tcount\texpected\tenr\tsame_bio_label\n")
    for ann_id in reference_overlap_raw:
        for ref_label_id in reference_overlap_raw[ann_id]:
            for ann_label_id in reference_overlap_raw[ann_id][ref_label_id]:
                count = reference_overlap_raw[ann_id][ref_label_id][ann_label_id]
                ref_total_counts = sum([reference_overlap_raw[ann_id][ref_label_id][ann_label_id2] for ann_label_id2 in reference_overlap_raw[ann_id][ref_label_id]])
                ann_total_counts = sum([reference_overlap_raw[ann_id][ref_label_id2][ann_label_id] for ref_label_id2 in reference_overlap_raw[ann_id]])
                total_counts = sum([sum([reference_overlap_raw[ann_id][ref_label_id2][ann_label_id2]
                                         for ann_label_id2 in reference_overlap_raw[ann_id][ref_label_id2]])
                                    for ref_label_id2 in reference_overlap_raw[ann_id]])
                expected = (float(total_counts)
                            * (float(ref_total_counts) / total_counts)
                            * (float(ann_total_counts) / total_counts))
                ann_bio_label = ann_label_id.split("_")[0]
                ref_bio_label = ref_label_id.split("_")[0]
                same_bio_label =  ann_bio_label == ref_bio_label
                if expected > 0:
                    enr = numpy.log2(float(count) / expected)
                else:
                    print("expected = 0; ", count, ref_total_counts[ref_bio_label], ann_total_counts[ann_bio_label], ref_bio_label, ann_bio_label)
                    enr = 0.0
                f.write("{ann_id}\t{ref_label_id}\t{ann_label_id}\t{count}\t{expected}\t{enr}\t{same_bio_label}\n".format(**locals()))


script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)
data$enr_disc = cut(data$enr, breaks=c(-Inf, -1.5, -0.5, 0.5, 2.0, Inf))
palette="RdBu"
colors <- brewer.pal(palette, n=length(levels(data$enr_disc)))
color_mapping <- levels(data$enr_disc)
colors <- rev(colors)

p <- ggplot(data) +
  geom_tile(data=data[data$same_bio_label == "False", ], mapping=aes(x=ref_label_id, y=ann_label_id, fill=enr_disc)) +
  geom_tile(data=data[data$same_bio_label == "True", ], mapping=aes(x=ref_label_id, y=ann_label_id, fill=enr_disc), color="black") +
  facet_wrap(~ ann_id, scales="free") +
  scale_fill_manual(values=colors, breaks=color_mapping, drop=FALSE) +
  scale_color_manual(values=c("True"="white", "False"="black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("{key}.pdf", p, width=25, height=25, units="in")
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)

##################################
# track vs bio_label signal distribution heatmap
##################################

log_mem()
assay_name_translation_path = experiment_dir / "assay_name_translation.txt"
assay_name_translation = {}
with open(assay_name_translation_path, "r") as f:
    for line in f:
        line = line.split()
        assay_name_translation[line[0]] = line[1]


signal_dist_experiment_dir = experiment_dir / "../34_2016-07-27_combined_signal_distribution/nobackup/02_2016-10-01_all_pairs/"
target_pairs_path = signal_dist_experiment_dir / "target_pairs.txt"
target_pairs = []
assay_target_celltypes = {}
with open(target_pairs_path, "r") as f:
    for line in f:
        line = line.split()
        celltype = line[0].strip("/")
        assaytype = line[1].strip("/")
        target_pairs.append({"celltype": celltype, "assaytype": assaytype})
        if not assaytype in assay_target_celltypes:
            assay_target_celltypes[assaytype] = []
        assay_target_celltypes[assaytype].append(celltype)

bio_label_signal_dist = {assaytype:
                         {bio_label:
                          {"mean": 0, "n": 0}
                          for bio_label in all_bio_labels}
                         for assaytype in assay_target_celltypes.keys()
                         if len(assay_target_celltypes.keys()) >= 5}
for assaytype in bio_label_signal_dist.keys():
    for celltype in assay_target_celltypes[assaytype]:
        if not (celltype in bio_labels):
            logger.warning("Skipping {assaytype} signal_dist for {celltype} -- missing bio_label assignments".format(**locals()))
            continue
        signal_dist_path = signal_dist_experiment_dir / "{assaytype}/{celltype}/signal_dist/signal_distribution.tab".format(**locals())
        if signal_dist_path.exists():
            with open(signal_dist_path, "r") as f:
                for line in f:
                    if line.startswith("label"):
                        continue # header line
                    line = line.split()
                    int_label = line[0]
                    trackname = line[1]
                    mean = float(line[2])
                    sd = float(line[3])
                    n = float(line[4])
                    bio_label = bio_labels[celltype][str(int_label)]
                    old_n = bio_label_signal_dist[assaytype][bio_label]["n"]
                    old_mean = bio_label_signal_dist[assaytype][bio_label]["mean"]
                    new_n = old_n + n
                    new_mean = ((mean * n) + (old_mean * old_n)) / (old_n + n)
                    bio_label_signal_dist[assaytype][bio_label]["n"] = new_n
                    bio_label_signal_dist[assaytype][bio_label]["mean"] = new_mean


selected_tracks = [ 'FAIRESEQ', 'REPLISEQ', 'HISTONE.H2A.Z', 'HISTONE.H3K4ME1', 'DNASESEQ', 'TFBS.EP300', 'TFBS.CTCF', 'HISTONE.H3K36ME3', 'HISTONE.H3K27ME3', 'HISTONE.H3K18AC', 'HISTONE.H3K27AC', 'HISTONE.H3K9AC', 'TFBS.REST', 'TFBS.POLR2A', 'HISTONE.H3K9ME3', 'HISTONE.H4K20ME1', 'HISTONE.H3K79ME2', 'HISTONE.H3K79ME1', 'HISTONE.H3K4ME2', 'HISTONE.H3K4ME3']

key = "stats/combined_signal_dist"
data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("assaytype\tbio_label\tmean\tn\n")
    for assaytype in bio_label_signal_dist:
        if not (assaytype in selected_tracks): continue
        assaytype_str = assay_name_translation[assaytype]
        assaytype_means = [bio_label_signal_dist[assaytype][bio_label]["mean"] for bio_label in bio_label_signal_dist[assaytype] if bio_label_signal_dist[assaytype][bio_label]["n"] > 0]
        if len(assaytype_means) == 0:
            continue
        assaytype_max = max(assaytype_means)
        assaytype_min = min(assaytype_means)
        #if numpy.abs(assaytype_max - assaytype_min) < 1e-4:
            #print assaytype, assaytype_max, assaytype_min, assaytype_means
            #logger.warning("min = max for {assaytype} -- skipping".format(**locals()))
            #continue
        for bio_label in all_bio_labels:
            mean = bio_label_signal_dist[assaytype][bio_label]["mean"]
            mean_scale = (mean - assaytype_min) / (assaytype_max - assaytype_min)
            n = bio_label_signal_dist[assaytype][bio_label]["n"]
            if n > 0:
                f.write("{assaytype_str}\t{bio_label}\t{mean_scale}\t{n}\n".format(**locals()))

script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
require(reshape2)
require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)

data$bio_label <- ordered(data$bio_label, {bio_label_order_R})

data$mean_disc <- cut(data$mean, breaks=c(0, 0.3, 0.5, 0.7, 1.0), include.lowest=TRUE)
palette="Reds"
colors <- brewer.pal(palette, n=length(levels(data$mean_disc)))
color_mapping <- levels(data$mean_disc)

mat_data <- t(acast(data, bio_label ~ assaytype, value.var="mean"))
ord <- hclust( dist(mat_data, method = "euclidean") )$order
data$assaytype_ord <- factor(data$assaytype, levels=rownames(mat_data)[ord])

p <- ggplot(data, aes(x=assaytype_ord, y=bio_label, fill=mean_disc)) +
  geom_tile() +
  theme_classic() +
  scale_fill_manual(values=colors, name="Mean signal\nvalue", breaks=color_mapping, drop=FALSE) +
  scale_x_discrete(name="Assay type") +
  scale_y_discrete(name="Label") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("{key}.pdf", p, width=6, height=3, units="in")
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)










##################################
# segtools-aggregation
##################################
log_mem()
if False:
    segtools_aggregation_dir = Path("stats/aggregation")
    for celltype in common_celltypes:
        target_dir = segtools_aggregation_dir / celltype
        if not target_dir.exists():
            target_dir.makedirs()
        aggregation_data_path = experiment_dir / "../28_2016-06-02_gencode_aggregation" / celltype / "feature_aggregation.tab"
        aggregation_data_path.copy(target_dir / "feature_aggregation.tab")
        mnem_fname = bio_label_dir / celltype / "mnemonics.txt"
        cmd = ["segtools-aggregation", "_", "_",
               "-m", mnem_fname,
               "--mode=gene",
               "--normalize",
               "--outdir={target_dir}".format(**locals()),
               "--replot"
              ]
        logger.info(" ".join(cmd))
        subprocess.check_call(cmd)



##################################
# Combined gene aggregation
##################################

log_mem()
segtools_features_path = project_dir / "results/39_2016-09-10_ann_features"
component_order = []
combined_gene_aggregation_num_bp = {} # {label: count}
combined_gene_aggregation_data = {} # {component: {offset: {label: count}}}
aggfile_num_features = None
aggfile_spacers = None
key = "stats/gene_aggregation"
combined_gene_agg_data_fname = Path("{key}.tab".format(**locals()))
if not combined_gene_agg_data_fname.exists():
    for celltype in list(common_celltypes): # XXX
        logger.info("Starting gene aggregation data for celltype: {celltype}...".format(**locals()))
        agg_stats_path = segtools_features_path / celltype / "aggregation/GENCODE/feature_aggregation.tab"
        aggregation_label_order = []
        with open(agg_stats_path, "r") as f:
            ann_num_bp = {} # {label: bp}
            for line in f:
                if line[0] == "#":
                    line = line.split()
                    line = line[1:] # "#" character
                    for entry_index, entry in enumerate(line):
                        entry = entry.split("=")
                        if entry[0] == "num_features":
                            aggfile_num_features = entry[1]
                        elif entry[0] == "spacers":
                            aggfile_spacers = entry[1]
                        else:
                            int_label = entry[0]
                            num_bp = int(entry[1])
                            ann_num_bp[int_label] = num_bp
                            bio_label = bio_labels[celltype][int_label]
                            if not bio_label in combined_gene_aggregation_num_bp:
                                combined_gene_aggregation_num_bp[bio_label] = 0
                            combined_gene_aggregation_num_bp[bio_label] += num_bp
                else:
                    line = line.split("\t")
                    if line[0] == "group":
                        aggregation_label_order = map(lambda x: x.strip(), line[3:])
                    else:
                        group = line[0]
                        component = line[1]
                        offset = line[2]
                        counts = line[3:]
                        if not (component in combined_gene_aggregation_data):
                            component_order.append(component)
                            combined_gene_aggregation_data[component] = {}
                        if not (offset in combined_gene_aggregation_data[component]):
                            combined_gene_aggregation_data[component][offset] = {}
                        for i, count in enumerate(counts):
                            count = int(count)
                            int_label = aggregation_label_order[i]
                            bio_label = bio_labels[celltype][int_label]
                            #count = float(count) / ann_num_bp[int_label] # XXX
                            if not bio_label in combined_gene_aggregation_data[component][offset]:
                                combined_gene_aggregation_data[component][offset][bio_label] = 0
                            combined_gene_aggregation_data[component][offset][bio_label] += count
    #for label in combined_gene_aggregation_num_bp:
        #combined_gene_aggregation_num_bp[label] = len(common_celltypes) # XXX
    with open(combined_gene_agg_data_fname, "w") as f:
        # metadata line
        bio_label_order = list(combined_gene_aggregation_num_bp.keys())
        f.write("# ")
        f.write("num_features={aggfile_num_features} ".format(**locals()))
        f.write("spacers={aggfile_spacers} ".format(**locals()))
        for bio_label in bio_label_order:
            count = combined_gene_aggregation_num_bp[bio_label]
            f.write("{bio_label}={count} ".format(**locals()))
        f.write("\n")
        # header line
        f.write("group\tcomponent\toffset\t")
        f.write("\t".join([bio_label for bio_label in bio_label_order]))
        f.write("\n")
        # data lines
        for component in component_order:
            offset_order = sorted(map(int, combined_gene_aggregation_data[component].keys()))
            for offset in offset_order:
                offset = str(offset)
                f.write("genes\t{component}\t{offset}\t".format(**locals()))
                f.write("\t".join([str(combined_gene_aggregation_data[component][offset][bio_label]) for bio_label in bio_label_order]))
                f.write("\n")

script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
segtools.r.dirname <-
  system2("python",
          c("-c", "'import segtools; print segtools.get_r_dirname()'"),
          stdout = TRUE)

source(file.path(segtools.r.dirname, 'common.R'))
source(file.path(segtools.r.dirname, 'aggregation.R'))

save.gene.aggregations('.',
                       '{key}_feature_aggregation.splicing',
                       '{key}_feature_aggregation.translation',
                       '{key}.tab',
                       normalize = TRUE,
                       mnemonic_file = '',
                       clobber = TRUE,
                       significance = FALSE)
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)

##################################
# Indiv signal distribution
##################################

signal_dist_data = {celltype: {} for celltype in celltypes}
for celltype in celltypes:
    logger.info("Reading signal dist for celltype {}".format(celltype))
    combined_signal_dist_dir = project_dir / "results/34_2016-07-27_combined_signal_distribution/nobackup/02_2016-10-01_all_pairs"
    for assaytype_dir in combined_signal_dist_dir.listdir():
        assaytype = assaytype_dir.basename().strip()
        if (assaytype_dir).isdir():
            signal_dist_file = assaytype_dir / celltype / "signal_dist" / "signal_distribution.tab"
            if signal_dist_file.exists():
                signal_dist_data[celltype][assaytype] = pandas.read_csv(signal_dist_file, sep="\t")

indiv_signal_dist_dir = Path("indiv_signal_dist")
for celltype in celltypes:
    logger.info("Starting plotting signal dist for celltype {}".format(celltype))
    target_dir = indiv_signal_dist_dir / celltype
    if not target_dir.exists():
        target_dir.makedirs()
    key = "indiv_signal_dist/{celltype}/signal_distribution".format(**locals())
    data_fname = workdir / "{key}.tab".format(**locals())
    with open(data_fname, "w") as f:
        f.write("assaytype\tlabel_str\tmean\tmean_norm\n")
        for assaytype in signal_dist_data[celltype]:
            assaytype_max = signal_dist_data[celltype][assaytype]["mean"].max()
            assaytype_min = signal_dist_data[celltype][assaytype]["mean"].min()
            for i in range(signal_dist_data[celltype][assaytype].shape[0]):
                int_label = signal_dist_data[celltype][assaytype]["label"][i]
                if celltype in bio_labels:
                    bio_label = bio_labels[celltype][str(int_label)]
                else:
                    raise Exception()
                    #bio_label = "FeaturesMissing"
                label_str = "{bio_label}_{int_label}".format(**locals())
                mean = signal_dist_data[celltype][assaytype]["mean"][i]
                assaytype_name = assay_name_translation[assaytype]
                if (assaytype_max - assaytype_min) >= 1e-4:
                    mean_norm = (mean - assaytype_min) / (assaytype_max - assaytype_min)
                else:
                    mean_norm = 0.0
                f.write("{assaytype_name}\t{label_str}\t{mean}\t{mean_norm}\n".format(**locals()))
    script_fname = workdir / "{key}.R".format(**locals())
    script = \
"""
require(Cairo)
require(ggplot2)
require(reshape2)
require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)
data$mean_text <- sprintf("%.1f", data$mean)
#data$mean_norm <- (data$mean - min(data$mean)) / (max(data$mean) - min(data$mean))
data$mean_disc <- cut(data$mean_norm, breaks=c(0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1), include.lowest=TRUE)
palette="Reds"
colors <- brewer.pal(palette, n=length(levels(data$mean_disc)))
color_mapping <- levels(data$mean_disc)

mat_data <- t(acast(data, label_str ~ assaytype, value.var="mean_norm"))
ord <- hclust( dist(mat_data, method = "euclidean") )$order
data$assaytype_ord <- factor(data$assaytype, levels=rownames(mat_data)[ord])

mat_data <- t(acast(data, assaytype ~ label_str, value.var="mean_norm"))
ord <- hclust( dist(mat_data, method = "euclidean") )$order
data$label_str_ord <- factor(data$label_str, levels=rownames(mat_data)[ord])

plot_width <- length(levels(data$label_str))*0.4 + 3
plot_height <- length(levels(data$assaytype))*0.3 + 1.5

  #geom_text(size=4) +
p <- ggplot(data) +
  aes(x=label_str_ord, y=assaytype_ord, label=mean_text, fill=mean_disc) +
  geom_tile() +
  scale_fill_manual(values=colors, name="Mean\nasinh(signal)", breaks=color_mapping, drop=FALSE) +
  scale_x_discrete(name="Assay type") +
  scale_y_discrete(name="Annotation label") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("{key}.pdf", p, width=plot_width, height=plot_height, units="in", limitsize=FALSE)
""".format(**locals())
    with open(script_fname, "w") as f:
        f.write(script)
    cmd = ["Rscript", script_fname]
    logger.info(" ".join(cmd))
    subprocess.check_call(cmd)

##################################
# Combined name mnemonic files
##################################
for celltype in celltypes:
    label_index_mnem_file = Path("classification") / celltype / "state_mnemonics.txt"
    with open(label_index_mnem_file, "w") as f:
        f.write("old\tnew\n")
        for int_label, bio_label in bio_labels[celltype].items():
            label_str = "{bio_label}_{int_label}".format(**locals())
            f.write("{int_label}\t{label_str}\n".format(**locals()))

##################################
# Indiv length dist
##################################

indiv_signal_dist_dir = Path("indiv_length_dist")
if not indiv_signal_dist_dir.exists():
    indiv_signal_dist_dir.makedirs()

for celltype in celltypes:
    logger.info("Starting length dist for {}".format(celltype))
    length_dist_dir = project_dir / "results/10_2015-12-05_annotation_stats/nobackup/" / celltype / "length_distribution/"
    target_dir = indiv_signal_dist_dir / celltype
    if target_dir.exists():
        logger.info("{target_dir} exists -- skipping".format(**locals()))
    else:
        shutil.copytree(length_dist_dir, target_dir)
        mnem_path = Path("classification") / celltype / "state_mnemonics.txt"
        cmd = ["segtools-length-distribution",
               "--replot",
               "--clobber",
               "--mnemonic-file={}".format(mnem_path),
               "--outdir={}".format(target_dir),
               "_"
               ]
        logger.info(" ".join(cmd))
        subprocess.check_call(cmd)


##################################
# Relabel anns
##################################
exit() # XXX

log_mem()
relabeled_dir = Path("relabeled")
if not relabeled_dir.exists():
    relabeled_dir.makedirs()

annotations_path = project_dir /  "results/2015-08-13_annotations"
for celltype in bio_labels:
    logger.info(celltype)
    mnem_path = Path("classification") / celltype / "mnemonics.txt"
    ann_path = annotations_path / celltype / "identify/segway.bed.gz"
    out_path = relabeled_dir / "{celltype}.bed.gz".format(**locals())
    cmd = "segtools-relabel {ann_path} {mnem_path} | gzip > {out_path}".format(**locals())
    logger.info(cmd)
    subprocess.check_call(cmd, shell=True)



##################################
# Combined gene aggregation -- my own plotting script
##################################

test_gene_agg_celltypes = ["SIGMOID_COLON", "A549"]

# Read gene_aggregation_data
annotations_path = project_dir /  "results/2015-08-13_annotations"
gene_aggregation_num_bp = {} # {celltype: {label: count}}
gene_aggregation_data = {} # {celltype: {component: {offset: {label: count}}}}
aggfile_num_features = None
aggfile_spacers = None
if not Path("indiv").exists():
    Path("indiv").makedirs()
#for celltype in list(common_celltypes):
for celltype in test_gene_agg_celltypes: # XXX
    gene_aggregation_num_bp[celltype] = {}
    gene_aggregation_data[celltype] = {}
    logger.info("Starting gene aggregation data for celltype: {celltype}...".format(**locals()))
    agg_stats_path = annotations_path / celltype / "aggregation/GENCODE/feature_aggregation.tab"
    aggregation_label_order = []
    ann_num_bp = {}
    with open(agg_stats_path, "r") as f:
        for line in f:
            if line[0] == "#":
                line = line.split()
                line = line[1:] # "#" character
                for entry_index, entry in enumerate(line):
                    entry = entry.split("=")
                    if entry[0] == "num_features":
                        aggfile_num_features = entry[1]
                    elif entry[0] == "spacers":
                        aggfile_spacers = entry[1]
                    else:
                        int_label = entry[0]
                        num_bp = int(entry[1])
                        ann_num_bp[int_label] = num_bp
                        #bio_label = bio_labels[celltype][int_label]
                        bio_label = int_label # XXX
                        if not bio_label in gene_aggregation_num_bp[celltype]:
                            gene_aggregation_num_bp[celltype][bio_label] = 0
                        gene_aggregation_num_bp[celltype][bio_label] += num_bp
            else:
                line = line.split("\t")
                if line[0] == "group":
                    aggregation_label_order = map(lambda x: x.strip(), line[3:])
                else:
                    group = line[0]
                    component = line[1]
                    offset = line[2]
                    counts = line[3:]
                    if not (component in gene_aggregation_data[celltype]):
                        gene_aggregation_data[celltype][component] = {}
                    if not (offset in gene_aggregation_data[celltype][component]):
                        gene_aggregation_data[celltype][component][offset] = {}
                    for i, count in enumerate(counts):
                        count = int(count)
                        int_label = aggregation_label_order[i]
                        #bio_label = bio_labels[celltype][int_label]
                        bio_label = int_label # XXX
                        #count = float(count) / ann_num_bp[int_label] # XXX
                        if not bio_label in gene_aggregation_data[celltype][component][offset]:
                            gene_aggregation_data[celltype][component][offset][bio_label] = 0
                        gene_aggregation_data[celltype][component][offset][bio_label] += count


not_used_components = [
    "initial 3' UTR (566 bp)",
    "internal 3' UTR (134 bp)",
    "3' UTR introns (3921 bp)",
    "internal 5' UTR (118 bp)",
    "terminal 3' UTR (642 bp)",
    'initial CDS (180 bp)',
    "initial 5' UTR (169 bp)",
    "5' UTR introns (12317 bp)",
    'terminal CDS (210 bp)',
    "terminal 5' UTR (112 bp)",
]
component_order = [
    "5' flanking: 10000 bp",
    'initial exon (265 bp)',
    'initial intron (11049 bp)',
    'internal exons (144 bp)',
    'internal introns (5308 bp)',
    'terminal exon (645 bp)',
    'terminal intron (5744 bp)',
    "3' flanking: 10000 bp",
]
num_components = len(component_order)

# calculate position_indices to combine everything into one series
arbitrary_celltype = gene_aggregation_data.keys()[0]
offset_order = {component:
                map(str, sorted(map(int, gene_aggregation_data[arbitrary_celltype][component].keys())))
                for component in gene_aggregation_data[arbitrary_celltype]}
#position_indices = {} # {component: {offset: position_index}}
component_offset_order = []
next_position_index = 0
for component in component_order:
    position_indices[component] = {}
    component_offsets = offset_order[component]
    num_offsets = len(component_offsets)
    assert ((num_offsets % 25) == 0)
    for i,j in enumerate(range(0, num_offsets, num_offsets/25)):
        offset = component_offsets[j]
        component_offset_order.append({"component": component, "offset": offset})
        next_position_index += 1


# Plot raw counts for test cell types
key = "test_gene_agg/raw"
data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("celltype\tint_label\tposition_index\tcomponent\toffset\tcounts\n")
    for celltype in test_gene_agg_celltypes:
        for int_label in bio_labels[celltype]:
            for position_index,x in enumerate(component_offset_order):
                component = x["component"]
                offset = x["offset"]
                counts = gene_aggregation_data[celltype][component][offset][int_label]
                f.write("{celltype}\t{int_label}\t{position_index}\t{component}\t{offset}\t{counts}\n".format(**locals()))

script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
#require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)

for (celltype in levels(data$celltype)) {{
    celltype_data <- data[data$celltype==celltype,]
    p <- ggplot(celltype_data, aes(x=position_index, y=counts)) +
      geom_line(stat="identity") +
      geom_hline(yintercept=0) +
      facet_wrap(~ int_label) +
      theme_bw()
    for (x in seq(1,{num_components})) {{
      p <- p + geom_vline(xintercept=25*x, linetype="dashed")
    }}
  ggsave(paste("{key}_", celltype, ".pdf", sep=""), p, width=12, height=10, units="in", limitsize=FALSE)
}}
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)


# calculate enrichments
# {celltype: {component: [position_index: {int_label: enrichment}]}}
gene_aggregation_enrichment = {}
for celltype in gene_aggregation_data:
    gene_aggregation_enrichment[celltype] = []
    for position_index,x in enumerate(component_offset_order):
        component = x["component"]
        offset = x["offset"]
        gene_aggregation_enrichment[celltype].append({})
        for int_label in bio_labels[celltype].keys():
            obs_count = gene_aggregation_data[celltype][component][offset][int_label]
            label_freq = float(gene_aggregation_num_bp[celltype][int_label]) / sum(gene_aggregation_num_bp[celltype].values())
            pos_counts = sum(gene_aggregation_data[celltype][component][offset].values())
            expected_count = pos_counts * label_freq
            val = numpy.log2(float(obs_count) / expected_count)
            gene_aggregation_enrichment[celltype][position_index][int_label] = val


key = "indiv/gene_aggregation"
data_fname = workdir / "{key}.tab".format(**locals())
with open(data_fname, "w") as f:
    f.write("celltype\tposition_index\tbio_label\tlabel\tenrichment\n")
    for celltype in gene_aggregation_enrichment:
        for position_index,x in enumerate(component_offset_order):
            component = x["component"]
            offset = x["offset"]
            for int_label in bio_labels[celltype].keys():
                bio_label = bio_labels[celltype][int_label]
                label_key = "{int_label}_{bio_label}_{celltype}".format(**locals())
                #enrichment = numpy.mean([gene_aggregation_enrichment[celltype][component][offset][int_label] for celltype in gene_aggregation_data])
                enrichment = gene_aggregation_enrichment[celltype][position_index][int_label]
                #position_index = position_indices[component][offset]
                line = "{celltype}\t{position_index}\t{bio_label}\t{label_key}\t{enrichment}\n".format(**locals())
                f.write(line)


script_fname = workdir / "{key}.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
#require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)

for (bio_label in levels(data$bio_label)) {{
    bio_label_data <- data[data$bio_label==bio_label,]
    plot_height <- length(unique(bio_label_data$label))*1
    p <- ggplot(bio_label_data, aes(x=position_index, y=enrichment)) +
      geom_line(stat="identity") +
      geom_hline(yintercept=0) +
      facet_grid(label ~ .) +
      theme_bw()
  ggsave(paste("{key}_", bio_label, ".pdf", sep=""), p, width=2.5, height=plot_height, units="in", limitsize=FALSE)
}}
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)


test_gene_agg_dir = path("test_gene_agg")
if not test_gene_agg_dir.exists():
    test_gene_agg_dir.makedirs()

# copy orig gene agg polots
test_gene_agg_orig_dir = test_gene_agg_dir / "orig"
if not test_gene_agg_orig_dir.exists():
    test_gene_agg_orig_dir.makedirs()

for celltype in test_gene_agg_celltypes:
    orig_plot_path = Path("/net/noble/vol2/home/maxwl/joint/results/2015-08-13_annotations/{celltype}/aggregation/GENCODE/feature_aggregation.splicing.png".format(**locals()))
    orig_plot_path.copy(test_gene_agg_orig_dir / "{celltype}.png".format(**locals()))

# indiv celltype plots
celltypes_str = "c(" + ",".join(map(lambda x: '"{x}"'.format(**locals()), test_gene_agg_celltypes)) + ")"
script_fname = workdir / "{key}_celltype.R".format(**locals())
script = \
"""
require(Cairo)
require(ggplot2)
#require(RColorBrewer)

data <- read.delim("{data_fname}", header=TRUE)

for (celltype in {celltypes_str}) {{
    celltype_data <- data[data$celltype==celltype,]
    p <- ggplot(celltype_data, aes(x=position_index, y=enrichment)) +
      geom_line(stat="identity") +
      geom_hline(yintercept=0) +
      facet_wrap(~ label) +
      theme_bw()
  ggsave(paste("{key}_celltype_", celltype, ".pdf", sep=""), p, width=12, height=10, units="in", limitsize=FALSE)
}}
""".format(**locals())

with open(script_fname, "w") as f:
    f.write(script)

cmd = ["Rscript", script_fname]
logger.info(" ".join(cmd))
subprocess.check_call(cmd)


# rerun segtools with no enr calculation
segtools_aggregation_dir = Path("test_gene_agg/segtools_rerun")
for celltype in test_gene_agg_celltypes:
    target_dir = segtools_aggregation_dir / celltype
    if not target_dir.exists():
        target_dir.makedirs()
    aggregation_data_path = experiment_dir / "../28_2016-06-02_gencode_aggregation" / celltype / "feature_aggregation.tab"
    aggregation_data_path.copy(target_dir / "feature_aggregation.tab")
    #mnem_fname = bio_label_dir / celltype / "mnemonics.txt"
           #"-m", mnem_fname,
    cmd = ["segtools-aggregation", "_", "_",
           "--mode=gene",
           "--outdir={target_dir}".format(**locals()),
           "--replot"
          ]
    logger.info(" ".join(cmd))
    subprocess.check_call(cmd)




