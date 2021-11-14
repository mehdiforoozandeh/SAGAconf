
import sys
import os
import argparse
import subprocess
import gzip
import math
import random
from path import Path
#from pathlib import Path
from collections import defaultdict
#import bedtools
import shutil
import os

import numpy



#################################################
# Classes for reading segtools directories
#################################################

# Stores the information in a segtools-signal-distribution (SD)
# file. Data structures are:
# SDAnnotation.signal_distributions: {label: {"mean": mean, "sd": sd, "n": n}}
# SDAnnotation.label_names = []
# SDAnnotation.track_names = []
class SDAnnotation:
    def __init__(self, file_path):
        self.file_path = file_path
        self.signal_distributions = {}
        self.label_names = set()
        self.track_names = set()
        with open(file_path, 'r') as f:
            header = f.readline()
            for line in f:
                line_data = line.split("\t")
                label_name = line_data[0]
                track_name = line_data[1]
                mean = float(line_data[2])
                sd = float(line_data[3])
                n = int(line_data[4])
                self.label_names.add(label_name)
                self.track_names.add(track_name)
                if not label_name in self.signal_distributions:
                    self.signal_distributions[label_name] = {}
                self.signal_distributions[label_name][track_name] = {"mean": mean, "sd": sd, "n": n}
        self.label_names = list(self.label_names)
        self.track_names = list(self.track_names)

def agg_parse_header(header):
    label_name_bases = header.split()
    labels = {}
    label_names = []
    if label_name_bases[0] == "#":
        for label in label_name_bases[1:]:
            label_data = label.split("=")
            label_name = label_data[0]
            if "num_features" == label_name:
                continue
            if "spacers" == label_name:
                continue
            label_bases = int(label_data[1])
            label = AggLabel(label_name, label_bases)
            labels[label_name] = label
            label_names.append(label_name)
    return labels, label_names

def agg_parse(file_path):
    with open(file_path, 'r') as f:
        header = f.readline().rstrip()
        labels, label_names = agg_parse_header(header)
        second_header = f.readline()
        label_names = second_header.rstrip().split("\t")[3:]
        for line in f:
            line_data = line.rstrip().split("\t")
            group = line_data[0]
            component = line_data[1]
            offset = line_data[2]
            label_counts = line_data[3:]
            assert len(label_counts) == len(label_names)
            #label_names = [ int(label_name) for label_name in label_names] # sort to match counts
            #label_names.sort() # sort to match counts
            #label_names = [ str(label_name) for label_name in label_names] # sort to match counts
            for label_name, label_count in zip(label_names, label_counts):
                labels[label_name].add_raw_count(component,int(label_count))
    return labels, label_names

class AggLabel:
    def __init__(self, label_name, label_bases):
        self.label_name = label_name
        self.label_bases = label_bases
        self.component_raw_counts = {}
        self.component_enrichment = {}
    def add_raw_count(self, component, label_count):
        if component not in self.component_raw_counts:
            self.component_raw_counts[component] = []
        self.component_raw_counts[component].append(label_count)
    def num_bases(self):
        return self.label_bases
    def raw_counts(self):
        return self.component_raw_counts
    def set_enrichment(self,component_enrichment):
        self.component_enrichment = component_enrichment
    def enrichment(self,*args):
        if len(args) == 0:
            return sum([self.component_enrichment[component] for component in self.component_enrichment],[])
        return self.component_enrichment[args[0]]
    def raw_enrichment(self,component):
        return self.component_raw_counts[component]
    def name(self):
        return self.label_name
    def component_enrichments(self):
        return self.component_enrichment

# Stores information from a segtools-aggregation file.
# Data structures:
# AggAnnotation.labels: {label: AggLabel}
# AggLabel.component_enrichment: 
class AggAnnotation:
    def __init__(self, file_path):
        self.file_path = file_path
        self.labels, self.label_names = agg_parse(file_path)
        self.set_enrichment()
    def __iter__(self):
        return iter(self.labels.values())
    def set_enrichment(self):
        genome_bases = 0.0
        component_sum_counts = {}
        for label in self.labels:
            genome_bases += self.labels[label].num_bases()
            component_raw_counts = self.labels[label].raw_counts()
            for component in component_raw_counts:
                if component not in component_sum_counts:
                    component_sum_counts[component] = numpy.zeros(len(component_raw_counts[component]))
                component_sum_counts[component] += numpy.array(component_raw_counts[component])
        for label in self.labels:
            component_raw_counts = self.labels[label].raw_counts()
            component_enrichment = {}
            for component in component_raw_counts:
                if component not in component_enrichment:
                    component_enrichment[component] = []
                for raw_count, sum_count in zip(component_raw_counts[component],component_sum_counts[component]):
                    f_obs = ((raw_count + 1)/ sum_count)
                    f_rand = (self.labels[label].num_bases() / genome_bases)
                    enr  = math.log((f_obs/f_rand),2)
                    component_enrichment[component].append(enr)
            self.labels[label].set_enrichment(component_enrichment)
    def enrichment(self,label_name,*args):
        if len(args) == 0:
            return self.labels[label_name].enrichment()
        assert len(args) == 1
        if args[0] == "components":
            return self.labels[label_name].component_enrichments()
        else:
            return self.labels[label_name].enrichment(args[0])
    def raw_enrichment(self,label_name,component):
        return self.labels[label_name].raw_enrichment(component)
    def get_labels(self):
        return self.labels
    def get_label_names(self):
        return self.label_names



############################################
# features_from_segtools_dir
############################################

# Map column names from segtools files to our preferred feature names.
def feature_name_map(feature_name):
    if feature_name == "H3K9ME3": return "(01) H3K9me3"
    if feature_name == "H3K27ME3": return "(02) H3K27me3"
    if feature_name == "H3K36ME3": return "(03) H3K36me3"
    if feature_name == "H3K4ME3":  return "(04) H3K4me3"
    if feature_name == "H3K27AC": return "(05) H3K27ac"
    if feature_name == "H3K4ME1":  return "(06) H3K4me1"
    if feature_name == "H3K9me3": return "(01) H3K9me3"
    if feature_name == "H3K27me3": return "(02) H3K27me3"
    if feature_name == "H3K36me3": return "(03) H3K36me3"
    if feature_name == "H3K4me3":  return "(04) H3K4me3"
    if feature_name == "H3K27ac": return "(05) H3K27ac"
    if feature_name == "H3K4me1":  return "(06) H3K4me1"
    if feature_name == "5' flanking (1000-10000 bp)":  return "(07) 5' flanking (1000-10000 bp)"
    if feature_name == "5' flanking (1-1000 bp)": return "(08) 5' flanking (1-1000 bp)"
    if feature_name.startswith("initial exon"): return "(09) initial exon"
    if feature_name.startswith("initial intron"): return "(10) initial intron"
    if feature_name.startswith("internal exons"): return "(11) internal exons"
    if feature_name.startswith("internal introns"): return "(12) internal introns"
    if feature_name.startswith("terminal exon"): return "(13) terminal exon"
    if feature_name.startswith("terminal intron"): return "(14) terminal intron"
    if feature_name == "3' flanking (1-1000 bp)": return "(15) 3' flanking (1-1000 bp)"
    if feature_name == "3' flanking (1000-10000 bp)": return "(16) 3' flanking (1000-10000 bp)"
    else:
        raise Exception("Unrecognized feature name {}".format(feature_name))

########################################
# features_from_segtools_dir
# Take as input two segtools files, from gene-aggregation and signal-distribution respectively.
# Read both files and output a vector of features (ann_features). 
# Also output ann_label_bases (number of bases covered by each label) and
# feature_names (list of feature names), which are used in a few places.
########################################
def features_from_segtools_dir(gencode_path, histone_path, segtools_trackname_mapping):
    feature_names = set()
    ann_features = {} # {label: {feature_name: val} }
    ann_label_bases = {}
    #celltype = ann["celltype"]
    #dataset_key = ann["dataset_key"]
    #concatenation_key = ann["concatenation_key"]
    #gencode_path = summary_dirpath / "aggregation/GENCODE/feature_aggregation.tab"
    #gencode_path = summary_dirpath / "feature_aggregation.tab"
    gencode_labels = set()
    gencode = AggAnnotation(gencode_path)
    for label, gencode_label_info in gencode.labels.items():
        gencode_labels.add(label)
        if not (label in ann_features):
            ann_features[label] = {}
        ann_label_bases[label] = gencode_label_info.num_bases()
        for component_name, component_enrichments in gencode_label_info.component_enrichments().items():
            if ("UTR" in component_name) or ("CDS" in component_name): continue
            # Split 5' and 3' flanking regions into two parts. 
            if component_name.startswith("5' flanking"):
                feature_name = "5' flanking (1-1000 bp)"
                feature_name = feature_name_map(feature_name)
                ann_features[label][feature_name] = numpy.mean(component_enrichments[-int(0.1*len(component_enrichments)):])
                feature_names.add(feature_name)
                feature_name = "5' flanking (1000-10000 bp)"
                feature_name = feature_name_map(feature_name)
                ann_features[label][feature_name] = numpy.mean(component_enrichments[:-int(0.1*len(component_enrichments))])
                feature_names.add(feature_name)
            elif component_name.startswith("3' flanking"):
                feature_name = "3' flanking (1-1000 bp)"
                feature_name = feature_name_map(feature_name)
                ann_features[label][feature_name] = numpy.mean(component_enrichments[:int(0.1*len(component_enrichments))])
                feature_names.add(feature_name)
                feature_name = "3' flanking (1000-10000 bp)"
                feature_name = feature_name_map(feature_name)
                ann_features[label][feature_name] = numpy.mean(component_enrichments[int(0.1*len(component_enrichments)):])
                feature_names.add(feature_name)
            else:
                feature_name = feature_name_map(component_name)
                ann_features[label][feature_name] = numpy.mean(component_enrichments)
                feature_names.add(feature_name)
            assert numpy.isfinite(ann_features[label][feature_name])
    #for histone in histone_features:
        #histone_path = summary_dirpath / "signal_distribution/HISTONE.{histone}/signal_distribution.tab".format(**locals())
    #histone_path = summary_dirpath / "signal_distribution.tab".format(**locals())
    #if not histone_path.exists():
        #logger.warning("!!! Missing {histone_path}!!".format(**locals()))
        #raise Exception("!!! Missing {histone_path}!!".format(**locals()))
        #for label in ann_features:
            #feature_name = histone
            #feature_name = feature_name_map(feature_name)
            #ann_features[label][feature_name] = float("nan")
            #feature_names.add(feature_name)
    #else:
    histone_signal = SDAnnotation(histone_path)
    
    #assert(len(histone_signal.track_names) == 9)
    # signal_labels = set()
    # #histone_features = [ "H3K9ME3", "H3K27ME3", "H3K36ME3", "H3K4ME3", "H3K27AC", "H3K4ME1" ]
    # for label_name in histone_signal.label_names:
    #     signal_labels.add(label_name)
    #     for i,track_name in enumerate(histone_signal.track_names):
    #         # feature_name = histone_features[i]
    #         feature_name = segtools_trackname_mapping[track_name]
    #         feature_name = feature_name_map(feature_name)
    #         ann_features[label_name][feature_name] = histone_signal.signal_distributions[label_name][track_name]["mean"]
    #         assert numpy.isfinite(ann_features[label_name][feature_name])
    #         feature_names.add(feature_name)
    
    signal_labels = set()
    histone_features = [ "H3K9ME3", "H3K27ME3", "H3K36ME3", "H3K4ME3", "H3K27AC", "H3K4ME1" ]
    for label_name in histone_signal.label_names:  # this walks on the clusterer labels
        #print(label_name)
        signal_labels.add(label_name)
        sixcounter = 0
        for i,track_name in enumerate(histone_signal.track_names): # for each label, we want to retrieve the track info as features
            # feature_name = histone_features[i]
            feature_name = segtools_trackname_mapping[track_name] # get the feature name 
            #print(feature_name)
            if ((feature_name.upper() in histone_features)): # to skip the DNase-seq one
                sixcounter += 1
                feature_name = feature_name_map(feature_name)
                #print(sixcounter, feature_name)
                ann_features[label_name][feature_name] = histone_signal.signal_distributions[label_name][track_name]["mean"]
                assert numpy.isfinite(ann_features[label_name][feature_name])
                feature_names.add(feature_name)
            
        # print("found {} out of 6 necessary histone tracks".format(sixcounter))
        
    if not signal_labels == gencode_labels:
        raise Exception("signal_labels and gencode_labels don't match: signal_labels = {}. gencode_labels = {}".format(signal_labels, gencode_labels))
    return ann_features, ann_label_bases, feature_names

#################################3
# features_frame_to_matrix
# Convert features in pandas to numpy matrix (?).
#################################3
def features_frame_to_matrix(features_frame, feature_names):
    num_examples = features_frame.shape[0]
    mat = numpy.empty(shape=(num_examples, len(feature_names)), dtype="float")
    for feature_index, feature_name in enumerate(feature_names):
        for example_index in range(num_examples):
            feature_val = features_frame.loc[example_index, feature_name]
            # XXX Is normalization necessary? Requires remembering mean/stds from training.
            # norm_feature_val = float(feature_val - feature_means[feature_name]) / feature_stdevs[feature_name]
            norm_feature_val = feature_val
            mat[example_index, feature_index] = norm_feature_val
    return mat
