
import collections
from turtle import color, width
from scipy.ndimage import gaussian_filter1d
from statistics import mean
import functools
import multiprocessing as mp
from matplotlib.cm import get_cmap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from _cluster_matching import *
import plotly.graph_objects as go
import plotly.express as px
from sklearn.preprocessing import PolynomialFeatures, MinMaxScaler
from sklearn.isotonic import IsotonicRegression
from sklearn.pipeline import make_pipeline
from scipy.interpolate import UnivariateSpline, interp1d
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from datetime import datetime
from scipy.stats import norm
import math


class Agreement(object):
    def __init__(self, loci_1, loci_2, savedir, filter_nan=True, assume_uniform_coverage=True):
        if filter_nan:
            loci_1 = loci_1.dropna()
            loci_1 = loci_1.reset_index(drop=True)
            loci_2 = loci_2.dropna()
            loci_2 = loci_2.reset_index(drop=True)

        self.epsilon = 1e-3
        self.loci_1 = loci_1
        self.loci_2 = loci_2
        self.savedir = savedir
        self.num_labels = len(self.loci_1.columns)-3

        self.max_posteri1 = list(self.loci_1.iloc[:, 3:].idxmax(axis=1))
        self.max_posteri2 = list(self.loci_2.iloc[:, 3:].idxmax(axis=1))
        self.assume_uniform_coverage = assume_uniform_coverage
        self.expected_agreement = {}

    def per_label_agreement(self):
        label_dict = {}
        label_coverage = {}

        for c in self.loci_1.columns[3:]:
            label_dict[c] = [0, 0] # [agreement_count, disagreement_count]
            label_coverage[c] = 0
            if self.assume_uniform_coverage == False:
                self.expected_agreement[c] = 0
            else:
                self.expected_agreement[c] = float(1) / float(self.num_labels)

        for i in range(self.loci_1.shape[0]):
            if self.assume_uniform_coverage == False:
                self.expected_agreement[self.max_posteri2[i]] += 1 / self.loci_1.shape[0]

            if self.max_posteri1[i] == self.max_posteri2[i]:
                label_dict[self.max_posteri1[i]][0] += 1 

            else:
                label_dict[self.max_posteri1[i]][1] += 1 

        for k in label_dict.keys():
            label_coverage[k] = float(label_dict[k][0]+label_dict[k][1]) / self.loci_1.shape[0]

            if (label_dict[k][0]+label_dict[k][1]) != 0:
                label_dict[k] = label_dict[k][0] / (label_dict[k][0]+label_dict[k][1])
            else:
                label_dict[k] = 0

        self.label_agreements = label_dict
        self.label_coverage = label_coverage
        return self.label_agreements

    def general_agreement(self):
        # self.per_label_agreement()
        # l = list(self.label_agreements.values())

        self.overall_agreement = 0

        for l in self.label_agreements.keys():
            self.overall_agreement += self.label_agreements[l] * self.label_coverage[l]

        return self.overall_agreement

    def per_label_OE_ratio(self, log_transform=True):
        self.label_OE = {}
        for k in self.label_agreements.keys():
            if log_transform:
                self.label_OE[k] = np.log(
                    self.label_agreements[k] / (self.expected_agreement[k] + self.epsilon)
                    )
            else:
                self.label_OE[k] = self.label_agreements[k] / (self.expected_agreement[k] + self.epsilon)

        return self.label_OE

    def general_OE_ratio(self, log_transform=True):
        if log_transform:
            self.overall_OE = np.log(self.overall_agreement/ ((float(1) / float(self.num_labels)) + self.epsilon))
        else:
            self.overall_OE =  self.overall_agreement/ (float(1) / float(self.num_labels) + self.epsilon)

        return self.overall_OE  

    def per_label_cohens_kappa(self):
        self.label_CK = {}

        for k in self.label_agreements.keys():
            self.label_CK[k] = \
                (self.label_agreements[k] - (self.expected_agreement[k] + self.epsilon)) / (1 - (self.expected_agreement[k] + self.epsilon))
                
        return self.label_CK

    def general_cohens_kappa(self):
        self.overall_CK = \
            (self.overall_agreement - (float(1) / float(self.num_labels) + self.epsilon)) / (1 - (float(1) / float(self.num_labels) + self.epsilon))
        return self.overall_CK
    
    def plot_agreement(self):
        x = []
        height = []
        for k, v in self.label_agreements.items():
            x.append(k.replace('posterior', 'label'))
            height.append(v)
        
        x.append('Overall')
        height.append(self.overall_agreement)
        plt.figure(figsize=(11,9))
        plt.bar(x, height, width=0.4, color='black', alpha=0.5)
        plt.xticks(rotation=45, fontsize=7)
        plt.title('Agreement between two replicates')
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.xlabel('Label')
        plt.ylabel('Agreement')
        plt.tight_layout()
        plt.savefig('{}/agreement.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/agreement.svg'.format(self.savedir), format='svg')
        plt.clf()
        with open("{}/agreement.txt".format(self.savedir), 'w') as pltxt:
            """
            title, xaxis, yaxis, x, y
            """
            pltxt.write(
                "{}\n{}\n{}\n{}\n{}".format(
                    'Agreement between two replicates',
                    "Label", 
                    "Agreement",
                    x, height
                )
            )

    def plot_OE(self):
        # self.per_label_OE_ratio()
        x = []
        height = []
        for k, v in self.label_OE.items():
            x.append(k.replace('posterior', 'label'))
            height.append(v)

        self.general_OE_ratio()
        x.append('Overall')
        height.append(self.overall_OE)
        plt.figure(figsize=(11,9))
        plt.bar(x, height, width=0.4, color='black', alpha=0.5)
        plt.xticks(rotation=45, fontsize=7)
        plt.xlabel('Label')
        plt.title('Observed Agreement/Expected Agreement plot')
        plt.ylabel('log(O/E)')
        plt.tight_layout()
        plt.savefig('{}/oe_agreement.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/oe_agreement.svg'.format(self.savedir), format='svg')
        plt.clf()
        with open("{}/oe_agreement.txt".format(self.savedir), 'w') as pltxt:
            """
            title, xaxis, yaxis, x, y
            """
            pltxt.write(
                "{}\n{}\n{}\n{}\n{}".format(
                    'Observed Agreement/Expected Agreement plot',
                    "Label", 
                    "log(O/E)",
                    x, height
                )
            )

    def plot_CK(self):
        # self.per_label_cohens_kappa()
        x = []
        height = []
        for k, v in self.label_CK.items():
            x.append(k.replace('posterior', 'label'))
            height.append(v)

        self.general_cohens_kappa()
        x.append('Overall')
        height.append(self.overall_CK)
        plt.figure(figsize=(11,9))
        plt.bar(x, height, width=0.4, color='black', alpha=0.5)
        plt.xticks(rotation=45, fontsize=7)
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.xlabel('Label')
        plt.title("Cohen's Kappa")
        plt.ylabel("Cohen's Kappa Score")
        plt.tight_layout()
        plt.savefig('{}/cohenskappa.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/cohenskappa.svg'.format(self.savedir), format='svg')
        plt.clf()
        with open("{}/cohenskappa.txt".format(self.savedir), 'w') as pltxt:
            """
            title, xaxis, yaxis, x, y
            """
            pltxt.write(
                "{}\n{}\n{}\n{}\n{}".format(
                    "Cohen's Kappa ",
                    "Label", 
                    "Cohen's Kappa Score",
                    x, height
                )
            )
        
class posterior_calibration(object):
    def __init__(self, loci_1, loci_2, savedir, filter_nan=True, oe_transform=True):
        if filter_nan:
            loci_1 = loci_1.dropna()
            loci_1 = loci_1.reset_index(drop=True)
            loci_2 = loci_2.dropna()
            loci_2 = loci_2.reset_index(drop=True)

        self.loci_1 = loci_1
        self.loci_2 = loci_2
        self.num_labels = len(self.loci_1.columns)-3

        self.savedir = savedir
        del loci_1
        del loci_2
        
        self.oe_transform = oe_transform

        self.MAPestimate1 = self.loci_1.iloc[:,3:].idxmax(axis=1)
        self.MAPestimate2 = self.loci_2.iloc[:,3:].idxmax(axis=1)

    def perlabel_visualize_calibration(self, bins, label_name, scatter=False):
        if scatter:
            plt.scatter(
                x=np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), 
                y=bins[:,5], c='black', s=10000/len(bins))
            
            polyreg = IsotonicRegression(
                    y_min=float(np.array(bins[:, 5]).min()), y_max=float(np.array(bins[:, 5]).max()), 
                    out_of_bounds="clip")
            
            polyreg.fit(
                np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)), 
                bins[:, 5])
            
            plt.plot(
                np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), 
                polyreg.predict(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])), 
                '--', c='r', linewidth=3)

        else:
            plt.plot([(bins[i,0]+bins[i,1])/2 for i in range(len(bins))], bins[:,5], label=label_name)

        plt.title("Reproduciblity Plot {}".format(label_name))
        xlabel = "posterior in Replicate 1"
        if self.oe_transform:
            ylabel = "log(O/E) of Similarly Labeled Bins in replicate 2"
        else:
            ylabel = "log(Ratio) of Similarly Labeled Bins in replicate 2"
       
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.tight_layout()
        plt.savefig('{}/caliberation_{}.pdf'.format(self.savedir, label_name.replace("|","")), format='pdf')
        plt.savefig('{}/caliberation_{}.svg'.format(self.savedir, label_name.replace("|","")), format='svg')
        plt.clf()

        with open("{}/caliberation_{}.txt".format(self.savedir, label_name.replace("|","")), 'w') as pltxt:
            """
            title, xaxis, yaxis, x, y, polyreg
            """
            xlabel = "Posterior in Replicate 1"
            if self.oe_transform:
                ylabel = "O/E of Similarly Labeled Bins in replicate 2"
            else:
                ylabel = "Ratio of Similarly Labeled Bins in replicate 2"

            pltxt.write(
                "{}\n{}\n{}\n{}\n{}\n{}".format(
                    "Reproduciblity Plot {}".format(label_name),
                    xlabel, 
                    ylabel,
                    list(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])), 
                    list(bins[:,5]),
                    list(polyreg.predict(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])))
                )
            )

    def general_visualize_calibration(self, bins_dict, subplot=True):
        if subplot:
            list_binsdict = list(bins_dict.keys())
            num_labels = self.num_labels
            n_cols = math.floor(math.sqrt(num_labels))
            n_rows = math.ceil(num_labels / n_cols)

            fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
            label_being_plotted = 0
            
            for i in range(n_rows):
                for j in range(n_cols):
                    label_name = list_binsdict[label_being_plotted]
                    bins = bins_dict[label_name]

                    polyreg = IsotonicRegression(
                        y_min=float(np.array(bins[:, 5]).min()), y_max=float(np.array(bins[:, 5]).max()), 
                        out_of_bounds="clip")
                    
                    polyreg.fit(
                        np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)), 
                        bins[:, 5])
                    
                    r2 = r2_score(
                        bins[:, 5], 
                        polyreg.predict(
                            np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1))))
                    
                    axs[i,j].plot(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), 
                        polyreg.predict(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])),
                        c="black")

                    axs[i,j].set_title("{}_r2={:.2f}".format(label_name, float(r2)), fontsize=7)

                    label_being_plotted+=1
        

            xlabel = "Posterior in Replicate 1"
            if self.oe_transform:
                ylabel = "log(O/E) of Similarly Labeled Bins in replicate 2"
            else:
                ylabel = "log(Ratio) of Similarly Labeled Bins in replicate 2"


            plt.tight_layout()
            plt.savefig('{}/clb_{}.pdf'.format(self.savedir, "subplot"), format='pdf')
            plt.savefig('{}/clb_{}.svg'.format(self.savedir, "subplot"), format='svg')
            plt.clf()

        colors = [i for i in get_cmap('tab20').colors]
        ci = 0
        for label_name, bins in bins_dict.items():
            
            polyreg = IsotonicRegression(
                    y_min=float(np.array(bins[:, 5]).min()), y_max=float(np.array(bins[:, 5]).max()), 
                    out_of_bounds="clip")
            
            polyreg.fit(
                np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)), 
                bins[:, 5])
            
            r2 = r2_score(
                    bins[:, 5], 
                    polyreg.predict(
                        np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1))))
            
            plt.plot(
                np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), 
                polyreg.predict(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])), 
                label="{}_r2={:.2f}".format(label_name, float(r2)), c=colors[ci])
            ci+=1

        xlabel = "posterior in Replicate 1"
        if self.oe_transform:
            ylabel = "log(O/E) of Similarly Labeled Bins in replicate 2"
        else:
            ylabel = "log(Ratio) of Similarly Labeled Bins in replicate 2"


        plt.legend(loc='upper center', bbox_to_anchor=(0.45, -0.05),
            fancybox=True, ncol=4, fontsize=5)

        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.tight_layout()
        plt.savefig('{}/caliberation_{}.pdf'.format(self.savedir, "general"), format='pdf')
        plt.savefig('{}/caliberation_{}.svg'.format(self.savedir, "general"), format='svg')
        plt.clf()
        ################################################
        
        with open("{}/caliberation_{}.txt".format(self.savedir, "general"), 'w') as pltxt:
            """
            title, xaxis, yaxis, x, y(all labels)
            """
            xlabel = "posterior in Replicate 1"
            if self.oe_transform:
                ylabel = "log(O/E) of Similarly Labeled Bins in replicate 2"
            else:
                ylabel = "log(Ratio) of Similarly Labeled Bins in replicate 2"


            pltxt.write(
                "{}\n{}\n{}\n{}".format(
                    "Reproduciblity Plot {}".format(label_name),
                    xlabel, 
                    ylabel,
                    list(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])), 
                )
            )
            for label_name, bins in bins_dict.items():
                pltxt.write(label_name+ "\t" + str(list(bins[:, 5])) + '\n')
            
    def perlabel_calibration_function(self, method="isoton_reg", num_bins=500, return_caliberated_matrix=True, scale=True, scale_columnwise=False, strat_size="def"):
        perlabel_function = {}
        bins_dict = {}
        new_matrix = [self.loci_1.iloc[:,:3]]
        for k in range(self.num_labels):
            kth_label = self.loci_1.iloc[:,3:].columns[k]
            bins = []
            if strat_size == "def":
                strat_size = int(len(self.loci_1)/num_bins)

            posterior_vector_1 = self.loci_1.iloc[:,3:][kth_label]

            vector_pair = pd.concat([posterior_vector_1, self.MAPestimate2], axis=1)

            vector_pair.columns=["pv1", "map2"]
            vector_pair = vector_pair.sort_values("pv1").reset_index(drop=True)

            for b in range(0, vector_pair.shape[0], strat_size):
                subset_vector_pair = vector_pair.iloc[b:b+strat_size, :]
                # [bin_start, bin_end, num_values_in_bin, num_agreement, num_mislabeled, ratio_agreement]
                observed = len(subset_vector_pair.loc[subset_vector_pair["map2"] == kth_label])
                
                if self.oe_transform:
                    expected = ((1/self.num_labels) * len(subset_vector_pair))

                    if observed!=0:
                        oe = np.log((observed)/(expected))

                        bins.append([
                            float(subset_vector_pair["pv1"][subset_vector_pair.index[0]]), 
                            float(subset_vector_pair["pv1"][subset_vector_pair.index[-1]]),
                            len(subset_vector_pair),
                            len(subset_vector_pair.loc[subset_vector_pair["map2"] == kth_label]),
                            len(subset_vector_pair) - len(subset_vector_pair.loc[subset_vector_pair["map2"] == kth_label]),
                            oe])
                else:
                    bins.append([
                        float(subset_vector_pair["pv1"][subset_vector_pair.index[0]]), 
                        float(subset_vector_pair["pv1"][subset_vector_pair.index[-1]]),
                        len(subset_vector_pair),
                        len(subset_vector_pair.loc[subset_vector_pair["map2"] == kth_label]),
                        len(subset_vector_pair) - len(subset_vector_pair.loc[subset_vector_pair["map2"] == kth_label]),
                        float(observed)/float(len(subset_vector_pair))]) 

            bins = np.array(bins)
            bins_dict[kth_label] = bins
            self.perlabel_visualize_calibration(bins, kth_label, scatter=True)

            if method=="spline":
                f = UnivariateSpline(
                    x= np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]),
                    y=bins[:, 5], ext=3)

                r2 = r2_score(
                    bins[:, 5], 
                    f(
                        np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1))
                    ))
                
                print("R2 score for {}: {}".format(kth_label, r2))

            else:
                if method=="isoton_reg":
                    polyreg = IsotonicRegression(
                    y_min=float(np.array(bins[:, 5]).min()), y_max=float(np.array(bins[:, 5]).max()), 
                    out_of_bounds="clip")

                if method == "poly_reg":
                    make_pipeline(PolynomialFeatures(3), LinearRegression())
                
                polyreg.fit(
                np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)), #x_train is the mean of each bin
                bins[:, 5]) # y_train is oe/ratio_correctly_labeled at each bin

                r2 = r2_score(
                    bins[:, 5], 
                    polyreg.predict(
                        np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1))))
                
                print("R2 score for {}: {}".format(kth_label, r2))
            
            if return_caliberated_matrix:
                x = np.array(self.loci_1.iloc[:,3:][kth_label])
                if method == "spline":
                    calibrated_array = f(np.array(x))
                else:
                    calibrated_array = polyreg.predict(
                        np.reshape(np.array(x), (-1,1)))

                new_matrix.append(pd.DataFrame(calibrated_array, columns=[kth_label]))
                
            else:
                if method == "spline":
                    perlabel_function[k] = [f, r2]

                else:
                    perlabel_function[k] = [polyreg, r2]

        self.general_visualize_calibration(bins_dict)

        if return_caliberated_matrix:
            new_matrix = pd.concat(new_matrix, axis=1)
            
            if scale:
                scaler = MinMaxScaler(feature_range=(0,1))
                if scale_columnwise:
                    new_matrix.iloc[:, 3:] = scaler.fit_transform(new_matrix.iloc[:, 3:])

                else:
                    remember_shape = new_matrix.iloc[:,3:].shape
                    remember_columns = new_matrix.iloc[:,3:].columns
                    new_matrix.iloc[:, 3:] = pd.DataFrame(
                        np.reshape(scaler.fit_transform(np.reshape(np.array(
                            new_matrix.iloc[:, 3:]), (-1,1))), remember_shape), 
                                        columns=remember_columns)
    
            return new_matrix

        else:
            self.perlabel_function = perlabel_function
            return self.perlabel_function

class Posterior_dist(object):
    def __init__(self, loci, savedir):
        self.loci = loci
        self.savedir= savedir

    def plot_posterior_histogram(self, num_bins=10):
        num_labels = int(self.loci.shape[1])-3
        n_cols = math.floor(math.sqrt(num_labels))
        n_rows = math.ceil(num_labels / n_cols)

        fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
        label_being_plotted = 0
        for i in range(n_rows):
            for j in range(n_cols):
                x = self.loci.iloc[:, 3:].iloc[:, label_being_plotted]
                x = x.loc[(x!=0)].reset_index(drop=True)

                # bars = []
                # for b in range(num_bins):
                #     brange = [float(b)/num_bins, float(b+1)/num_bins]
                #     bars.append(
                #         len(list(
                #             x.loc[
                #                 (x > brange[0])&
                #                 (x <= brange[1])]
                #         ))/len(x)
                #     )
                # # bars = np.log(bars)
                # axs[i,j].bar(
                #     x=[float(b)/num_bins for b in range(num_bins)],
                #     height = bars,
                #     width = float(1)/num_bins
                # )

                sns.histplot(ax=axs[i,j], x=x, bins=10, stat="probability")
                axs[i,j].set_title("{}".format(self.loci.columns[3+label_being_plotted]), fontsize=7)
                
                label_being_plotted+=1

        plt.tight_layout()
        plt.savefig('{}/posterior_dir_histogram.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/posterior_dir_histogram.svg'.format(self.savedir), format='svg')
        plt.clf()

    def plot_posterior_track(self):
        num_labels = int(self.loci.shape[1])-3
        n_cols = math.floor(math.sqrt(num_labels))
        n_rows = math.ceil(num_labels / n_cols)

        fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
        label_being_plotted = 0
        for i in range(n_rows):
            for j in range(n_cols):
                x = self.loci.iloc[:, 3:].iloc[:, label_being_plotted]
                axs[i,j].plot(list(range(len(x))), x, linewidth=0.05)
                axs[i,j].set_title("{}".format(self.loci.columns[3+label_being_plotted]))
                label_being_plotted+=1

        plt.tight_layout()
        plt.savefig('{}/posterior_dir_track.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/posterior_dir_track.svg'.format(self.savedir), format='svg')
        plt.clf()

class correspondence_curve(object):
    '''
    implementation of the method described at
    https://arxiv.org/pdf/1110.4705.pdf
    '''
    def __init__(self, loci_1, loci_2, savedir, filter_nan=True):
        if filter_nan:
            loci_1 = loci_1.dropna()
            loci_1 = loci_1.reset_index(drop=True)
            loci_2 = loci_2.dropna()
            loci_2 = loci_2.reset_index(drop=True)

        self.loci_1 = loci_1
        self.loci_2 = loci_2
        del loci_1, loci_2
        self.savedir = savedir
        self.num_labels = len(self.loci_1.columns)-3

        # self.MAPestimate1 = self.loci_1.iloc[\
        #     :, 3:].idxmax(axis=1)
        # self.MAPestimate2 = self.loci_2.iloc[\
        #     :, 3:].idxmax(axis=1)

        # self.max_poster1 = self.loci_1.iloc[\
        #     :, 3:].max(axis=1)
        # self.max_poster2 = self.loci_2.iloc[\
        #     :, 3:].max(axis=1)

    def rank_general(self):
        self.annot_1 = []
        for i in range(self.loci_1.shape[0]):
            self.annot_1.append([
                self.loci_1['chr'][i], self.loci_1['start'][i], self.loci_1['end'][i], self.MAPestimate1[i], self.max_poster1[i]])

        self.annot_2 = []
        for i in range(self.loci_2.shape[0]):
            self.annot_2.append([
                self.loci_2['chr'][i], self.loci_2['start'][i], self.loci_2['end'][i], self.MAPestimate2[i], self.max_poster2[i]])


        self.annot_1 = pd.DataFrame(self.annot_1, columns=['chr', 'start', 'end', 'MAP', 'posterior'])
        self.annot_2 = pd.DataFrame(self.annot_2, columns=['chr', 'start', 'end', 'MAP', 'posterior'])

        self.is_ranked = False

        self.annot_1 = self.annot_1.sort_values(by=['posterior'], axis=0, ascending=False)
        self.annot_2 = self.annot_2.sort_values(by=['posterior'], axis=0, ascending=False)

        self.loci_1 = self.loci_1.reindex(self.annot_1.index).reset_index()
        self.loci_2 = self.loci_2.reindex(self.annot_2.index).reset_index()
        self.is_ranked = True

    def psi(self, t, per_label=True):
        if per_label:
            perlabel_psi = {}
            post_cols = list(self.loci_1.columns)[3:]
            for k in post_cols:
                k_sorted_loci_1 = self.loci_1.sort_values(by=[k], axis=0, ascending=False).reset_index()
                k_sorted_loci_2 = self.loci_2.sort_values(by=[k], axis=0, ascending=False).reset_index()
                
                sigma = 0
                upper_tn_1 = set(k_sorted_loci_1.loc[:int(t*len(self.loci_1)), 'index'])
                upper_tn_2 = set(k_sorted_loci_2.loc[:int(t*len(self.loci_2)), 'index'])

                for i in upper_tn_1:
                    if i in upper_tn_2:
                        sigma +=1

                perlabel_psi[k] = float(sigma) / len(self.loci_1)
            return perlabel_psi

        else:
            if not self.is_ranked:
                self.rank_general()

            sigma = 0
            upper_tn_1 = list(self.loci_1.loc[:int(t*len(self.loci_1)), 'index'])
            upper_tn_2 = list(self.loci_2.loc[:int(t*len(self.loci_2)), 'index'])

            for i in upper_tn_1:
                if i in upper_tn_2:
                    if self.annot_1['MAP'][i] == self.annot_2['MAP'][i]:
                        sigma += 1
            
            return sigma / len(self.loci_1)
    
    def plot_1vs2(self):
        post_cols = list(self.loci_1.columns)[3:]

        # colors = np.array([i for i in get_cmap('tab20').colors])
        for k in post_cols:
            k_sorted_loci_1 = self.loci_1.sort_values(by=[k], axis=0, ascending=False).reset_index()
            k_sorted_loci_2 = self.loci_2.sort_values(by=[k], axis=0, ascending=False).reset_index()
            bb= []
            kk = []
            for b in range(5):
                subset1 = list(k_sorted_loci_1.loc[:int(((b+1)/5)*len(k_sorted_loci_1)), "index"])
                subset2 = list(k_sorted_loci_2.loc[:int(((b+1)/5)*len(k_sorted_loci_2)), "index"])
                x = []
                y = []
                for i in range(len(subset1)):
                    if subset1[i] in subset2:
                        x.append(float(k_sorted_loci_1.loc[
                            (k_sorted_loci_1["index"] == subset1[i]),k
                        ]))
                        y.append(float(k_sorted_loci_2.loc[
                            (k_sorted_loci_2["index"] ==  subset1[i]),k
                        ]))

                print(len(x), len(subset1))
                plt.scatter(x, y, s=5, alpha=1-float(b/5), c="black", label = 'top {} percent'.format(float((b+1)/5)*100))

            plt.legend()
            plt.show()
                

    def plot_curve(self, num_t_steps=100, plot_labels=True, plot_general=True, merge_plots=False, subplot=True):
        self.psi_over_t = {}
        if plot_labels:
            post_cols = list(self.loci_1.columns)[3:]
            for k in post_cols:
                self.psi_over_t[k] = []

            for t in range(0, num_t_steps + 1):
                pl_psi = self.psi(t=float(t)/num_t_steps, per_label=True)

                for k in pl_psi.keys():
                    self.psi_over_t[k].append(pl_psi[k])
            
            self.per_label_psi_over_t = self.psi_over_t

        if plot_general:
            self.psi_over_t['general'] = []
            for t in range(0, num_t_steps + 1):
                t = float(t)/num_t_steps
                self.psi_over_t['general'].append(self.psi(t=t, per_label=False))
        
        # start plotting
        t_list = []
        for i in range(0, num_t_steps + 1):
            t_list.append(i / num_t_steps)

        if merge_plots:
            plt.plot(t_list, t_list, '--', label='Perfect Reproducibility')
            for p in self.psi_over_t.keys():
                URIs = self.psi_over_t[p]
                plt.plot(t_list, URIs, label=str(p))

            plt.title('Correspondence Curve')
            plt.xlabel('t')
            plt.ylabel("PSI")
            plt.legend()
            plt.tight_layout()
            plt.savefig('{}/cc_merged.pdf'.format(self.savedir), format='pdf')
            plt.savefig('{}/cc_merged.svg'.format(self.savedir), format='svg')
            plt.clf()

            with open("{}/cc_merged.txt".format(self.savedir), 'w') as pltxt:
                """
                title, xaxis, yaxis, x(t_list), y
                """

                pltxt.write(
                    "{}\n{}\n{}\n{}\n{}".format(
                        "Correspondence Curve",
                        "t", 
                        "PSI",
                        t_list
                    )
                )
                for p in self.psi_over_t.keys():
                    URIs = self.psi_over_t[p]
                    pltxt.write(str(p)+ "\t" + str(list(URIs)))

        else:
            if subplot:
                psiovert_list = list(self.psi_over_t.keys())
                num_labels = self.num_labels
                n_cols = math.floor(math.sqrt(num_labels))
                n_rows = math.ceil(num_labels / n_cols)

                fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
                label_being_plotted = 0
                
                for i in range(n_rows):
                    for j in range(n_cols):
                        p = psiovert_list[label_being_plotted]
                        URIs = self.psi_over_t[p]
                        coeffs = np.polyfit(t_list, URIs, 5)
                        ffit = np.poly1d(coeffs)
                        fderiv = ffit.deriv()
                        URIprime = fderiv(t_list)

                        axs[i,j].plot(t_list, t_list, '--', label='Perfect Reproducibility', linewidth=0.75)
                        axs[i,j].plot(t_list, URIprime, label = 'Derivative Correspondence', linewidth=0.75)
                        axs[i,j].plot(t_list, URIs, label = 'Correspondence', linewidth=0.75)
                        axs[i,j].set_title(str(p), fontsize=7)
                        # axs[i,j].set(xlabel="t", ylabel="PSI")

                        label_being_plotted+=1
                
                plt.tight_layout()
                # plt.legend(loc=0)
                plt.savefig('{}/cc_{}.pdf'.format(self.savedir, "subplot"), format='pdf')
                plt.savefig('{}/cc_{}.svg'.format(self.savedir, "subplot"), format='svg')
                plt.clf()

            for p in self.psi_over_t.keys():
                plt.plot(t_list, t_list, '--', label='Perfect Reproducibility')
                URIs = self.psi_over_t[p]

                coeffs = np.polyfit(t_list, URIs, 5)
                ffit = np.poly1d(coeffs)
                fderiv = ffit.deriv()
                URIprime = fderiv(t_list)
                plt.plot(t_list, URIprime, label = 'Derivative Correspondence')

                plt.plot(t_list, URIs, label = 'Correspondence')
                plt.title('Correspondence Curve for '+str(p))
                plt.xlabel('t')
                plt.ylabel("PSI")
                plt.legend()
                plt.tight_layout()
                plt.savefig('{}/cc_{}.pdf'.format(self.savedir, str(p)), format='pdf')
                plt.savefig('{}/cc_{}.svg'.format(self.savedir, str(p)), format='svg')
                plt.clf()
                with open("{}/cc_{}.txt".format(self.savedir, str(p)), 'w') as pltxt:
                    """
                    title, xaxis, yaxis, x(t_list), uriprime, uri
                    """

                    pltxt.write(
                        "{}\n{}\n{}\n{}\n{}\t{}".format(
                            "Correspondence Curve for " +str(p),
                            "t", 
                            "PSI",
                            t_list,
                            str(list(URIprime)),
                            str(list(URIs))
                        )
                    )

def plot_bidir_bar_chart(metr1, metr2, type, savedir):
    data = pd.DataFrame([metr1, metr2]).transpose()
    data.columns = ["metr1", "metr2"]
    color_red = '#fd625e'
    color_blue = '#01b8aa'
    index = data.index
    column0 = data['metr1']
    column1 = data['metr2']
    title0 = '{}: Annotation 1 vs. Annotation 2'.format(type)
    title1 = '{}: Annotation 2 vs. Annotation 1'.format(type)
    fig, axes = plt.subplots(figsize=(15,10), ncols=2, sharey=True)
    axes[0].barh(index, column0, align='center', color=color_red, zorder=10)
    axes[0].set_title(title0, fontsize=18, pad=15, color=color_red)
    axes[1].barh(index, column1, align='center', color=color_blue, zorder=10)
    axes[1].set_title(title1, fontsize=18, pad=15, color=color_blue)
    axes[0].invert_xaxis() 
    plt.gca().invert_yaxis()
    axes[0].set(yticks=data.index, yticklabels=data.index)
    axes[0].yaxis.tick_left()
    if type == "Raw_Agreement" or type == "Cohens_Kappa":
        axes[1].set_xticks(np.arange(0, 1.1, step=0.1))
        axes[0].set_xticks(np.arange(0, 1.1, step=0.1))
    plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
    fig.tight_layout()
    plt.savefig('{}/bidir_{}.pdf'.format(savedir, type), format='pdf')
    plt.savefig('{}/bidir_{}.svg'.format(savedir, type), format='svg')

    plt.clf()
    sns.reset_orig
    plt.style.use('default')

class sankey(object):
    def __init__(self, loci_1, corrected_loci_2, savedir):
        self.loci_1 = loci_1
        self.loci_2 = corrected_loci_2
        self.num_labels = len(self.loci_1.columns)-3
        self.savedir = savedir

    def sankey_diag(self):
        print('creating sankey diag')
        confmat = confusion_matrix(self.loci_1, self.loci_2, self.num_labels, OE_transform=False)
        label_1 = [i.replace('posterior','replicate1_label') for i in confmat.columns]
        label_2 = [i.replace('posterior','replicate2_label') for i in confmat.columns]
        label = label_1 + label_2

        source = []
        target = []
        value = []  
        for i in range(confmat.shape[0]):
            for j in range(confmat.shape[1]):
                source.append(i)
                target.append(confmat.shape[0] + j)
                value.append(confmat.iloc[i,j])

        color_list = []
        opacity = 0.7
        for i in range(len(confmat.columns)):
            clist = px.colors.qualitative.Prism
            color_list.append(clist[i%len(clist)].replace("rgb","rgba").replace(")", ", {})".format(opacity)))

        color_link = []
        for i in range(len(confmat.columns)):
             color_link = color_link + color_list

        # print(color_link)
        link = dict(source = source, target=target, value=value, color=color_link)
        node = dict(label = label, pad = 10, thickness = 20)
        data = go.Sankey(link=link, node=node)

        fig = go.Figure(data)
        fig.update_layout(
            hovermode = 'y', title = "Agreement Sankey", 
            font=dict(size = 10))
        fig.write_image("{}/sankey.pdf".format(self.savedir))
        fig.write_image("{}/sankey.svg".format(self.savedir))
        fig.write_html("{}/sankey.html".format(self.savedir))
    
    def heatmap(self, new_column=None):
        if new_column != None:
            pass
        else:
            new_column = self.loci_1.columns[3:]
            new_column = self.loci_2.columns[3:]
        confmat = confusion_matrix(self.loci_1, self.loci_2, self.num_labels, OE_transform=True)
        p = sns.heatmap(
            confmat.astype(int), annot=True, fmt="d",
            linewidths=0.01,  cbar=False,
            yticklabels=new_column,
            xticklabels=new_column,
            # vmin=confmat.min().min(),
            # vmax=confmat.max().max()
            )
        sns.set(rc={'figure.figsize':(15,20)})
        p.tick_params(axis='x', rotation=30, labelsize=7)
        p.tick_params(axis='y', rotation=45, labelsize=7)

        plt.title('Label Matching Heatmap (log(O/E) overlap)')
        plt.xlabel('Replicate 1 Labels')
        plt.ylabel("Replicate 2 Labels")
        plt.tight_layout()
        plt.savefig('{}/heatmap.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/heatmap.svg'.format(self.savedir), format='svg')
        plt.clf()
        sns.reset_orig
        plt.style.use('default')

        confmat.to_csv("{}/heatmap.csv".format(self.savedir))
 
class TSS_enrichment(object):
    def __init__(self, loci, TSSdir, savedir):
        self.TSSs = pd.read_csv(TSSdir, sep="\t", header=None)
        self.TSSs.columns = ["chr", "coord", "strand"]
        self.loci = loci
        self.savedir = savedir

    def get_coverage(self):
        MAP = self.loci.iloc[:,3:].idxmax(axis=1)
        self.coverage = dict(zip(list(self.loci.columns[3:]), [0 for _ in range(len(self.loci.columns[3:]))]))

        for c in list(self.loci.columns[3:]):
            self.coverage[c] = len(MAP.loc[MAP == c]) / len(MAP)

    def calculate_overlap_with_tss(self, loci):
        overlaps = []
        for i in range(int(len(self.TSSs))):
            loci_subset = loci.loc[
                (loci["chr"] == self.TSSs["chr"][i]) &
                (loci["start"].astype("int") < int(self.TSSs["coord"][i])) &
                (int(self.TSSs["coord"][i]) < loci["end"].astype("int")), :
            ]

            if len(loci_subset) > 0:
                overlaps.append(loci_subset)

        if len(overlaps) > 0:
            return pd.concat(overlaps, axis=0).reset_index(drop=True)

    def get_overlap(self, m_p, n_p=8):
        if m_p:
            locis = [self.loci.iloc[
                p:p+int(len(self.loci)/n_p) ,:].reset_index(drop=True)
                    for p in range(0, len(self.loci), int(len(self.loci)/n_p))]

            pobj = mp.Pool(n_p)
            overlaps = pobj.map(self.calculate_overlap_with_tss, locis)

                # functools.partial(, 
                # TSSs=self.TSSs), )

            self.overlaps = pd.concat(overlaps, axis=0).reset_index(drop=True)

        else:
            overlaps = []
            
            for i in range(int(len(self.TSSs))):
                loci_subset = self.loci.loc[
                    (self.loci["chr"] == self.TSSs["chr"][i]) &
                    (self.loci["start"] <= int(self.TSSs["coord"][i])) &
                    (int(self.TSSs["coord"][i]) <= self.loci["end"]), :
                ]

                if len(loci_subset) > 0:
                    overlaps.append(loci_subset)
            self.overlaps = pd.concat(overlaps, axis=0).reset_index(drop=True)

        self.overlaps.iloc[:,3:] = self.overlaps.iloc[:,3:].astype("float16")

    def tss_enrich(self, log_transform=True, m_p=True):
        
        self.get_coverage()
        self.get_overlap(m_p=m_p)

        MAP = self.overlaps.iloc[:,3:].idxmax(axis=1)
        enrich_dict = dict(zip(list(self.overlaps.columns[3:]), [0 for _ in range(len(self.overlaps.columns[3:]))]))

        for i in range(len(MAP)):
            enrich_dict[MAP[i]] += 1/len(MAP)

        for k, v in enrich_dict.items():
            if log_transform:
                if v == 0:
                    enrich_dict[k] = np.log(
                        (v + 1/len(MAP)) / (self.coverage[k]+ 1/len(MAP))
                        )
                else:
                    enrich_dict[k] = np.log(
                        (v) / (self.coverage[k])
                        )
                        
            else:
                enrich_dict[k] = (v) / (self.coverage[k])

        with open('{}/enrichment_general.txt'.format(self.savedir),"w") as enr_depo:
                enr_depo.write(str(enrich_dict))

        plt.bar(list(enrich_dict.keys()), list(enrich_dict.values()), color="black", alpha=0.5)
        plt.ylabel("log(O/E) TSS enrichment")
        plt.xticks(rotation=45, fontsize=7)
        plt.tight_layout()
        plt.savefig('{}/general_TSS_enrichment.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/general_TSS_enrichment.svg'.format(self.savedir), format='svg')
        plt.clf()

        return enrich_dict
        
    def tss_enrich_vs_repr(self, num_bins=10, stratified_bins=True, strat_size="def"):
        """
        OPT: calibrate reprod (to get reproducibility scores)
        x-axis: bins of calibrated reproducibility score
        y-axis: enrichment of TSS in that bin
        """
        
        enr_dict = {}       
        if stratified_bins:
            for l in self.loci.columns[3:]:
                if strat_size == "def":
                    strat_size = int(len(self.overlaps)/40)
                
                enr_l = []
                posterior_vector = self.overlaps[l].astype("float").sort_values().reset_index(drop=True)

                for b in range(0, posterior_vector.shape[0], strat_size):
                    subset_vector_pair = posterior_vector.iloc[b:b+strat_size].reset_index(drop=True)
                    bin_range = [
                        float(subset_vector_pair.iloc[subset_vector_pair.index[0]]), 
                        float(subset_vector_pair.iloc[subset_vector_pair.index[-1]])
                    ]
                    _clb = float(len(self.loci.loc[
                        (bin_range[0] <= self.loci[l].astype("float"))&
                        (self.loci[l].astype("float") <= bin_range[1]), l]))

                    _tlb = float(len(self.overlaps.loc[
                        (bin_range[0] <= self.overlaps[l].astype("float"))&
                        (self.overlaps[l].astype("float") <= bin_range[1]), l]))
                    
                    # if _tlb == 0 or _clb == 0:
                    #     # to avoid log crash 0
                    #     _clb += 1/len(self.overlaps)
                    #     _tlb += 1/len(self.overlaps)
                    
                    if _tlb != 0 and _clb != 0:
                        coverage_l_b = float(_clb / (len(self.loci)))
                        tss_l_b = float(_tlb / (len(self.overlaps)))
                        
                        enr_l.append([
                            bin_range[0], bin_range[1], 
                            np.log(tss_l_b/coverage_l_b)]
                            )

                enr_dict[l] = enr_l
            
            colors = [i for i in get_cmap('tab20').colors]

            with open('{}/enrichment_vs_repr.txt'.format(self.savedir), "w") as enr_depo:
                enr_depo.write(str(enr_dict))
            #######################################################################################
            ci = 0
            for k in enr_dict.keys():
                plt.plot(
                    [(kk[0]+kk[1])/2 for kk in enr_dict[k]],
                    [kk[2] for kk in enr_dict[k]], label=k, c=colors[ci]
                )
                ci+=1

            plt.ylabel("log(O/E) TSS enrichment")
            plt.xlabel("Caliberated Posterior")
            plt.title("TSS enrichment VS. Reproducibility")
            plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", ncol=1)
            plt.tight_layout()
            plt.savefig('{}/rep_vs_TSS_enrichment.pdf'.format(self.savedir), format='pdf')
            plt.savefig('{}/rep_vs_TSS_enrichment.svg'.format(self.savedir), format='svg')
            plt.clf()
            
            #######################################################################################

            ci = 0
            
            for k in enr_dict.keys():
                x = [(kk[0]+kk[1])/2 for kk in enr_dict[k]]
                y = [kk[2] for kk in enr_dict[k]]

                try:
                    if len(x)<=3:
                        f = UnivariateSpline(x= np.array(x), y=y, k=1)
                    else:
                        f = UnivariateSpline(x= np.array(x), y=y, k=3)

                    r2 = r2_score(
                        y, 
                        f(
                            np.reshape(np.array(x), (-1,1))
                        ))
                    
                    # print("R2 score for {}: {}".format(k, r2))

                    plt.plot(
                        x, 
                        f(x),
                        label=k, c=colors[ci]
                    )
                    ci+=1
                except:
                    ci+=1
                    print("faced error with {}".format(k))

            plt.ylabel("log(O/E) TSS enrichment")
            plt.xlabel("Caliberated Posterior")
            plt.title("TSS enrichment VS. Reproducibility")
            plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", ncol=1)
            plt.tight_layout()
            plt.savefig('{}/rep_vs_TSS_enrichment_spline.pdf'.format(self.savedir), format='pdf')
            plt.savefig('{}/rep_vs_TSS_enrichment_spline.svg'.format(self.savedir), format='svg')
            plt.clf()
            
            #######################################################################################

            list_binsdict = list(enr_dict.keys())
            num_labels = len(list_binsdict)
            n_cols = math.floor(math.sqrt(num_labels))
            n_rows = math.ceil(num_labels / n_cols)

            fig, axs = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)
            label_being_plotted = 0
            ci = 0
            for i in range(n_rows):
                for j in range(n_cols):
                    k = list_binsdict[label_being_plotted]
                    x = [(kk[0]+kk[1])/2 for kk in enr_dict[k]]
                    y = [kk[2] for kk in enr_dict[k]]
                    
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
                    
            plt.tight_layout()
            plt.savefig('{}/rep_vs_TSS_enrichment_bars.pdf'.format(self.savedir), format='pdf')
            plt.savefig('{}/rep_vs_TSS_enrichment_bars.svg'.format(self.savedir), format='svg')
            plt.clf()
            #######################################################################################

        else:
            for l in self.loci.columns[3:]:
                bin_maps = []
                bin_length = (float(self.loci.iloc[:,3:].max().max()) - float(self.loci.iloc[:,3:].min().min())) / num_bins
                enr_l = []
                for b in range(
                    int(self.loci.iloc[:,3:].min().min()*10000), 
                    int(self.loci.iloc[:,3:].max().max()*10000), 
                    int(bin_length*10000)):

                    bin_range = [float(b/10000), float(b/10000)+bin_length]

                    if len(bin_maps) < num_bins:
                        bin_maps.append(bin_range)

                    # what is the coverage of label l with score in bin range? (in loci)
                    coverage_l_b = float(len(self.loci.loc[
                        (bin_range[0] <= self.loci[l].astype("float"))&
                        (self.loci[l].astype("float") <= bin_range[1]), l])) / len(self.loci)
                    
                    # what is the coverage of label l with score in bin range? (in overlap)
                    tss_l_b = float(len(self.overlaps.loc[
                        (bin_range[0] <= self.overlaps[l].astype("float"))&
                        (self.overlaps[l].astype("float") <= bin_range[1]), l])) / len(self.overlaps)

                    enr_l.append( 
                        np.log(
                            (tss_l_b) / (coverage_l_b)
                        )
                    )

                enr_dict[l] = enr_l

            colors = [i for i in get_cmap('tab20').colors]
            ci = 0
            for k in enr_dict.keys():
                plt.plot([(bin_maps[b][0] + bin_maps[b][1]) /2 for b in range(len(bin_maps))],
                list(enr_dict[k]), label=k, c=colors[ci])
                ci+=1

            plt.ylabel("log(O/E) TSS enrichment")
            plt.xlabel("Caliberated Posterior")
            plt.title("TSS enrichment VS. Reproducibility")
            plt.legend()
            plt.tight_layout()
            plt.savefig('{}/rep_vs_TSS_enrichment.pdf'.format(self.savedir), format='pdf')
            plt.savefig('{}/rep_vs_TSS_enrichment.svg'.format(self.savedir), format='svg')
            plt.clf()
