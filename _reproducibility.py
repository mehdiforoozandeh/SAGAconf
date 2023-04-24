
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
import scipy
from scipy.special import expit

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
    def __init__(self, loci_1, loci_2, window_size, savedir, allow_w=False, plot_raw=True, filter_nan=True, oe_transform=True):
        if filter_nan:
            loci_1 = loci_1.dropna()
            loci_1 = loci_1.reset_index(drop=True)
            loci_2 = loci_2.dropna()
            loci_2 = loci_2.reset_index(drop=True)

        self.loci_1 = loci_1
        self.loci_2 = loci_2
        self.num_labels = len(self.loci_1.columns)-3

        self.resolution = loci_1["end"][0] - loci_1["start"][0]

        self.window_bin = math.ceil(window_size/self.resolution)

        self.enr_ovr = overlap_matrix(loci_1, loci_2, type="IoU")

        self.per_label_matches = {}
        for k in list(loci_1.columns[3:]):
            sorted_k_vector = self.enr_ovr.loc[k,:].sort_values(ascending=False)

            good_matches = [sorted_k_vector.index[0]]
            self.per_label_matches[k] = good_matches[0]

        self.plot_raw = plot_raw
        self.allow_w = allow_w
        self.savedir = savedir
        del loci_1
        del loci_2
        
        self.oe_transform = oe_transform

        self.MAPestimate1 = self.loci_1.iloc[:,3:].idxmax(axis=1)
        self.MAPestimate2 = self.loci_2.iloc[:,3:].idxmax(axis=1)
        self.coverage_1 = {k:len(self.MAPestimate1.loc[self.MAPestimate1 == k]) / len(self.MAPestimate1) for k in self.loci_1.columns[3:]}
        self.coverage_2 = {k:len(self.MAPestimate2.loc[self.MAPestimate2 == k]) / len(self.MAPestimate2) for k in self.loci_2.columns[3:]}

    def perlabel_visualize_calibration(self, bins, label_name, strat_size, scatter=False):
        

        if scatter:
            
            polyreg = IsotonicRegression(
                    y_min=float(np.array(bins[:, 5]).min()), y_max=float(np.array(bins[:, 5]).max()), 
                    out_of_bounds="clip")
            
            polyreg.fit(
                np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)), 
                bins[:, 5])
            
            x = np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])
            y = polyreg.predict(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]))

            if self.plot_raw:
                expected = strat_size * self.coverage_2[self.per_label_matches[label_name]]
                y = (np.exp(y) * expected) / strat_size
                bins[:,5] = (np.exp(bins[:,5]) * expected) / strat_size
                x = expit(x) #sigmoid back from logit
                plt.plot(x, [float(expected/strat_size) for i in range(len(x))], '--', c= "green", linewidth=3)

            plt.scatter(
                x=x, 
                y=bins[:,5], c='black', s=10000/len(bins))
            plt.plot(x, y, '--', c='r', linewidth=3)

        else:
            plt.plot([(bins[i,0]+bins[i,1])/2 for i in range(len(bins))], bins[:,5], label=label_name)

        plt.title("Reproduciblity Plot {}".format(label_name))
        xlabel = "posterior in Replicate 1"

        ylabel = "Similarly Labeled Bins in replicate 2"
       
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

    def general_visualize_calibration(self, bins_dict, strat_size, subplot=True):
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
                    bins = bins_dict[label_name].copy()

                    polyreg = IsotonicRegression(
                        y_min=float(np.array(bins[:, 5]).min()), y_max=float(np.array(bins[:, 5]).max()), 
                        out_of_bounds="clip")
                    
                    polyreg.fit(
                        np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)), 
                        bins[:, 5])
                    
                    x = np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])
                    y = polyreg.predict(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]))
                    r2_y =polyreg.predict(np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)))
                    if self.plot_raw:
                        expected = strat_size * self.coverage_2[self.per_label_matches[label_name]]
                        

                        y = (np.exp(y) * expected ) / strat_size
                        r2_y = (np.exp(r2_y) * expected ) / strat_size
                        bins[:, 5] = (np.exp(bins[:, 5]) * expected ) / strat_size

                        x = expit(x)
                        axs[i,j].plot(x, [float(expected/strat_size) for i in range(len(x))], '--', c="r")


                    r2 = r2_score(bins[:, 5], r2_y)
                        
                    axs[i,j].plot(x, y,c="black")

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
        strat_size = int(len(self.loci_1)/num_bins)

        for k in range(self.num_labels):
            kth_label = self.loci_1.iloc[:,3:].columns[k]
            match_label = self.per_label_matches[kth_label]
            bins = []
                
            posterior_vector_1 = self.loci_1.iloc[:,3:][kth_label]

            vector_pair = pd.concat([posterior_vector_1, self.MAPestimate2], axis=1)

            vector_pair.columns=["pv1", "map2"]

            if self.allow_w:
                vector_pair = vector_pair.sort_values("pv1").reset_index()

            else:
                vector_pair = vector_pair.sort_values("pv1").reset_index(drop=True)

            MAP2 = list(self.MAPestimate2)
            for b in range(0, vector_pair.shape[0], strat_size):
                # [bin_start, bin_end, num_values_in_bin, num_agreement, num_mislabeled, ratio_agreement]
                if self.allow_w:
                    subset_vector_pair = vector_pair.iloc[b:b+strat_size, :]
                    subset_vector_pair = subset_vector_pair.reset_index(drop=True)
                    
                    observed = 0
                    for i in range(len(subset_vector_pair)):
            
                        leftmost = max(0, subset_vector_pair["index"][i] - self.window_bin)
                        rightmost = min(len(MAP2), subset_vector_pair["index"][i] + self.window_bin)

                        neighbors_i = MAP2[leftmost:rightmost]

                        if match_label in neighbors_i:
                            observed += 1
                    
                else:
                    subset_vector_pair = vector_pair.iloc[b:b+strat_size, :]
                    observed = len(subset_vector_pair.loc[subset_vector_pair["map2"] == match_label])
                    
                
                if self.oe_transform:
                    expected = self.coverage_2[match_label] * len(subset_vector_pair)

                    if observed==0:
                        oe = np.log(
                            (observed + 1)/
                            (expected + 1))

                    else:
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
            self.perlabel_visualize_calibration(bins.copy(), kth_label, strat_size, scatter=True)

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

        self.general_visualize_calibration(bins_dict, strat_size)

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

                # print(len(x), len(subset1))
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

class TSS_enrichment(object):
    def __init__(self, loci, TSSdir, savedir):
        self.TSSs = pd.read_csv(TSSdir, sep="\t", header=None)
        self.TSSs.columns = ["chr", "coord", "strand"]
        self.loci = loci
        self.resolution = int(loci["end"][0]) - int(loci["start"][0])
        self.savedir = savedir
        self.MAPestimate = list(self.loci.iloc[:,3:].idxmax(axis=1))

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

            self.overlaps = pd.concat(overlaps, axis=0).reset_index(drop=True)

        else:
            nploci_chr = np.array(self.loci["chr"])
            nploci_start = np.array(self.loci["start"])
            overlaps = []
            t0 = datetime.now()
            nptss = np.array(self.TSSs)
            
            for i in range(int(len(self.TSSs))):
                loci_subset_chr = np.where(nploci_chr == nptss[i, 0])
                loci_subset_start = np.where(
                    (0 <(nploci_start - nptss[i, 1])) &
                     ((nploci_start - nptss[i, 1]) < self.resolution))

                ind_both_conds = np.intersect1d(loci_subset_chr, loci_subset_start)

                if len(ind_both_conds) > 0:
                    loci_subset = self.loci.iloc[ind_both_conds, :]
                    overlaps.append(loci_subset)

            self.overlaps = pd.concat(overlaps, axis=0).reset_index(drop=True)
            # print(self.overlaps)
            print("taking the TSS overlap took: ", datetime.now() - t0)

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
        plt.xticks(rotation=90, fontsize=7)
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
    
    def tss_enr_vs_posterior_rank(self, num_bins=30):
        enr_dict = {}       
        # print(self.overlaps)

        strat_size = int(len(self.overlaps) / num_bins)
        for l in self.loci.columns[3:]:
            
            enr_l = []
            
            """
            split the loci into subsets
            find highest and lowest p values of each subset

            _clb would be the number of bins in that subset
            _tlb would be the number of intersecting bins between that subset and self.overlap
            MAP_frac would be the fraction of bins in the subset for which MAP == l
            
            expected value for all subsets is constant and is len(self.overlap) / len(loci)
            observed is the _tlb / _clb
            """

            posterior_vector = self.overlaps[l].astype("float").sort_values().reset_index(drop=True)
            for b in range(0, posterior_vector.shape[0], strat_size):

                subset_vector_pair = posterior_vector.iloc[b:b+strat_size].reset_index(drop=True)
                bin_range = [
                    float(subset_vector_pair.iloc[subset_vector_pair.index[0]]), 
                    float(subset_vector_pair.iloc[subset_vector_pair.index[-1]])
                ]
                c = self.loci.loc[
                    (bin_range[0] <= self.loci[l].astype("float"))&
                    (self.loci[l].astype("float") <= bin_range[1]), ["chr", "start", "end", l]]

                _clb = float(len(c))

                t = self.overlaps.loc[
                    (bin_range[0] <= self.overlaps[l].astype("float"))&
                    (self.overlaps[l].astype("float") <= bin_range[1]), :]
                
                _tlb = float(len(t))
                
                MAP_frac = len([m for m in c.index if self.MAPestimate[m] == l]) / len(c)
                
                if _tlb != 0 and _clb != 0:
                    coverage_l_b = float(_clb / (len(self.loci)))
                    tss_l_b = float(_tlb / (len(self.overlaps)))
                    
                    enr_l.append([
                        bin_range[0], bin_range[1], 
                        np.log(tss_l_b/coverage_l_b), 
                        MAP_frac])

            enr_dict[l] = pd.DataFrame(enr_l, columns=["left", "right", "enr", "map_frac"])

        list_binsdict = list(enr_dict.keys())
        num_labels = len(list_binsdict)
        n_cols = math.floor(math.sqrt(num_labels))
        n_rows = math.ceil(num_labels / n_cols)

        fig, axs = plt.subplots(n_rows, n_cols, sharex=False, sharey=True, figsize=(15,10))
        label_being_plotted = 0

        for i in range(n_rows):
            for j in range(n_cols):

                l = list_binsdict[label_being_plotted]
                enr = enr_dict[l]
                # enr =  enr.drop_duplicates(subset=['enr'], keep=False).reset_index(drop=True)

                enr_to_plot_MAP = enr.loc[
                    enr["map_frac"]>0, :] # removing bins where MAP != l
                
                enr_to_plot_notMAP = enr.loc[(enr["map_frac"]==0), :]

                axs[i,j].bar(
                    [i for i in enr_to_plot_notMAP.index], enr_to_plot_notMAP.enr, 
                    align='center', width=1, edgecolor='red', color='red', alpha=0.5)

                axs[i,j].bar(
                    [i for i in enr_to_plot_MAP.index], enr_to_plot_MAP.enr, 
                    align='center', width=1, edgecolor='green', color='green', alpha=0.5)

                try:
                    f = UnivariateSpline(x= np.array([i for i in range(len(enr))]), y=enr.enr, k=3)
                    axs[i,j].plot(
                        [i for i in range(len(enr))], 
                        f([i for i in range(len(enr))]),
                        c="black", linewidth=1.5)
                
                except:
                    pass
                
                axs[i,j].plot(
                    [i for i in range(len(enr))], [0 for i in range(len(enr))], 
                    color='black', linewidth=0.5)
                
                axs[i,j].set_xticks([])
                axs[i,j].set_title(l, fontsize=7)
                label_being_plotted +=1

        plt.tight_layout()
        plt.savefig('{}/rep_vs_TSS_enr_posterior_rank.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/rep_vs_TSS_enr_posterior_rank.svg'.format(self.savedir), format='svg')
        plt.clf()

        P_correlations = {}
        S_correlations = {}
        for k in list_binsdict:
            x = np.array([i for i in enr_dict[k].index])
            y = np.array(enr_dict[k]["enr"]) 
            try:
                P_correlations[k] = scipy.stats.pearsonr(x, y)[0]
                S_correlations[k] = scipy.stats.spearmanr(x, y)[0]
            except:
                pass

        ##########################################################################################

        plt.bar(P_correlations.keys(), P_correlations.values(), color="black", alpha=0.5)
        plt.ylabel("Pearson's Correlation")
        plt.xticks(rotation=90, fontsize=8)
        plt.tight_layout()

        plt.savefig(self.savedir+"/poster_tss_pearson_correl.pdf", format='pdf')
        plt.savefig(self.savedir+"/poster_tss_pearson_correl.svg", format='svg')
        
        sns.reset_orig
        plt.close("all")
        plt.style.use('default')

        ##########################################################################################

        plt.bar(S_correlations.keys(), S_correlations.values(), color="black", alpha=0.5)
        plt.ylabel("Spearman's Correlation")
        plt.xticks(rotation=90, fontsize=8)
        plt.tight_layout()

        plt.savefig(self.savedir+"/poster_tss_spearman_correl.pdf", format='pdf')
        plt.savefig(self.savedir+"/poster_tss_spearman_correl.svg", format='svg')
        
        sns.reset_orig
        plt.close("all")
        plt.style.use('default')

        with open(self.savedir+"/poster_tss_spearman_correl.txt", "w") as f:
            f.write(str(S_correlations))

def best_match_overlap(loci1, loci2):
    enr_ovr = overlap_matrix(loci1, loci2, type="IoU")
    best_match_overlap = {i:enr_ovr.loc[i, :].max() for i in enr_ovr.index}
    return best_match_overlap

def nbyn_joint_prob_with_binned_posterior(loci1, loci2, n_bins=10, IoU=False):
    num_labels = len(loci1.columns[3:])

    p_min, p_max = loci1.iloc[:,3:].min().min(), loci1.iloc[:,3:].max().max()
    step_size = float((p_max - p_min)/n_bins)

    p_arange = np.arange(p_min, p_max, step_size)

    xlabels = []
    for i in loci1.columns[3:]:
        for j in p_arange:
            xlabels.append(i + "|" + "{:.2f}".format(j) + "|" + "{:.2f}".format(j+step_size))

    ylabels = []
    for i in loci2.columns[3:]:
        for j in range(n_bins):
            ylabels.append(i + "|" + "{:.2f}".format(j) + "|" + "{:.2f}".format(j+step_size))

    joint = np.zeros((len(xlabels), len(ylabels))) #pd.DataFrame(np.zeros((len(xlabels), len(ylabels))), columns=ylabels, index=xlabels)

    for i in range(joint.shape[0]):
        parsed_i = xlabels[i].split("|")
        r1_label, r1_bin_start, r1_bin_end = parsed_i[0], float(parsed_i[1]), float(parsed_i[2])

        r1_bin_subset = loci1.loc[
            (r1_bin_start <= loci1[r1_label]) & (loci1[r1_label] < r1_bin_end)
            , :]
        
                
        for j in range(joint.shape[1]):
            parsed_j = ylabels[j].split("|")
            r2_label, r2_bin_start, r2_bin_end = parsed_j[0], float(parsed_j[1]), float(parsed_j[2])

            r2_bin_subset = loci2.loc[
                (r2_bin_start <= loci2[r2_label]) & (loci2[r2_label] < r2_bin_end)
                , :]
            
            r1_inds = list(r1_bin_subset.index)
            r2_inds = list(r2_bin_subset.index)

            overlap = len(set(r1_inds).intersection(r2_inds)) / len(loci1)
            
            if IoU:
                if overlap * len(r1_bin_subset) * len(r2_bin_subset) != 0:
                    joint[i, j] = overlap / ((len(r1_bin_subset)/len(loci1)) + (len(r2_bin_subset)/len(loci1)) - overlap)
                else:
                    joint[i, j] = 0

            else:
                joint[i, j] = len(overlap) 

    joint = pd.DataFrame(joint/(num_labels*num_labels), columns=ylabels, index=xlabels)
    
    return joint

def joint_prob_with_binned_posterior(loci1, loci2, n_bins=50, conditional=False, stratified=True):
    num_labels = len(loci1.columns[3:])
    MAP2 = loci2.iloc[:,3:].idxmax(axis=1)
    # coverage2 = {k:len(MAP2.loc[MAP2==k])/len(loci2) for k in loci2.columns[3:]}

    if stratified:
        joint = np.zeros((((n_bins+1) * num_labels), num_labels))
        xlabels = []
        ylabels = list(loci2.columns[3:])

        bi = 0
        strat_size = int(len(loci1)/n_bins)
        for r1_label in loci1.columns[3:]:
            for n in range(n_bins+1):
                xlabels.append("{} - bin #{}".format(r1_label, n+1))

            posterior_vector_1 = loci1.loc[:, r1_label]

            vector_pair = pd.concat([posterior_vector_1, MAP2], axis=1)

            vector_pair.columns=["pv1", "map2"]

            vector_pair = vector_pair.sort_values("pv1").reset_index(drop=True)

            for b in range(0, vector_pair.shape[0], strat_size):
                subset_pair_vector = vector_pair.iloc[b:b+strat_size,:]
                
                for r2_label in range(len(loci2.columns[3:])):
                    intersection_l_r = len(
                        subset_pair_vector.loc[
                        subset_pair_vector["map2"]==loci2.columns[3:][r2_label],:
                        ]) / (len(loci1) * num_labels)
                    
                    P_r = (len(subset_pair_vector) / len(loci1))

                    if conditional:
                        if intersection_l_r * P_r > 0:
                            joint[bi, r2_label] =  (intersection_l_r / P_r ) * num_labels

                        else:
                            joint[bi, r2_label] =  0
                        
                    else:
                        joint[bi, r2_label] =  intersection_l_r 

                bi +=1
        
        joint = pd.DataFrame(joint, columns=ylabels, index=xlabels)
        for ii in range((n_bins * num_labels)+1):
            i = joint.index[ii]
            binnumber = int(i.split("-")[1].replace(" bin #", ""))
            if binnumber > n_bins:
                joint.iloc[ii-1, :] = joint.iloc[ii-1, :] + joint.iloc[ii, :]
                joint = joint.drop(i)


    else:
        p_min, p_max = loci1.iloc[:,3:].min().min(), loci1.iloc[:,3:].max().max()
        step_size = float((p_max - p_min)/n_bins)

        p_arange = np.arange(p_min, p_max, step_size)

        xlabels = []
        for i in loci1.columns[3:]:
            for j in p_arange:
                xlabels.append(i + "|" + str(j) + "|" + str(j+step_size))

        ylabels = []
        for i in loci2.columns[3:]:
            ylabels.append(i)

        joint = np.zeros((len(xlabels), len(ylabels)))

        for i in range(joint.shape[0]):
            parsed_i = xlabels[i].split("|")
            r1_label, r1_bin_start, r1_bin_end = parsed_i[0], float(parsed_i[1]), float(parsed_i[2])

            r1_bin_subset = loci1.loc[(r1_bin_start <= loci1[r1_label]) & (loci1[r1_label] < r1_bin_end), :] 
            r1_bin_subset_MAP2 = MAP2.loc[r1_bin_subset.index].reset_index(drop=True)

            for j in range(joint.shape[1]):
                r2_label = ylabels[j] 
                P_r = (len(r1_bin_subset) / len(loci1))
                # P_l = coverage2[r2_label]
                intersection_l_r = len(r1_bin_subset_MAP2.loc[r1_bin_subset_MAP2 == r2_label]) / (len(loci1) * num_labels)
                
                if conditional:
                    if intersection_l_r * P_r > 0:
                        joint[i,j] = (intersection_l_r / P_r) * num_labels

                    else:
                        joint[i,j] = 0
                    
                else:
                    joint[i,j] = intersection_l_r 
    
        xlabels = []
        for i in loci1.columns[3:]:
            for j in p_arange:
                xlabels.append(i + "|" + "{:.2f}".format(j) + "|" + "{:.2f}".format(j+step_size))
            
        joint = pd.DataFrame(joint, columns=ylabels, index=xlabels)
    
    return joint

def joint_prob_MAP_with_posterior(loci1, loci2, n_bins=50, conditional=False, stratified=True):
    num_labels = len(loci1.columns[3:])
    MAP2 = loci2.iloc[:,3:].idxmax(axis=1)
    MAP1 = loci1.iloc[:,3:].idxmax(axis=1)
    

    joint = np.zeros((((n_bins+1) * num_labels), num_labels))
    xlabels = []
    ylabels = list(loci2.columns[3:])

    bi = 0
    strat_size = int(len(loci1)/n_bins)
    for r1_label in loci1.columns[3:]:
        for n in range(n_bins+1):
            xlabels.append("{} - bin #{}".format(r1_label, n+1))

        posterior_vector_1 = loci1.loc[:, r1_label]

        vector_pair = pd.concat([posterior_vector_1, MAP2, MAP1], axis=1)

        vector_pair.columns=["pv1", "map2", "map1"]

        vector_pair = vector_pair.sort_values("pv1").reset_index(drop=True)

        for b in range(0, vector_pair.shape[0], strat_size):
            subset_pair_vector = vector_pair.iloc[b:b+strat_size,:]
            
            for r2_label in range(len(loci2.columns[3:])):
                intersection_l_r = len(
                    subset_pair_vector.loc[
                    (subset_pair_vector["map2"]==loci2.columns[3:][r2_label])&(subset_pair_vector["map1"] == r1_label),:
                    ]) / (len(loci1) * num_labels)
                
                P_r = (len(subset_pair_vector) / len(loci1))

                if conditional:
                    if intersection_l_r * P_r > 0:
                        joint[bi, r2_label] =  (intersection_l_r / P_r ) * num_labels

                    else:
                        joint[bi, r2_label] =  0
                    
                else:
                    joint[bi, r2_label] =  intersection_l_r 

            bi +=1
    
    joint = pd.DataFrame(joint, columns=ylabels, index=xlabels)
    for ii in range((n_bins * num_labels)+1):
        i = joint.index[ii]
        binnumber = int(i.split("-")[1].replace(" bin #", ""))
        if binnumber > n_bins:
            joint.iloc[ii-1, :] = joint.iloc[ii-1, :] + joint.iloc[ii, :]
            joint = joint.drop(i)
    
    return joint
    
def normalized_mutual_information(loci_1, loci_2, soft=True):
    """
    get raw overlaps  ->  P(A,B) JOINT
    get coverages     ->  P(A), P(B)
    """

    "soft joint prob"
    "soft coverage"
    if soft:
        joint = soft_joint_prob(loci_1, loci_2)
        coverage1, coverage2 = soft_coverage(loci_1, loci_2)
    else:
        joint = joint_overlap_prob(loci_1, loci_2, w=0, symmetric=False)
        MAP1 = loci_1.iloc[:,3:].idxmax(axis=1)
        MAP2 = loci_2.iloc[:,3:].idxmax(axis=1)
        coverage1 = {k:len(MAP1.loc[MAP1 == k])/len(loci_1) for k in loci_1.columns[3:]}
        coverage2 = {k:len(MAP2.loc[MAP2 == k])/len(loci_2) for k in loci_2.columns[3:]}
    
    # entropies
    H_A = 0
    for a in coverage1.keys():
        if coverage1[a] > 0:
            H_A += coverage1[a] * np.log(coverage1[a])

    H_A = -1 * H_A

    H_B = 0
    for b in coverage2.keys():
        if coverage2[b] > 0:
            H_B += coverage2[b] * np.log(coverage2[b])
        
    H_B = -1 * H_B

    # mutual information
    MI = 0
    for a in coverage1.keys():
        for b in coverage2.keys(): 
            
            if (joint.loc[a, b]) != 0:
                MI += joint.loc[a, b] * np.log(
                    (joint.loc[a, b]) / 
                    (coverage1[a] * coverage2[b])
                    )

    NMI = (2*MI)/(H_A + H_B)

    return NMI

def NMI_from_matrix(joint, return_MI=False):
    coverage1 = {k:sum(joint.loc[k,:]) for k in joint.index}
    coverage2 = {k:sum(joint.loc[:,k]) for k in joint.columns}

    # entropies
    H_A = 0
    for a in coverage1.keys():
        if coverage1[a] > 0:
            H_A += coverage1[a] * np.log(coverage1[a])

    H_A = -1 * H_A

    H_B = 0
    for b in coverage2.keys():
        if coverage2[b] > 0:
            H_B += coverage2[b] * np.log(coverage2[b])
        
    H_B = -1 * H_B

    # mutual information
    MI = 0
    for a in coverage1.keys():
        for b in coverage2.keys(): 
            
            if (joint.loc[a, b]) != 0:
                MI += joint.loc[a, b] * np.log(
                    (joint.loc[a, b]) / 
                    (coverage1[a] * coverage2[b])
                    )

    # NMI = (2*MI)/(H_A + H_B)

    # print(MI, H_A, H_B)
    NMI = (MI)/(H_B)

    if return_MI:
        return MI
    else: 
        return NMI

def overlap_heatmap(matrix):
    if matrix.shape[0] <20:
        p = sns.heatmap(
            matrix.astype(float), annot=True, fmt=".2f",
            linewidths=0.01,  cbar=False)

        sns.set(rc={'figure.figsize':(15,20)})
        p.tick_params(axis='x', rotation=30, labelsize=7)
        p.tick_params(axis='y', rotation=30, labelsize=7)

        # plt.title('Overlap Metric')
        plt.xlabel('Replicate 2 Labels')
        plt.ylabel("Replicate 1 Labels")
        plt.tight_layout()
        plt.show()
        plt.clf()
        sns.reset_orig
        plt.style.use('default')

    else:
        p = sns.heatmap(
            matrix.astype(float), annot=False,
            linewidths=0.001,  cbar=True)

        sns.set(rc={'figure.figsize':(20,15)})
        p.tick_params(axis='x', rotation=90, labelsize=7)
        p.tick_params(axis='y', rotation=0, labelsize=7)

        plt.tight_layout()
        plt.show()
        plt.clf()
        sns.reset_orig
        plt.style.use('default')
