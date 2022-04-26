from cProfile import label
from csv import excel
from mimetypes import init
from operator import gt
from re import L, M
from statistics import mean
from textwrap import fill
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d
from _cluster_matching import *
import plotly.graph_objects as go
import plotly.express as px
from sklearn.preprocessing import PolynomialFeatures
from sklearn.isotonic import IsotonicRegression
from sklearn.pipeline import make_pipeline
from scipy.interpolate import interp1d, UnivariateSpline, splrep
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from datetime import datetime


class Agreement(object):
    def __init__(self, loci_1, loci_2, savedir, filter_nan=True):
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

        self.expected_agreement = 1 / float(self.num_labels)

    def per_label_agreement(self):
        label_dict = {}
        label_coverage = {}

        for c in self.loci_1.columns[3:]:
            label_dict[c] = [0, 0] # [agreement_count, disagreement_count]
            label_coverage[c] = 0

        for i in range(self.loci_1.shape[0]):

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
        self.per_label_agreement()
        # l = list(self.label_agreements.values())

        self.overall_agreement = 0

        for l in self.label_agreements.keys():
            self.overall_agreement += self.label_agreements[l] * self.label_coverage[l]

        return self.overall_agreement

    def per_label_OE_ratio(self, log_transform=True):
        self.label_OE = {}
        for k in self.label_agreements.keys():
            if log_transform:
                self.label_OE[k] = np.log(self.label_agreements[k] / (self.expected_agreement + self.epsilon))
            else:
                self.label_OE[k] = self.label_agreements[k] / (self.expected_agreement + self.epsilon)

        return self.label_OE

    def general_OE_ratio(self, log_transform=True):
        if log_transform:
            self.overall_OE = np.log(self.overall_agreement/ (self.expected_agreement + self.epsilon))
        else:
            self.overall_OE =  self.overall_agreement/ (self.expected_agreement + self.epsilon)

        return self.overall_OE  

    def per_label_cohens_kappa(self):
        self.label_CK = {}

        for k in self.label_agreements.keys():
            self.label_CK[k] = \
                (self.label_agreements[k] - (self.expected_agreement + self.epsilon)) / (1 - (self.expected_agreement + self.epsilon))
        
        return self.label_CK

    def general_cohens_kappa(self):
        self.overall_CK = \
            (self.overall_agreement - (self.expected_agreement + self.epsilon)) / (1 - (self.expected_agreement + self.epsilon))
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
        plt.savefig('{}/agreement.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/agreement.svg'.format(self.savedir), format='svg')
        plt.clf()

    def plot_OE(self):
        self.per_label_OE_ratio()
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
        plt.savefig('{}/oe_agreement.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/oe_agreement.svg'.format(self.savedir), format='svg')
        plt.clf()

    def plot_CK(self):
        self.per_label_cohens_kappa()
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
        plt.savefig('{}/cohenskappa.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/cohenskappa.svg'.format(self.savedir), format='svg')
        plt.clf()

class posterior_calibration(object):
    def __init__(self, loci_1, loci_2, savedir, log_transform=True, ignore_overconf=True, filter_nan=True, oe_transform=True):
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
        
        self.log_transform = log_transform
        self.oe_transform = oe_transform

        if ignore_overconf:
            max_posteri1 = self.loci_1.iloc[:,3:].max(axis=1)
            max_posteri2 = self.loci_2.iloc[:,3:].max(axis=1)
            to_drop = []
            for i in range(len(max_posteri1)):
                if max_posteri1[i] == 1 or max_posteri2[i] == 1:
                    to_drop.append(i)

            print("filtered results for overconfident bins: ignored {}/{}".format(len(to_drop), len(self.loci_1)))
            self.loci_1 = self.loci_1.drop(to_drop, axis=0)
            self.loci_1 = self.loci_1.reset_index(drop=True)

            self.loci_2 = self.loci_2.drop(to_drop, axis=0)
            self.loci_2 = self.loci_2.reset_index(drop=True)

        self.MAPestimate1 = self.loci_1.iloc[:,3:].idxmax(axis=1)
        self.MAPestimate2 = self.loci_2.iloc[:,3:].idxmax(axis=1)

        if log_transform:
            self.loci_1.iloc[:,3:] = -1 * pd.DataFrame(
                np.log(1 - self.loci_1.iloc[:,3:]) , 
                columns=self.loci_1.iloc[:,3:].columns
            )
            self.loci_2.iloc[:,3:] = -1 * pd.DataFrame(
                np.log(1 - self.loci_2.iloc[:,3:]) , 
                columns=self.loci_2.iloc[:,3:].columns
            )

    def perlabel_visualize_calibration(self, bins, label_name, scatter=False):
        if scatter:
            plt.scatter(
                x=np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), 
                y=bins[:,5], c='black', s=1000/len(bins))
            
            polyreg = IsotonicRegression(
                    y_min=float(np.array(bins[:, 5]).min()), y_max=float(np.array(bins[:, 5]).max()), 
                    out_of_bounds="clip")
            
            polyreg.fit(
                np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)), 
                bins[:, 5])
            
            plt.plot(
                np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), 
                polyreg.predict(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])), 
                '--', c='r')

        else:
            plt.plot([(bins[i,0]+bins[i,1])/2 for i in range(len(bins))], bins[:,5], label=label_name)

        plt.title("Reproduciblity Plot {}".format(label_name))
        if self.oe_transform:
            plt.ylabel("O/E of Similarly Labeled Bins in replicate 2")
        else:
            plt.ylabel("Ratio of Similarly Labeled Bins in replicate 2")

        if self.log_transform:
            plt.xlabel("- log(1-posterior) in Replicate 1")
            plt.tick_params(axis='x', which='both', bottom=False,top=False,labelbottom=False)

        else:
            plt.xlabel("Posterior in Replicate 1")
       
        plt.savefig('{}/caliberation_{}.pdf'.format(self.savedir, label_name), format='pdf')
        plt.savefig('{}/caliberation_{}.svg'.format(self.savedir, label_name), format='svg')
        plt.clf()

    def general_visualize_calibration(self, bins_dict):
        for label_name, bins in bins_dict.items():
            plt.plot([(bins[i,0]+bins[i,1])/2 for i in range(len(bins))], bins[:, 5], label=label_name)

        plt.legend()
        if self.oe_transform:
            plt.ylabel("o/e reproducibility")
        else:
            plt.ylabel("reproducibility")

        if self.log_transform:
            plt.xlabel("-log(1-p)")
        else:
            plt.xlabel("p")

        plt.savefig('{}/caliberation_{}.pdf'.format(self.savedir, "general"), format='pdf')
        plt.savefig('{}/caliberation_{}.svg'.format(self.savedir, "general"), format='svg')
        plt.clf()
            

    def perlabel_calibration_function(self, method="isoton_reg", degree=3, num_bins=10, return_caliberated_matrix=True):
        perlabel_function = {}
        bins_dict = {}
        new_matrix = [self.loci_1.iloc[:,:3]]
        for k in range(self.num_labels):
            # try:
            kth_label = self.loci_1.iloc[:,3:].columns[k]
            bins = []
            posterior_vector_1 = np.array(self.loci_1.iloc[:,3:][kth_label])

            bin_length = (float(np.max(posterior_vector_1)) - float(np.min(posterior_vector_1))) / num_bins
            for j in range(int(np.min(posterior_vector_1)*1000), int(np.max(posterior_vector_1)*1000), int(bin_length*1000)):
                bins.append([float(j/1000), float(j/1000) + int(bin_length*1000)/1000, 0, 0, 0, 0])
                
                # [bin_start, bin_end, num_values_in_bin, num_agreement_in_bin, num_mislabeled, 
                # ratio_correctly_labeled]
            
            for b in range(len(bins)):
                
                for i in range(len(posterior_vector_1)):
                    if bins[b][0] <= posterior_vector_1[i] <= bins[b][1]:
                        bins[b][2] += 1

                        if kth_label == self.MAPestimate2[i]:
                            bins[b][3] += 1
                        else:
                            bins[b][4] += 1

            for b in bins:
                # if b[2] != 0:
                if self.oe_transform:
                    #enrichment (o/e) correctly labeled
                    expected = ((1/self.num_labels) * b[2]) #(num in bins * 1/numlabels) 
                    b[5] = (b[3] + 1) / (expected + 1)

                else:
                    b[5] = b[3] / b[2] # raw correctly labeled

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
                    make_pipeline(PolynomialFeatures(degree), LinearRegression())
                
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
            return pd.concat(new_matrix, axis=1)

        else:
            self.perlabel_function = perlabel_function
            return self.perlabel_function

    
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
    
    def plot_curve(self, num_t_steps=100, plot_labels=True, plot_general=True, merge_plots=False):
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
            plt.savefig('{}/cc_merged.pdf'.format(self.savedir), format='pdf')
            plt.savefig('{}/cc_merged.svg'.format(self.savedir), format='svg')
            plt.clf()

        else:
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
                plt.savefig('{}/cc_{}.pdf'.format(self.savedir, str(p)), format='pdf')
                plt.savefig('{}/cc_{}.svg'.format(self.savedir, str(p)), format='svg')
                plt.clf()

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
    
    def heatmap(self):
        confmat = confusion_matrix(self.loci_1, self.loci_2, self.num_labels, OE_transform=False)
        print(confmat)
        p = sns.heatmap(
            confmat.astype(int), annot=True, fmt="d",
            linewidths=0.01,  cbar=False,
            yticklabels=self.loci_1.columns[3:],
            xticklabels=self.loci_2.columns[3:],
            # vmin=confmat.min().min(),
            # vmax=confmat.max().max()
            )
        sns.set(rc={'figure.figsize':(15,20)})
        p.tick_params(axis='x', rotation=30, labelsize=7)
        p.tick_params(axis='y', rotation=45, labelsize=7)

        plt.title('Label Matching Heatmap')
        plt.xlabel('Replicate 1 Labels')
        plt.ylabel("Replicate 2 Labels")
        plt.savefig('{}/heatmap.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/heatmap.svg'.format(self.savedir), format='svg')
        plt.clf()

class validate_EXT():
    def __init__(self):
        """
        1. read segtools feature aggregation tab file
        2. read enrichment of each label +- X bp s upstream and downstream of “initial exon start site”
        3. open a df containing these columns
            1. the label
            2. the enrichment
            3. label's reproducibility (agreement)
            4. label's calibrated posterior
        """

    def read_feat_agg_enrichment(self, agg_tab_file, label_agreement_dict, TSS_offset=2):
        aggr_df = []
        with open(agg_tab_file,'r') as aggfile:
            lines = aggfile.readlines()
            for i in range(1, len(lines)):
                sp = lines[i].split('\t')
                sp[-1] = sp[-1].replace("\n", "")
                aggr_df.append(sp)

                if "initial exon" in sp[1] and sp[2]=="0":
                    initial_exon_idx = i-2

        if "E1" in aggr_df[0]:
            aggr_df[0][3:] = [str(i) for i in range(len(aggr_df[0])-3)]

        self.aggr_df = pd.DataFrame(
            aggr_df[1:], columns=aggr_df[0])
        
        TSS = (self.aggr_df.iloc[
            int(initial_exon_idx-(TSS_offset/2)):int(initial_exon_idx+(TSS_offset/2))
            , :])

        print(TSS)
        enrich_df = []
        for k, v in label_agreement_dict.items():
            if "posterior" in k:
                k = k.replace("posterior", "")
            
            enrich_df.append([k, v, mean(list(TSS.loc[:, k].astype('int')))])

        self.enrich_df = pd.DataFrame(enrich_df, columns=["label", "agreement", "enrichment"])

    def TSS_vs_agreement_plot(self):
        plt.scatter(
            x=self.enrich_df["agreement"],
            y=self.enrich_df["enrichment"]
        )
        for i in range(self.enrich_df.shape[0]):
            plt.annotate(
                self.enrich_df['label'][i],
                (
                    self.enrich_df["agreement"][i],
                    self.enrich_df["enrichment"][i]
                )
            )
        plt.xlabel("Agreement")
        plt.xlabel("Label enrichment around TSS")
        plt.show()

    def TSS_from_gtf(self, gtf_file):
        with open(gtf_file, 'r') as gtff:
            lines = gtff.readlines()
            self.TSS_coords = []
            for l in lines:
                if l[0] != "#":
                    cols = l.split(";")
                    cols[0] = cols[0].split("\t")
                    cols = cols[0] + cols[1:]
                    # assure TSS
                    if cols[2] == "gene":
                        if cols[6] == "+":
                            #for + strands TSS=Startsite
                            tss_i = cols[3]
                        elif cols[6] == '-':
                            #for - strands TSS=endsite
                            tss_i = cols[4]
                        self.TSS_coords.append(
                            [cols[0], int(tss_i), cols[10].replace("gene_name ","")[2:-1]])
        self.TSS_coords =  pd.DataFrame(self.TSS_coords, columns=["chr", "TSS", "gene_name"])   
        return self.TSS_coords

    def TSS_enrich_vs_reproducibility(self, loci_1, num_bins=10):
        t0 = datetime.now()
        label_enrich = {}
        num_of_tss_occurence = 0
        bin_length = float(loci_1.iloc[:,3:].max().max() - loci_1.iloc[:,3:].max().min()) / num_bins

        for k in loci_1.iloc[:,3:].columns:
            k_MAP = loci_1.loc[
                (loci_1.iloc[:,3:].idxmax(axis=1) == k), :
            ] 

            bins = {}
            for b in range(int(loci_1.iloc[:,3:].max().min()*10000), int(loci_1.iloc[:,3:].max().max()*10000), int(bin_length*10000)):
                bin_range = [float(b/10000), float(b/10000)+bin_length]

                within_bin = k_MAP.loc[
                    (bin_range[0] <= k_MAP[k]) & (k_MAP[k] <= bin_range[1]), 
                    :]
                
                chrs = np.unique(within_bin[['chr']].values)
                tss_occurence = []

                for c in chrs:

                    chr_c_tss_coords = self.TSS_coords.loc[(self.TSS_coords['chr'] == c), :]
                    chr_c_tss_coords = chr_c_tss_coords.reset_index(drop=True)

                    chr_c_loci = within_bin.loc[(within_bin['chr'] == c), :]
                    chr_c_loci = chr_c_loci.reset_index(drop=True)
                    
                    for i in range(len(chr_c_loci)):

                        tss_occ = chr_c_tss_coords.loc[
                            (chr_c_loci["start"][i] <= chr_c_tss_coords["TSS"]) & 
                            (chr_c_tss_coords["TSS"] <= chr_c_loci["end"][i]),
                            :
                        ]
                        
                        if len(tss_occ) > 0:
                            tss_occurence.append(np.array(chr_c_loci.iloc[i,:]))

                if len(tss_occurence) > 0:
                    tss_occurence = pd.DataFrame(tss_occurence, columns=chr_c_loci.columns)
                    num_of_tss_occurence += len(tss_occurence)

                    bins[tuple(bin_range)] = [len(within_bin), len(tss_occurence)]
            
            label_enrich[k] = bins

        print(label_enrich)

        for k, v in label_enrich.items():
            for kk in v.keys():
                expected_tss_occ = float(v[kk][0])/len(loci_1) * num_of_tss_occurence
                v[kk] = v[kk][1] / expected_tss_occ
        
        print(label_enrich)
        
        for k, v in label_enrich.items():
            plt.plot(
                [(kk[0] + kk[1])/2 for kk in v.keys()],
                [v[kk] for kk in v.keys()],
                label=k
            )

        print(datetime.now() - t0)
        plt.legend()
        plt.xlabel("Calibrated Confidence Score")
        plt.ylabel("O/E # of segment at TSS")
        plt.show()

