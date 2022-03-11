import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d
from _cluster_matching import *
import plotly.graph_objects as go
import plotly.express as px

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

        self.max_posteri1 = self.loci_1.iloc[:, 3:].idxmax(axis=1)
        self.max_posteri2 = self.loci_2.iloc[:, 3:].idxmax(axis=1)

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
            return 
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
        plt.bar(x, height, width=0.4, color='black', alpha=0.5)
        plt.xticks(rotation=45)
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
        plt.bar(x, height, width=0.4, color='black', alpha=0.5)
        plt.xticks(rotation=45)
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
        plt.bar(x, height, width=0.4, color='black', alpha=0.5)
        plt.xticks(rotation=45)
        plt.yticks(np.arange(0, 1.1, step=0.1))
        plt.xlabel('Label')
        plt.title("Cohen's Kappa")
        plt.ylabel("Cohen's Kappa Score")
        plt.savefig('{}/cohenskappa.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/cohenskappa.svg'.format(self.savedir), format='svg')
        plt.clf()


class Reprodroducibility_vs_posterior(object):
    def __init__(self, loci_1, loci_2, savedir, log_transform=True, ignore_overconf=True, filter_nan=True):
        if filter_nan:
            loci_1 = loci_1.dropna()
            loci_1 = loci_1.reset_index(drop=True)
            loci_2 = loci_2.dropna()
            loci_2 = loci_2.reset_index(drop=True)

        self.loci_1 = loci_1
        self.loci_2 = loci_2
        self.savedir = savedir
        self.num_labels = len(self.loci_1.columns)-3
        self.log_transform = log_transform

        self.post_mat_1 = self.loci_1.iloc[:, 3:]
        self.post_mat_2 = self.loci_2.iloc[:, 3:]

        if ignore_overconf:
            max_posteri1 = self.post_mat_1.max(axis=1)
            max_posteri2 = self.post_mat_2.max(axis=1)
            to_drop = []
            for i in range(len(max_posteri1)):
                if max_posteri1[i] == 1 or max_posteri2[i] == 1:
                    to_drop.append(i)

            print("filtered results for overconfident bins: ignored {}/{}".format(len(to_drop), len(self.post_mat_1)))
            self.post_mat_1 = self.post_mat_1.drop(to_drop, axis=0)
            self.post_mat_1 = self.post_mat_1.reset_index(drop=True)

            self.post_mat_2 = self.post_mat_2.drop(to_drop, axis=0)
            self.post_mat_2 = self.post_mat_2.reset_index(drop=True)
        
        self.MAPestimate1 = self.post_mat_1.idxmax(axis=1)
        self.MAPestimate2 = self.post_mat_2.idxmax(axis=1)

        if log_transform:
            self.post_mat_1 = -1 * pd.DataFrame(
                np.log(1 - self.post_mat_1) , 
                columns=self.post_mat_1.columns
            )
            self.post_mat_2 = -1 * pd.DataFrame(
                np.log(1 - self.post_mat_2) , 
                columns=self.post_mat_2.columns
            )

    def per_label_count_independent(self, num_bins=100):
        for k in range(self.post_mat_1.shape[1]):
            try:
                bins = []
                posterior_vector_1 = np.array(self.post_mat_1[self.post_mat_1.columns[k]])
                bin_length = (float(np.max(posterior_vector_1)) - float(np.min(posterior_vector_1))) / num_bins
                for j in range(int(np.min(posterior_vector_1)*1000), int(np.max(posterior_vector_1)*1000), int(bin_length*1000)):
                    bins.append([float(j/1000), float(j/1000) + int(bin_length*1000)/1000, 0, 0, 0, 0, 0])
                    
                    # [bin_start, bin_end, num_values_in_bin, num_agreement_in_bin, num_mislabeled, 
                    # ratio_correctly_labeled, ratio_mislabeled]
                
                for b in range(len(bins)):
                    
                    for i in range(len(posterior_vector_1)):
                        if bins[b][0] <= posterior_vector_1[i] <= bins[b][1]:
                            bins[b][2] += 1

                            if self.post_mat_1.columns[k] == self.MAPestimate2[i]:
                                bins[b][3] += 1
                            else:
                                bins[b][4] += 1

                for b in range(len(bins)):
                    if bins[b][2] != 0:
                        bins[b][5] = bins[b][3] / bins[b][2]
                        bins[b][6] = (bins[b][4]) / bins[b][2]

                bins = np.array(bins)
                plt.scatter(x=bins[:,0], y=bins[:,5], c='black', s=100)
                ysmoothed = gaussian_filter1d(bins[:,5], sigma=3)
                plt.plot(bins[:,0], ysmoothed, c='r')
                plt.title("Reproduciblity Plot label {}".format(k))
                plt.ylabel("Ratio of Similarly Labeled Bins in replicate 2")

                if self.log_transform:
                    plt.xlabel("- log(1-posterior) in Replicate 1")

                else:
                    plt.xlabel("Posterior in Replicate 1")
                
                plt.tick_params(axis='x', which='both', bottom=False,top=False,labelbottom=False)
                plt.savefig('{}/caliberation_label{}.pdf'.format(self.savedir, k), format='pdf')
                plt.savefig('{}/caliberation_label{}.svg'.format(self.savedir, k), format='svg')
                plt.clf()

            except:
                print('could not process label {}'.format(k))

    def general_count_independent(self, num_bins=100):
        # try:
        bins = []
        max_post, min_post = self.post_mat_1.max(axis=1).max(axis=0), self.post_mat_1.min(axis=1).min(axis=0)
        bin_length = (float(max_post) - float(min_post)) / num_bins

        for j in range(int(min_post*1000), int(max_post*1000), int(bin_length*1000)):
            bins.append([float(j/1000), float(j/1000) + int(bin_length*1000)/1000, 0, 0, 0, 0, 0])
            # [bin_start, bin_end, num_values_in_bin, num_agreement_in_bin, num_mislabeled, ratio_correctly_labeled, ratio_mislabeled]
        
        for b in range(len(bins)):
            for k in range(self.post_mat_1.shape[1]):
                for i in range(self.post_mat_1.shape[0]):

                    if bins[b][0] <= self.post_mat_1.iloc[i, k] <= bins[b][1]:
                        bins[b][2] += 1

                        if self.post_mat_1.columns[k] == self.MAPestimate2[i]:
                            bins[b][3] += 1
                        else:
                            bins[b][4] += 1

        for b in range(len(bins)):
            if bins[b][2] != 0:
                bins[b][5] = bins[b][3] / bins[b][2]
                bins[b][6] = (bins[b][4]) / bins[b][2]

        bins = np.array(bins)
        plt.scatter(x=bins[:,0], y=bins[:,5], c='black', s=100)
        ysmoothed = gaussian_filter1d(bins[:,5], sigma=3)
        plt.plot(bins[:,0], ysmoothed, c='r')
        plt.title("Reproduciblity Plot All Labels")
        plt.ylabel("Ratio of Similarly Labeled Bins in replicate 2")

        if self.log_transform:
            plt.xlabel("- log(1-posterior) in Replicate 1")

        else:
            plt.xlabel("Posterior in Replicate 1")
        
        plt.tick_params(axis='x', which='both', bottom=False,top=False,labelbottom=False)
        plt.savefig('{}/caliberation_general.pdf'.format(self.savedir), format='pdf')
        plt.savefig('{}/caliberation_general.svg'.format(self.savedir), format='svg')
        plt.clf()
        
        # except:
        #     print('could not process label {}'.format(k))

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
        self.savedir = savedir
        self.num_labels = len(self.loci_1.columns)-3

        MAPestimate1 = self.loci_1.iloc[\
            :, 3:].idxmax(axis=1)
        MAPestimate2 = self.loci_2.iloc[\
            :, 3:].idxmax(axis=1)

        max_poster1 = self.loci_1.iloc[\
            :, 3:].max(axis=1)
        max_poster2 = self.loci_2.iloc[\
            :, 3:].max(axis=1)

        self.annot_1 = []
        for i in range(loci_1.shape[0]):
            self.annot_1.append([
                loci_1['chr'][i], loci_1['start'][i], loci_1['end'][i], MAPestimate1[i], max_poster1[i]])

        self.annot_2 = []
        for i in range(loci_2.shape[0]):
            self.annot_2.append([
                loci_2['chr'][i], loci_2['start'][i], loci_2['end'][i], MAPestimate2[i], max_poster2[i]])


        self.annot_1 = pd.DataFrame(self.annot_1, columns=['chr', 'start', 'end', 'MAP', 'posterior'])
        self.annot_2 = pd.DataFrame(self.annot_2, columns=['chr', 'start', 'end', 'MAP', 'posterior'])

        self.is_ranked = False

    def rank_general(self):
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
                t = float(t)/num_t_steps
                pl_psi = self.psi(t=t, per_label=True)

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