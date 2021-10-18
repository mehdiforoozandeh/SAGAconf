import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from _utils import *
from scipy.ndimage.filters import gaussian_filter1d

class visualize(object):
    def __init__(self, loci_1, corrected_loci_2, num_labels, savedir=None):
        self.loci_1 = loci_1
        self.loci_2 = corrected_loci_2
        self.num_labels = int(num_labels)
        self.savedir = bool(savedir)
    
    def confusion_matrix_heatmap(self):
        confmat = confusion_matrix(self.loci_1, self.loci_2, self.num_labels)
        labels = [int(i.replace('posterior',''))+1 for i in confmat.columns]
        sns.heatmap(
            confmat.astype(int), annot=True, fmt="d",
            linewidths=0.01, center=0, cbar=False,
            yticklabels=labels,
            xticklabels=labels,
            vmin=0,
            vmax=200)
            
        # plt.colorbar()
        plt.title('Label Matching Heatmap')
        plt.xlabel('Replicate 1 Labels')
        plt.ylabel("Replicate 2 Labels")
        plt.show()

    def sankey_diag(self):
        pass

    def filter_overconf(self):
        # corrected_loci1 = self.loci_1.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]]
        # corrected_loci2 = self.loci_2.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]]

        max_posteri1 = self.loci_1.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]].max(axis=1)
        max_posteri2 = self.loci_2.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]].max(axis=1)

        to_drop = []
        for i in range(len(max_posteri1)):
            if max_posteri1[i] == 1 or max_posteri2[i] == 1:
                to_drop.append(i)

        print("filtered results for overconfident bins: ignored {}/{}".format(len(to_drop), len(self.loci_1)))
        self.loci_1 = self.loci_1.drop(to_drop, axis=0)
        self.loci_1 = self.loci_1.reset_index(drop=True)

        self.loci_2 = self.loci_2.drop(to_drop, axis=0)
        self.loci_2 = self.loci_2.reset_index(drop=True)
    
    def set_threshold(self, margin_of_tolerance=1e-80):
        c_capped = 0
        c_floored = 0
        for i in range(len(self.loci_1)):
            for j in range(self.num_labels):

                if self.loci_1.loc[i, 'posterior'+str(j)] > 1-margin_of_tolerance:
                    self.loci_1.loc[i, 'posterior'+str(j)] = 1.0
                    c_capped +=1

                elif self.loci_1.loc[i, 'posterior'+str(j)] < margin_of_tolerance:
                    self.loci_1.loc[i, 'posterior'+str(j)] = 0.0
                    c_floored +=1

                if self.loci_2.loc[i, 'posterior'+str(j)] > 1-margin_of_tolerance:
                    self.loci_2.loc[i, 'posterior'+str(j)] = 1.0
                    c_capped +=1

                elif self.loci_2.loc[i, 'posterior'+str(j)] < margin_of_tolerance:
                    self.loci_2.loc[i, 'posterior'+str(j)] = 0.0
                    c_floored +=1

        print('number of capped posteriors: {}\nnumber of floored posteriors: {}'.format(c_capped, c_floored))

    def raw_score_histogram(self, n_bins, log_transform=True):
        corrected_loci1 = self.loci_1.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]]
        corrected_loci2 = self.loci_2.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]]

        if log_transform:
            corrected_loci1 = pd.DataFrame(-1 * np.log(np.array(corrected_loci1)), columns=corrected_loci1.columns)
            corrected_loci2 = pd.DataFrame(-1 * np.log(np.array(corrected_loci2)), columns=corrected_loci2.columns)

        for label in corrected_loci1.columns:
            try:
                plt.hist(corrected_loci1[label], bins=n_bins)
                plt.hist(corrected_loci2[label], bins=n_bins)
                plt.title('raw log(posterior) histogram')
                plt.xlabel('log(posterior)')
                plt.ylabel('count')

                if self.savedir:
                    plt.savefig(self.savedir+'/raw_hist_{}.png'.format(label))
                else:
                    plt.show()
            except:
                print('encountered error on label {}'.format(label))

    def raw_scores_plot(self):
        for i in range(self.num_labels):
            fig, axs = plt.subplots(2)
            fig.suptitle('Raw Posterior Scores _ Label {}'.format(i))
            axs[0].plot(range(len(self.loci_1['posterior'+str(i)])), self.loci_1['posterior'+str(i)], c='black')
            axs[1].plot(range(len(self.loci_2['posterior'+str(i)])), self.loci_2['posterior'+str(i)], c='black')

            for ax in axs.flat:
                ax.set(xlabel='Genomic Bin (Locus)', ylabel='Posterior Probability')
                
            plt.show()
    def raw_scores_plot_1vs2_ranked(self, log_transform=True):
        for i in range(self.num_labels):
            if log_transform:
                corrected_loci1 = np.argsort(np.log(1 - np.array(self.loci_1['posterior'+str(i)])))
                corrected_loci2 = np.argsort(np.log(1 - np.array(self.loci_2['posterior'+str(i)])))
                plt.scatter(corrected_loci1, corrected_loci2, s=0.5, c='black')
                plt.title('- log(Posteriors) Rep 1 vs Rep 2 _ label {}'.format(i))

            else:
                plt.scatter(np.argsort(self.loci_1['posterior'+str(i)]), np.argsort(self.loci_2['posterior'+str(i)]), s=0.5, c='black')
                plt.title('Raw Posteriors Rep 1 vs Rep 2 _ label {}'.format(i))
            
            plt.show()

    def raw_scores_plot_1vs2(self, union=False, log_transform=True, oneminusp=False):
        print(self.loci_1.shape)
        for i in range(self.num_labels):
            if i == 11:
                if log_transform:
                    if oneminusp:
                        corrected_loci1 = -1 * np.log(1 - np.array(self.loci_1['posterior'+str(i)]))
                        corrected_loci2 = -1 * np.log(1 - np.array(self.loci_2['posterior'+str(i)]))

                    else:
                        corrected_loci1 = -1 * np.log(np.array(self.loci_1['posterior'+str(i)]))
                        corrected_loci2 = -1 * np.log(np.array(self.loci_2['posterior'+str(i)]))

                    plt.scatter(corrected_loci1, corrected_loci2, s=0.8, c='black')
                    plt.title('Rep 1 vs Rep 2 _ label {}'.format(i))
                    plt.tick_params(axis='x', which='both', bottom=False,top=False,labelbottom=False)
                    plt.yticks([])  # Disable yticks
                    # plt.tick_params(axis='y', which='both', bottom=False,top=False,labelbottom=False)

                else:
                    plt.scatter(self.loci_1['posterior'+str(i)], self.loci_2['posterior'+str(i)], s=0.8, c='black')
                    # plt.title('Raw Posteriors Rep 1 vs Rep 2 _ label {}'.format(i))
                
                if union == False:
                    if oneminusp:
                        plt.xlabel('- log(1-p) Replicate 1')
                        plt.ylabel('- log(1-p) Replicate 2')
                    else:
                        plt.xlabel('- log(p) Replicate 1')
                        plt.ylabel('- log(p) Replicate 2')
                    
                    plt.show()

        if union:
            plt.xlabel('Replicate 1')
            plt.ylabel('Replicate 2')
            plt.show()
        
    def reproducibility_histogram(self, num_bins, Dens=False, log_transform=True):
        corrected_loci1 = self.loci_1.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]]
        corrected_loci2 = self.loci_2.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]]

        maxl1 = corrected_loci1.idxmax(axis=1)
        maxl2 = corrected_loci2.idxmax(axis=1)

        if log_transform:
            corrected_loci1 = pd.DataFrame(-1 * np.log(np.array(corrected_loci1)), columns=corrected_loci1.columns)
            corrected_loci2 = pd.DataFrame(-1 * np.log(np.array(corrected_loci2)), columns=corrected_loci2.columns)

        mis_labeled_list = []

        for i in range(corrected_loci1.shape[0]):
            if maxl1[i] != maxl2[i]:
                mis_labeled_list.append(corrected_loci1[maxl1[i]][i])
        
        print(len(mis_labeled_list)/len(corrected_loci1))
        # mis_labeled_list =[
        #     float(mis_labeled_list[i]/corrected_loci1.shape[0]) for i in range(len(mis_labeled_list))
        #     ]

        plt.hist(mis_labeled_list, bins= num_bins, alpha=0.6, histtype='stepfilled', density=Dens, c='black')

        if Dens:
            plt.ylabel("Probability density of segments with different label in replicate 2")
        else:
            plt.ylabel("Number of segments with different label in replicate 2")
        
        if log_transform:
            plt.xlabel("- Log-posterior Probability ")
        else:
            plt.xlabel("Posterior Probability ")

        plt.title("Irreproduciblity Plot")

        plt.show() 
        plt.clf()

    def count_independent_reproducibility(self, num_bins, irrep=False, log_transform=True):
        corrected_loci1 = self.loci_1.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]]
        corrected_loci2 = self.loci_2.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]]

        max_posteri_label1 = corrected_loci1.idxmax(axis=1)
        max_posteri_label2 = corrected_loci2.idxmax(axis=1)

        max_posteri1 = corrected_loci1.max(axis=1)
        max_posteri2 = corrected_loci2.max(axis=1)

        if log_transform:
            corrected_loci1 = pd.DataFrame(-1 * np.log(1- np.array(corrected_loci1)), columns=corrected_loci1.columns)
            corrected_loci2 = pd.DataFrame(-1 * np.log(1- np.array(corrected_loci2)), columns=corrected_loci2.columns)
            max_posteri1 = np.array(-1 * np.log(1- max_posteri1))
            max_posteri2 = np.array(-1 * np.log(1- max_posteri2))

        bins = []
        bin_length = (float(np.max(max_posteri1)) - float(np.min(max_posteri1))) / num_bins

        for i in range(int(np.min(max_posteri1)*1000), int(np.max(max_posteri1)*1000), int(bin_length*1000)):
            bins.append([float(i/1000), float(i/1000) + int(bin_length*1000)/1000, 0, 0, 0, 0]) 
            # [bin_start, bin_end, num_values_in_bin, num_mislabeled_in_bin, ratio_mislabeled, ratio_correctly_labeled]
        
        # print(bins)

        for i in range(corrected_loci1.shape[0]):
            for j in range(len(bins)):
                if bins[j][0] < float(max_posteri2[i]) < bins[j][1]:
                    bins[j][2] += 1
                
                    if max_posteri_label1[i] != max_posteri_label2[i]:
                        bins[j][3] += 1

        for k in range(len(bins)):
            if bins[k][2] != 0:
                bins[k][4] = bins[k][3] / bins[k][2]
                bins[k][5] = (bins[k][2] - bins[k][3]) / bins[k][2]

        bins = np.array(bins)

        if irrep:
            plt.bar(x=bins[:,0], height=bins[:,4], width=bin_length, alpha=0.8)
            plt.title("Irreproduciblity Plot")
            plt.ylabel("Ratio of Differently Labeled Bins in replicate 2")

        else:
            plt.bar(x=bins[:,0], height=bins[:,5], width=bin_length,alpha=0.8)
            plt.title("Reproduciblity Plot")
            plt.ylabel("Ratio of Similarly Labeled Bins in replicate 2")

        if log_transform:
            plt.xlabel("-Log-posterior Probability in Replicate 1")
        else:
            plt.xlabel("Posterior Probability in Replicate 1")
        
        plt.show()
    
    def per_label_CI_reproducibility(self, num_bins, irrep=False, oneminusp=False, log_transform=True, ignore_overconf=True):
        corrected_loci1 = self.loci_1.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]]
        corrected_loci2 = self.loci_2.loc[:, ['posterior'+str(x) for x in range(self.num_labels)]]

        if ignore_overconf:
            max_posteri1 = corrected_loci1.max(axis=1)
            max_posteri2 = corrected_loci2.max(axis=1)
            to_drop = []
            for i in range(len(max_posteri1)):
                if max_posteri1[i] == 1 or max_posteri2[i] == 1:
                    to_drop.append(i)

            print("filtered results for overconfident bins: ignored {}/{}".format(len(to_drop), len(corrected_loci1)))
            corrected_loci1 = corrected_loci1.drop(to_drop, axis=0)
            corrected_loci1 = corrected_loci1.reset_index(drop=True)

            corrected_loci2 = corrected_loci2.drop(to_drop, axis=0)
            corrected_loci2 = corrected_loci2.reset_index(drop=True)

        max_posteri_label1 = corrected_loci1.idxmax(axis=1)
        max_posteri_label2 = corrected_loci2.idxmax(axis=1)

        if log_transform:
            corrected_loci1 = -1 * pd.DataFrame(np.log(1 - corrected_loci1), columns=corrected_loci1.columns)
            corrected_loci2 = -1 * pd.DataFrame(np.log(1 - corrected_loci2), columns=corrected_loci2.columns)

        for k in range(corrected_loci1.shape[1]):
            try:
                bins = []
                posterior_vector_1 = np.array(corrected_loci1[corrected_loci1.columns[k]])
                posterior_vector_2 = np.array(corrected_loci2[corrected_loci2.columns[k]])
                bin_length = (float(np.max(posterior_vector_1)) - float(np.min(posterior_vector_1))) / num_bins

                for j in range(int(np.min(posterior_vector_1))*1000, int(np.max(posterior_vector_1))*1000, int(bin_length*1000)):
                    bins.append([float(j/1000), float(j/1000) + int(bin_length*1000)/1000, 0, 0, 0, 0, 0])
                    # [bin_start, bin_end, num_values_in_bin, num_agreement_in_bin, num_mislabeled, ratio_correctly_labeled, ratio_mislabeled]
                
                for b in range(len(bins)):
                    
                    for i in range(len(posterior_vector_1)):
                        if bins[b][0] <= posterior_vector_1[i] <= bins[b][1]:
                            bins[b][2] += 1

                            if max_posteri_label2[i] != max_posteri_label1[i]:
                                bins[b][4] += 1
                            else:
                                bins[b][3] += 1

                for b in range(len(bins)):
                    if bins[b][2] != 0:
                        bins[b][5] = bins[b][3] / bins[b][2]
                        bins[b][6] = (bins[b][4]) / bins[b][2]
                
                print(pd.DataFrame(bins))
                bins = np.array(bins)
                if irrep:
                    plt.scatter(x=bins[:,0], y=bins[:,6], c='black', s=100)
                    ysmoothed = gaussian_filter1d(bins[:,6], sigma=3)
                    plt.plot(bins[:,0], ysmoothed, c='r')
                    plt.title("Irreproduciblity Plot label {}".format(k))
                    plt.ylabel("Ratio of Differently Labeled Bins in replicate 2")

                else:
                    plt.scatter(x=bins[:,0], y=bins[:,5], c='black', s=100)
                    ysmoothed = gaussian_filter1d(bins[:,5], sigma=3)
                    plt.plot(bins[:,0], ysmoothed, c='r')
                    plt.title("Reproduciblity Plot label {}".format(k))
                    plt.ylabel("Ratio of Similarly Labeled Bins in replicate 2")

                if log_transform:
                    plt.xlabel("- log(1-posterior) in Replicate 1")
                else:
                    plt.xlabel("Posterior in Replicate 1")
                
                plt.tick_params(axis='x', which='both', bottom=False,top=False,labelbottom=False)
                plt.show()
            except:
                print('could not process label {}'.format(k))
            
