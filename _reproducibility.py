import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyidr

class IDR(object):
    def __init__(self, loci_1, corrected_loci_2, num_labels, log_transform=True):

        loci_1 = loci_1.loc[:, ['posterior'+str(x) for x in range(int(num_labels))]]
        corrected_loci2 = corrected_loci_2.loc[:, ['posterior'+str(x) for x in range(int(num_labels))]]

        if log_transform:
            self.loci_1 = pd.DataFrame(-1 * np.log(np.array(loci_1)), columns=loci_1.columns)
            self.loci_2 = pd.DataFrame(-1 * np.log(np.array(corrected_loci2)), columns=corrected_loci2.columns)
        else:
            self.loci_1 = loci_1
            self.loci_2 = corrected_loci_2
            
        self.num_labels = int(num_labels)

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

    def fit_idr(self, mean=0.01, sigma=0.001):
        self.reprod_proportion = {}
        self.reproducibility = {}

        for i in range(self.num_labels):
            print(i)
            idr_obj = pyidr.IDR(mean=mean, sigma=sigma)
            idr_obj.fit(np.array(self.loci_1['posterior'+str(i)]), np.array(self.loci_2['posterior'+str(i)]))
            local_idr = idr_obj.predict_global(np.array(self.loci_1['posterior'+str(i)]), np.array(self.loci_2['posterior'+str(i)]))

            self.reprod_proportion['posterior'+str(i)] = idr_obj.proportion
            self.reproducibility['posterior'+str(i)] = local_idr
            
    def rank(self):
        self.loci_1 = pd.DataFrame(-1 * np.argsort(np.array(self.loci_1), axis=1), columns=self.loci_1.columns)
        self.loci_2 = pd.DataFrame(-1 * np.argsort(np.array(self.loci_2), axis=1), columns=self.loci_2.columns)

    def idr_1vs2_plot(self):
        for i in range(self.num_labels):
            plt.scatter(
                self.loci_1['posterior'+str(i)], np.array(self.loci_2['posterior'+str(i)]), 
                s=1, c=self.reproducibility['posterior'+str(i)], cmap='gray')
            plt.title("Labee {}".format(i))
            plt.xlabel('Replicate 1')
            plt.ylabel('Replicate 2')
            plt.colorbar(label='IDR')
            plt.show()
            print(self.reprod_proportion['posterior'+str(i)])
    
    def idr_prop_perlabel(self):
        print(self.reprod_proportion)
        print(self.reprod_proportion.keys())
        print(self.reprod_proportion.values())

        plt.bar(
            x=[i.replace("posterior","") for i in list(self.reprod_proportion.keys())], 
            height=self.reprod_proportion.values(), color='grey')
        plt.show()