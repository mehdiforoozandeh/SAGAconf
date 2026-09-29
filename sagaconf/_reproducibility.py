import math
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.special import expit
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import PolynomialFeatures
from ._cluster_matching import overlap_matrix


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
                    if label_being_plotted < len(list_binsdict):
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

                        label_being_plotted += 1
        

            xlabel = "Posterior in Replicate 1"
            if self.oe_transform:
                ylabel = "log(O/E) of Similarly Labeled Bins in replicate 2"
            else:
                ylabel = "log(Ratio) of Similarly Labeled Bins in replicate 2"


            plt.tight_layout()
            plt.savefig('{}/clb_{}.pdf'.format(self.savedir, "subplot"), format='pdf')
            plt.savefig('{}/clb_{}.svg'.format(self.savedir, "subplot"), format='svg')
            plt.clf()

        # colors = [i for i in get_cmap('tab20').colors]
        # ci = 0
        # for label_name, bins in bins_dict.items():
            
        #     polyreg = IsotonicRegression(
        #             y_min=float(np.array(bins[:, 5]).min()), y_max=float(np.array(bins[:, 5]).max()), 
        #             out_of_bounds="clip")
            
        #     polyreg.fit(
        #         np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1)), 
        #         bins[:, 5])
            
        #     r2 = r2_score(
        #             bins[:, 5], 
        #             polyreg.predict(
        #                 np.reshape(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), (-1,1))))
            
        #     plt.plot(
        #         np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))]), 
        #         polyreg.predict(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])), 
        #         label="{}_r2={:.2f}".format(label_name, float(r2)), c=colors[ci])
        #     ci+=1

        # xlabel = "posterior in Replicate 1"
        # if self.oe_transform:
        #     ylabel = "log(O/E) of Similarly Labeled Bins in replicate 2"
        # else:
        #     ylabel = "log(Ratio) of Similarly Labeled Bins in replicate 2"


        # plt.legend(loc='upper center', bbox_to_anchor=(0.45, -0.05),
        #     fancybox=True, ncol=4, fontsize=5)

        # plt.ylabel(ylabel)
        # plt.xlabel(xlabel)
        # plt.tight_layout()
        # plt.savefig('{}/caliberation_{}.pdf'.format(self.savedir, "general"), format='pdf')
        # plt.savefig('{}/caliberation_{}.svg'.format(self.savedir, "general"), format='svg')
        # plt.clf()
        # ################################################
        
        # with open("{}/caliberation_{}.txt".format(self.savedir, "general"), 'w') as pltxt:
        #     """
        #     title, xaxis, yaxis, x, y(all labels)
        #     """
        #     xlabel = "posterior in Replicate 1"
        #     if self.oe_transform:
        #         ylabel = "log(O/E) of Similarly Labeled Bins in replicate 2"
        #     else:
        #         ylabel = "log(Ratio) of Similarly Labeled Bins in replicate 2"


        #     pltxt.write(
        #         "{}\n{}\n{}\n{}".format(
        #             "Reproduciblity Plot {}".format(label_name),
        #             xlabel, 
        #             ylabel,
        #             list(np.array([(bins[i, 0] + bins[i, 1])/2 for i in range(len(bins))])), 
        #         )
        #     )
        #     for label_name, bins in bins_dict.items():
        #         pltxt.write(label_name+ "\t" + str(list(bins[:, 5])) + '\n')
            
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


def NMI_from_matrix(joint, return_MI=False, posterior=False):
    coverage1 = {k:sum(joint.loc[k,:]) for k in joint.index}
    coverage2 = {k:sum(joint.loc[:,k]) for k in joint.columns}

    # entropies
    H_A = 0
    for a in coverage1.keys():
        if coverage1[a] > 0:
            H_A += coverage1[a] * np.log2(coverage1[a])

    H_A = -1 * H_A

    H_B = 0
    for b in coverage2.keys():
        if coverage2[b] > 0:
            H_B += coverage2[b] * np.log2(coverage2[b])
        
    H_B = -1 * H_B

    # mutual information
    MI = 0
    for a in coverage1.keys():
        for b in coverage2.keys(): 
            
            if (joint.loc[a, b]) != 0:
                MI += joint.loc[a, b] * np.log2(
                    (joint.loc[a, b]) / 
                    (coverage1[a] * coverage2[b])
                    )

    # print(MI, H_A, H_B)
    # if posterior:
    NMI = (MI)/(H_B)
    # else:
    #     NMI = (2*MI)/(H_A + H_B)

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
