import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Agreement(object):
    def __init__(self, loci_1, loci_2):
        self.loci_1 = loci_1
        self.loci_2 = loci_2

    def per_label_agreement(self):
        pass

    def general_agreement(self):
        pass

    def per_label_cohens_kappa(self):
        pass

    def general_cohens_kappa(self):
        pass

    def per_label_OE_ratio(self, log_transform=True):
        pass

    def general_OE_ratio(self, log_transform=True):
        pass


class Reprodroducibility_vs_posterior(object):
    def __init__(self, loci_1, loci_2):
        self.loci_1 = loci_1
        self.loci_2 = loci_2

    def per_label_count_independent():
        pass

    def general_count_independent():
        pass

class correspondence_curve(object):
    def __init__(self, loci_1, loci_2):
        self.loci_1 = loci_1
        self.loci_2 = loci_2

    def rank():
        pass
    
    def psi():
        pass

    def psi_prime():
        ''' derivation of psi '''
        pass