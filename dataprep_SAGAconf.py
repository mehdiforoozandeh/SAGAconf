# import argparse
import os
import multiprocessing as mp
from _utils import *

# parser = argparse.ArgumentParser()
# args = parser.parse_args()
posterior_dir = "tests/test_posteriors/segway/"
# posteriors = read_posteriordir()

binned_posterior = mp_inplace_binning(posterior_dir, 100, assert_coord_match=False, m_p=False)
print(binned_posterior)