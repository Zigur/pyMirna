from __future__ import division
from __future__ import print_function
import abc
from model.constants import *
import numpy as np
import pandas as pd
import eleven


__author__ = 'massi'


def execute_global_mean_normalization(self, cutoff=32.0):
    """Implements Global Mean Normalization
    for each sample in the plate"""
    # samples = self.data.set_index(SAMPLE_NAME, append=True) # add the column SAMPLE_NAME as index
    # samples = self.data.groupby(SAMPLE_NAME)
    sample_list = []
    # for name, sample in samples:
    for ind in self.data.index.levels[0].values:
        sample = self.data.loc[ind]
        # print(len(sample))
        reference_value = sample[sample[PREFILTERED_CT] < cutoff][PREFILTERED_CT].mean(numeric_only=True)
        sample[REFERENCE_VALUE] = reference_value   # pd.Series(np.ones(len(sample))*reference_value)
        sample[DELTA_CT] = sample[CT_MEAN] - sample[REFERENCE_VALUE]
        sample[EXP_DELTA_CT] = (2.0*sample[EFFICIENCY])**-sample[DELTA_CT]
        sample_list.append(sample)
    self.data = pd.concat(sample_list)


def execute_endogenous_control_normalization(self):
    pass