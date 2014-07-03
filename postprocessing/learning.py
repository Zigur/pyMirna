from __future__ import division
from __future__ import print_function
from model.constants import *
from sklearn import decomposition

__author__ = 'massi'

def get_pca(qpcr_experiment):
    
    data = qpcr_experiment.data[EXP_DELTA_CT]



