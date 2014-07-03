from __future__ import division
import abc
import numpy as np
import pandas as pd

__author__ = 'massi'

class SampleSet(object):
    __metaclass__=abc.ABCMeta

    @abc.abstractmethod
    def n_probes(self):
        """compute the number of indexes contained by the dataset"""
        return

    @abc.abstractmethod
    def n_samples(self):
        """returns the number of samples contained by the SampleSet"""
        return

    @abc.abstractmethod
    def get_sample(self, sample_name):
        """returns a sample given its name"""
        return


class QPCRSampleSet(SampleSet):

    def __init__(self, samples=[], set_name='Sample Set'):
        """initialize a new qPCR sample set"""
        self.__samples = samples
        self.__set_name = set_name