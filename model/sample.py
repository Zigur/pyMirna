import abc
import numpy as np
import pandas as pd
from constants import *
__author__ = 'massi'


class Sample(object):

    __metaclass__=abc.ABCMeta

    @abc.abstractmethod
    def size(self):
        """returns the number of probes contained in your sample"""
        return

    @abc.abstractmethod
    def mean(self):
        """returns the mean value of your sample"""
        return

    @abc.abstractmethod
    def std(self):
        """returns the standard deviation in your sample"""
        return


class QPCRSample(Sample):

    def __init__(self, data, sample_name):
        self.__data = data
        self.__sample_name = sample_name

    @classmethod
    def from_data_frame(cls, data_frame, sample_name=SAMPLE_NAME):
        return cls(data_frame, sample_name)

    @classmethod
    def from_excel(cls, file_path, file_format):
        e_reader = pd.ExcelFile(file_path)
        data = e_reader.parse(file_format.results_sheet, header=file_format.column_header)
        e_reader.close()

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, data):
        self.__data = data

    @data.deleter
    def data(self):
        del self.__data

    @property
    def sample_name(self):
        return self.__sample_name

    @data.setter
    def sample_name(self, sample_name):
        self.__sample_name = sample_name

    @data.deleter
    def data(self):
        del self.__sample_name

    def size(self):
        return len(self.__data)

    def mean(self):
        return np.nanmean(self.__data)

    def std(self):
        return np.nanstd(self.__data)
