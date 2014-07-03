from __future__ import print_function
from __future__ import division
import abc

__author__ = 'massi'


class Format(object):

    __metaclass__=abc.ABCMeta


class QPCRFileFormat(Format):

    def __init__(self, format_obj):
        """initialize the fields that define the format of the EXCEL file"""
        self.__chemistry = format_obj['chemistry']
        self.__experiment_name = format_obj['experiment_name']
        self.__instrument_name = format_obj['instrument_name']
        self.__instrument_type = format_obj['instrument_type']
        self.__column_header = int(format_obj['column_header'])
        self.__results_sheet = format_obj['results_sheet']
        self.__index_col = format_obj['well']
        self.__n_wells = int(format_obj['n_wells'])
        self.__na_values = format_obj['na_values']
        self.__well = format_obj['well']
        self.__well_position = format_obj['well_position']
        self.__sample_name = format_obj['sample_name']
        self.__target_name = format_obj['target_name']
        self.__cq = format_obj['Cq']
        self.__efficiency = format_obj['efficiency']

    @property
    def column_header(self):
        return self.__column_header

    @property
    def results_sheet(self):
        return self.__results_sheet

    @property
    def index_col(self):
        return self.__index_col

    @property
    def n_wells(self):
        return self.__n_wells

    @property
    def na_values(self):
        return self.__na_values

    @property
    def well(self):
        return self.__well

    @property
    def well_position(self):
        return self.__well_position

    @property
    def sample_name(self):
        return self.__sample_name

    @property
    def target_name(self):
        return self.__target_name

    @property
    def cq(self):
        return self.__cq

    @property
    def efficiency(self):
        return self.__efficiency

    @classmethod
    def from_text_file(cls, file_name):
        """read a colon separated text file to get the parameters of the
        EXCEL file format"""
        format_dict = dict(line.rstrip('\n').split(':', 1) for line in open(file_name))
        print(format_dict)
        return cls(format_dict)


