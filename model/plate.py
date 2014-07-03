import abc
import pandas as pd

__author__ = 'massi'


class GeneralPlate(object):

    __metaclass__=abc.ABCMeta

    @abc.abstractproperty
    def data(self):
        """property to access the underlying structure of the plate"""
        return

    @abc.abstractmethod
    def normalize(self):
        """given a normalizer object perform a normalization on the plate RAW data"""
        return

    pass


class QPCRPlate(GeneralPlate):

    def __init__(self, data, normalize_func=None):
        self.__data = data
        if normalize_func:
            self.normalize = normalize_func

    @property
    def data(self):
        return self.__data

    def normalize(self):
        pass

    # @property
    # def normalizer(self):
    #     return self.__normalizer
    #
    # @normalizer.setter
    # def normalizer(self, normalizer):
    #     assert isinstance(normalizer, object)
    #     self.__normalizer = normalizer

    @classmethod
    def from_excel(cls, file_path, file_format):
        e_reader = pd.ExcelFile(file_path)
        print(file_format.results_sheet in e_reader.sheet_names, e_reader.sheet_names)
        sheetname = file_format.results_sheet if file_format.results_sheet in e_reader.sheet_names \
            else e_reader.sheet_names[0]
        import pdb; pdb.set_trace()
        data = e_reader.parse(sheetname, header=file_format.column_header,
                              na_values=file_format.na_values)
        data = data[[file_format.well, file_format.well_position, file_format.sample_name, file_format.target_name,
                     file_format.cq, file_format.efficiency]]
        print(data)
        print(file_format.n_wells)
        data = data.ix[0:file_format.n_wells-1]
        data = data.set_index(file_format.index_col)
        e_reader.close()
        return cls(data)

