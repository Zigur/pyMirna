from __future__ import division
from __future__ import print_function
import abc
import types
import pandas as pd
import numpy as np
from constants import *
import os

__author__ = 'massi'


class GeneralExperiment(object):

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def add_plates(self, *plates):
        """add a new plate to the experiment"""
        return

    @abc.abstractmethod
    def remove_plate(self):
        """remove a plate from the experiment"""
        return

    @abc.abstractmethod
    def normalize(self):
        """given a normalizer object perform a normalization on the plate RAW data"""
        return

    @abc.abstractmethod
    def to_excel(self):
        "save the experiment data on an output stream (i.e. excel)"
        return


class QPCRExperiment(GeneralExperiment):

    def __init__(self, data=None, normalize_func=None, prefilters=None):
        if data is not None:
            self.__data = data
            self.__n_plates = data[PLATE_NAME].nunique()
            self.__data.reset_index(inplace=True)
            self.__data.set_index([PLATE_NAME, WELL], inplace=True)
        else:
            self.__data = pd.DataFrame()
            self.__n_plates = 0
        if normalize_func:
            self.normalize = types.MethodType(normalize_func, self, QPCRExperiment)
        if prefilters:
            self.__prefilters = prefilters

    def __str__(self):
        return self.__data.to_string()

    @property
    def data(self):
        return self.__data

    @data.setter
    def data(self, data):
        self.__data = data

    @property
    def n_plates(self):
        return self.__n_plates

    def add_plates(self, *plates):
        """add one or more plates to the experiment"""
        plates_data = [pd.concat([plate.data,
                                  pd.DataFrame(np.ones(len(plate.data.index))*(self.__n_plates + i),
                                               index=np.arange(1, len(plate.data.index)+1), columns=[PLATE_NAME])], 1)
                       for i, plate in enumerate(plates)]
        if self.__data.empty:
            self.__data = pd.concat(plates_data)
            self.__data.reset_index(inplace=True)
            self.__data.set_index([PLATE_NAME, WELL], inplace=True)
        else:
            self.__data = pd.concat([self.__data] + plates_data)
        self.__n_plates += len(plates_data)

    def remove_plate(self, plate_index):
        del self.__plates[plate_index]

    def execute_qc_threshold_filter(self, cutoff=32.0, discard=False):
        """filters all wells/probes that are above a specific (aka pre-determined)
        CT/Cq value cutoff value. Tipical cutoff values are 30, 32 or 35."""
        def threshold_filter(data_row):
            if data_row[PREFILTERED_CT] < cutoff:
                return data_row[PREFILTERED_CT]
            else:
                return cutoff if discard is False else np.NaN
        self.__data[PREFILTERED_CT] = self.__data[RAW_CT]
        self.__data[PREFILTERED_CT] = self.__data.apply(threshold_filter, 1)

    def compute_cq_mean(self):
        if PREFILTERED_CT not in self.__data.columns.values:
            self.__data[PREFILTERED_CT] = self.__data[RAW_CT]
        samples = self.__data.groupby(SAMPLE_NAME)
        res = []
        for name, sample in samples:
            ct_mean = sample[[TARGET_NAME, PREFILTERED_CT]].groupby(TARGET_NAME).mean()
            ct_std = sample[[TARGET_NAME, PREFILTERED_CT]].groupby(TARGET_NAME).std()
            sample[CT_MEAN] = sample.apply(lambda series: ct_mean.ix[series[TARGET_NAME]][PREFILTERED_CT], axis=1)
            sample[CT_STD] = sample.apply(lambda series: ct_std.ix[series[TARGET_NAME]][PREFILTERED_CT], axis=1)
            res.append(sample)
        self.__data = pd.concat(res)

    def normalize(self):
        pass

    def to_excel(self, e_reader, **kwd):
        self.__data.to_excel(e_reader, **kwd)

    def trial(self):
        pass

    @classmethod
    def from_collection(cls, root_dir, file_format, **kwd):
        data_set = []
        file_list = []
        for dir_name, dir_names, file_names in os.walk(root_dir):
            for file_name in file_names:
                if not file_name.startswith('.'):
                    file_list.append(os.path.join(dir_name, file_name))
        for i, file_elem in enumerate(file_list):
            try:
                e_reader = pd.ExcelFile(file_elem)
                sheet = file_format.results_sheet if file_format.results_sheet in e_reader.sheet_names else 0
                # import pdb; pdb.set_trace()
                data = e_reader.parse(sheet, header=file_format.column_header,
                                      na_values=file_format.na_values)
                data = data[[file_format.well, file_format.well_position,
                             file_format.sample_name, file_format.target_name,
                             file_format.cq, file_format.efficiency]].ix[0:file_format.n_wells-1]
                data = data.set_index(file_format.index_col)
                data_set.append(data)
                data_set[-1][PLATE_NAME] = len(data_set)-1
                e_reader.close()
            except ImportError:
                print("Error importing data from file ", file_elem)
        return cls(data=pd.concat(data_set), **kwd)


class FoldRegulationExperiment(GeneralExperiment):

    def __init__(self, data=None):
        self.__data = data.reset_index(drop=False).set_index([CONTROL_GROUP, TEST_GROUP, INDEX])

    def compute_fold_regulation(self):
        """given a fold_change dataset compute the fold regulation
        Fold Regulation is equal to FC if FC >= 1
        and is equal to -1/FC if FC <=1"""
        # if isinstance(fc_data_list, list):
        fr_list = [data_elem[FOLD_CHANGE].applymap(lambda x: -1.0/x if x < 1.0 else x)
                   for i, data_elem in enumerate(self.__data)]

        # return fc_data_list.applymap(lambda x: -1.0/x if x < 1.0 else x)

    @property
    def data(self):
        return self.__data

    def add_plates(self, *plates):
        pass

    def remove_plate(self):
        pass

    def normalize(self):
        raise NotImplementedError("You cannot normalize a Fold Regulation Dataset")

    def to_excel(self, e_reader, **kwd):
        self.__data.to_excel(e_reader, **kwd)

    def get_diffregulated(self, fr_threshold=10.0, qc_filter=[FC_QC_OK, FC_QC_A], keep_hsa_only=True):
        res = self.__data[((self.__data[FOLD_REGULATION] >= fr_threshold) |
                           (self.__data[FOLD_REGULATION] <= -fr_threshold)) &
                           (self.__data[FC_QUALITY].apply(lambda elem: elem in qc_filter))]
        if keep_hsa_only:
            res = res[res[TARGET_NAME].apply(lambda el: True if el.startswith('hsa-') else False)]
        res[REGULATION_TYPE] = res[FOLD_REGULATION].apply(lambda elem: REG_UP if elem >= fr_threshold else REG_DOWN)
        return res

    def get_upregulated(self, **kwd):
        res = self.get_diffregulated(**kwd)
        return res[res[REGULATION_TYPE] == REG_UP]

    def get_downregulated(self, **kwd):
        res = self.get_diffregulated(**kwd)
        return res[res[REGULATION_TYPE] == REG_DOWN]

    @classmethod
    def from_qpcrexperiment(cls, qpcr_experiment, control_groups=None, qe_threshold=NORMALIZATION_CUTOFF):
        """Build a FoldRegulationExperiment from a QPCRExperiment"""
        def fc_qe(test_row):
            """Evaluates the fold change quality"""
            if test_row.max() < qe_threshold:
                return FC_QC_OK
            elif test_row.min() < qe_threshold:
                return FC_QC_A
            else:
                return FC_QC_B
        qpcr_experiment.data.reset_index(inplace=True)
        if GROUP_NAME not in qpcr_experiment.data.columns.values:
            qpcr_experiment.data[GROUP_NAME] = qpcr_experiment.data[SAMPLE_NAME]
        qpcr_experiment.data.set_index([GROUP_NAME, WELL], inplace=True)
        group_names = []
        groups = []
        for ind in qpcr_experiment.data.index.levels[0].values:
            group_names.append(ind.split()[-1])
            groups.append(qpcr_experiment.data.loc[ind])
        fc_frame = pd.DataFrame()
        if control_groups is None:
            control_groups = groups
        for i, control_group in enumerate(control_groups):
            for j, test_group in enumerate(groups):
                if i is not j:
                    in_frame = pd.DataFrame(control_group[[WELL_POSITION, TARGET_NAME]])
                    # import pdb; pdb.set_trace()
                    in_frame[FOLD_CHANGE] = test_group[EXP_DELTA_CT]/control_group[EXP_DELTA_CT]
                    temp_frame = pd.concat([test_group[RAW_CT], control_group[RAW_CT]], axis=1)
                    in_frame[FC_QUALITY] = temp_frame.apply(fc_qe, axis=1)
                    # import pdb; pdb.set_trace()
                    # fc_list.append({TEST_SAMPLE: group_names[j], CONTROL_SAMPLE: group_names[i],
                    # FOLD_CHANGE: fc_frame})
                    in_frame[TEST_GROUP] = group_names[j]
                    in_frame[CONTROL_GROUP] = group_names[i]
                    fc_frame = pd.concat([fc_frame, in_frame])
        fc_frame[FOLD_REGULATION] = fc_frame[FOLD_CHANGE].apply(lambda x: -1.0/x if x < 1.0 else x)
        return cls(fc_frame[[WELL_POSITION, TARGET_NAME, TEST_GROUP, CONTROL_GROUP, FOLD_CHANGE, FOLD_REGULATION,
                             FC_QUALITY]])