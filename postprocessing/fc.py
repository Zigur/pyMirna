from __future__ import division
from __future__ import print_function
from model.constants import *
import pandas as pd

__author__ = 'massi'


class RegulationExperiment:
    pass


def permute_fold_change_computation(experiment_data, cutoff=32.0):
    def fc_qe(test_row):
        """Evaluates the fold change quality"""
        if test_row.max() < 32.0:
            return FC_QC_OK
        elif test_row.min() < 32.0:
            return FC_QC_A
        else:
            return FC_QC_B
    experiment_data.reset_index(inplace=True)
    experiment_data.set_index([SAMPLE_NAME, WELL], inplace=True)
    sample_names = []
    sample_list = []
    for ind in experiment_data.index.levels[0].values:
        sample_names.append(ind.split()[-1])
        sample_list.append(experiment_data.loc[ind])
    fc_list = []
    for i, sample_ref in enumerate(sample_list):
        for j, sample_test in enumerate(sample_list):
            if i is not j:
                fc_frame = pd.DataFrame(sample_ref[TARGET_NAME])
                # import pdb; pdb.set_trace()
                fc_frame[FOLD_CHANGE + '_' + sample_names[j] + '/' + sample_names[i]] = \
                    sample_test[EXP_DELTA_CT]/sample_ref[EXP_DELTA_CT]

                temp_frame = pd.concat([sample_test[RAW_CT], sample_ref[RAW_CT]], axis=1)
                fc_frame[FC_QUALITY + '_' + sample_names[j] + '/' + sample_names[i]] = temp_frame.apply(fc_qe, axis=1)
                # import pdb; pdb.set_trace()
                fc_list.append({TEST_SAMPLE: sample_names[j], CONTROL_SAMPLE: sample_names[i], FOLD_CHANGE: fc_frame})
    return fc_list


def compute_fold_regulation(fc_data_list):
    """given a fold_change dataset compute the fold regulation
    Fold Regulation is equal to FC if FC >= 1
    and is equal to -1/FC if FC <=1"""
    # if isinstance(fc_data_list, list):
    return [{TEST_SAMPLE: fc_elem[TEST_SAMPLE], CONTROL_SAMPLE: fc_elem[CONTROL_SAMPLE],
             FOLD_REGULATION: fc_elem[FOLD_CHANGE].applymap(lambda x: -1.0/x if x < 1.0 else x)}
            for fc_elem in fc_data_list]
    # return fc_data_list.applymap(lambda x: -1.0/x if x < 1.0 else x)


def filter_nonhsa_targets(data_list):
    """Filter all the targets that are neither humans nor miRNAs
    from a fold change/fold regulation dataset represented as a list"""
    if FOLD_CHANGE in data_list[0].keys():  # check whether it is a FC or a FR dataset
        data_type = FOLD_CHANGE
    elif FOLD_REGULATION in data_list[0].keys():
        data_type = FOLD_REGULATION
    else:
        data_type = 'Results'
    return [{TEST_SAMPLE: data_elem[TEST_SAMPLE], CONTROL_SAMPLE: data_elem[CONTROL_SAMPLE],
             data_type: data_elem[data_type][data_elem[data_type][TARGET_NAME]
                                 .apply(lambda el: True if el.startswith('hsa-') else False)]}
            for data_elem in data_list]


def filter_qc_targets(data_list, discarded=[FC_QC_B], qc_column=2):
    if FOLD_CHANGE in data_list[0].keys():  # check whether it is a FC or a FR dataset
        data_type = FOLD_CHANGE
    elif FOLD_REGULATION in data_list[0].keys():
        data_type = FOLD_REGULATION
    else:
        data_type = 'Results'
    return [{TEST_SAMPLE: data_elem[TEST_SAMPLE], CONTROL_SAMPLE: data_elem[CONTROL_SAMPLE],
             data_type: data_elem[data_type][data_elem[data_type][data_elem[data_type].columns.values[qc_column]].apply(
                 lambda el: False if el in discarded else True)]} for data_elem in data_list]


def get_regulated_targets(data_list, fr_threshold=2.0, fr_column=1):
    """return the upregulated and downregulated samples for each
    fold regulation analysis"""
    if FOLD_CHANGE in data_list[0].keys():  # check whether it is a FC or a FR dataset
        data_type = FOLD_CHANGE
        fr_threshold_down = 1.0/fr_threshold
    elif FOLD_REGULATION in data_list[0].keys():
        data_type = FOLD_REGULATION
        fr_threshold_down = -fr_threshold
    else:
        raise TypeError('The dataframes in data_list  do not contain a Fold Change/Fold Regulation attribute')
    return [{TEST_SAMPLE: data_elem[TEST_SAMPLE], CONTROL_SAMPLE: data_elem[CONTROL_SAMPLE],
             REG_UP: data_elem[data_type][data_elem[data_type][data_elem[data_type].columns
                                                                                   .values[fr_column]] > fr_threshold],
             REG_DOWN: data_elem[data_type][data_elem[data_type][data_elem[data_type].columns
                                                                                     .values[fr_column]] < fr_threshold_down]}
            for data_elem in data_list]


def get_info(data_list):
    """print a set of useful information about a fold change dataset"""
    info = [{TEST_SAMPLE: data_elem[TEST_SAMPLE], CONTROL_SAMPLE: data_elem[CONTROL_SAMPLE],
            N_UP: len(data_elem[REG_UP].index), N_DOWN: len(data_elem[REG_DOWN].index)} for data_elem in data_list]