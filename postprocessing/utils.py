from __future__ import division
from __future__ import print_function
import pandas as pd
from model.constants import *

__author__ = 'massi'


def regulated_to_excel(regulated_frame, file_path, float_format='%.3f'):
    """print a dataframe containing up and/or downregulated
    miRNA in test groups with respect to controls. Each page contains a control"""
    e_writer = pd.ExcelWriter(file_path)
    control_groups = regulated_frame.index.levels[0].values
    regulated_frame = regulated_frame.reset_index(drop=False)
    stats = pd.DataFrame()
    common_mirnas = pd.DataFrame()
    for group in control_groups:
        subframe = regulated_frame[regulated_frame[CONTROL_GROUP] == group].sort([TEST_GROUP, REGULATION_TYPE,
                                                                                  FC_QUALITY, FOLD_REGULATION],
                                                                                 ascending=[True, False, False, False],
                                                                                 )
        subframe = subframe.set_index([CONTROL_GROUP, TEST_GROUP, REGULATION_TYPE, FC_QUALITY, WELL_POSITION,
                                       TARGET_NAME], drop=True)
        subframe = subframe.drop(INDEX, 1)
        subframe.reset_index(TARGET_NAME).to_excel(e_writer, sheet_name="CONTROL " + str(group),
                                                   float_format=float_format)
        stats = subframe[FOLD_REGULATION].groupby(level=[CONTROL_GROUP, TEST_GROUP, REGULATION_TYPE]).count() \
            if stats.empty else pd.concat([stats, subframe[FOLD_REGULATION]
                 .groupby(level=[CONTROL_GROUP, TEST_GROUP, REGULATION_TYPE]).count()])
        common_mirnas = subframe.groupby(level=[CONTROL_GROUP, REGULATION_TYPE, FC_QUALITY, WELL_POSITION,
                                                TARGET_NAME]).count() if common_mirnas.empty else \
            pd.concat([common_mirnas, subframe.groupby(level=[CONTROL_GROUP, REGULATION_TYPE,
                                                              FC_QUALITY, WELL_POSITION, TARGET_NAME]).count()])
    pd.DataFrame(stats, columns=['Occurrences']).to_excel(e_writer, sheet_name='STATS', float_format=float_format)
    common_mirnas.reset_index(TARGET_NAME, inplace=True)
    common_mirnas = common_mirnas[common_mirnas[FOLD_REGULATION] > 1].drop(FOLD_CHANGE, 1)
    common_mirnas.columns = [TARGET_NAME, 'Occurrences']
    common_mirnas.sort(columns=['Occurrences'],ascending=False, inplace=True)
    common_mirnas.sort_index(ascending=[True, False, False, True], inplace=True)
    # common_mirnas.sort(columns='Occurrences', ascending=False, inplace=True)
    # common_mirnas.sort(columns=[CONTROL_GROUP, REGULATION_TYPE, FC_QUALITY, 'Occurrences'],
    #                    ascending=[True, False, False, False], inplace=True)
    common_mirnas.to_excel(e_writer, sheet_name='Common miRNA', float_format=float_format)
    e_writer.close()