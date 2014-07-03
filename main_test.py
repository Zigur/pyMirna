from __future__ import division
from __future__ import print_function
from ins.formats import QPCRFileFormat
from model.plate import QPCRPlate
from model.experiment import QPCRExperiment
from model.constants import *
from model.normalizer import execute_global_mean_normalization

from os import listdir
from os.path import isfile, join, dirname

__author__ = 'massi'


def main():
    in_path = '/home/massi/Projects/miRNA-Exo-Project/AnalisiGlico/Grezzi/'
    out_path = '/home/massi/Projects/miRNA-Exo-Project/AnalisiGlico/Elaborati/'
    file_name_1 = '2014-03-06 Plasma Glico campione2-ViiA7-export.xls'
    file_name_2 = '2014-03-06 Plasma glico 1964 2ml-ViiA7-export.xls'

    viia7_format = QPCRFileFormat.from_text_file(join(dirname(__file__),
                                                              'resources/Viia7_TaqMan_Human_MicroRNA_A_B.txt'))

    plate_1 = QPCRPlate.from_excel(in_path + file_name_1, viia7_format)
    plate_2 = QPCRPlate.from_excel(in_path + file_name_2, viia7_format)
    experiment = QPCRExperiment(normalize_func=execute_global_mean_normalization)
    experiment.add_plates(plate_1, plate_2)
    res = experiment.compute_cq_mean()
    # print(res)
    experiment.normalize()
    # print(experiment.data.tail(10))
    experiment.to_excel(out_path + 'prova-new.xls')
    return experiment

if __name__ == '__main__':
    main()