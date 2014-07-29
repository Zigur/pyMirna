from __future__ import division
from __future__ import print_function
from ins.formats import QPCRFileFormat
from model.plate import QPCRPlate
from model.experiment import QPCRExperiment
from model.constants import *
from model.normalizer import execute_global_mean_normalization
import pdb
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

def standard_patients_analysis():
    import time
    import pandas as pd
    from model.normalizer import execute_global_mean_normalization
    from model.experiment import QPCRExperiment
    # from model.constants import *
    from ins.formats import QPCRFileFormat
    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans
    import matplotlib.pyplot as plt
    import numpy as np
    # from mpl_toolkits.mplot3d import Axes3D
    coll = "/home/massi/Projects/miRNA/miRNA-plasmi2013/Dati grezzi mirna esosomi-Plasmi/excel-xls"
    form = QPCRFileFormat.from_text_file("resources/Viia7_TaqMan_Human_MicroRNA_A_B.txt")
    exp = QPCRExperiment.from_collection(coll, form, normalize_func=execute_global_mean_normalization)
    exp.execute_qc_threshold_filter()
    exp.compute_cq_mean()
    exp.normalize()
    exp.data.reset_index(inplace=True)
    exp.data.set_index([SAMPLE_NAME, WELL], inplace=True)
    sample_names = [sample_name for sample_name in exp.data.index.levels[0].values]
    ex = np.asarray([exp.data.loc[sample_name][EXP_DELTA_CT].values for sample_name in sample_names])
    pdb.set_trace()
    pca = PCA()
    pca.fit(ex)
    meaning_dim = 2   # len([el for el in pca.explained_variance_ if el > 0.1])
    pca = PCA(n_components=meaning_dim)
    pca.fit(ex)
    ex = pca.transform(ex)
    k_means = KMeans(n_clusters=3, n_init=12, n_jobs=-1)
    t0 = time.time()
    k_means.fit(ex)
    t_batch = time.time() - t0
    h = 10
    ex_km = k_means.predict(ex)
    pdb.set_trace()
    # Plot the decision boundary. For that, we will assign a color to each
    x_min, x_max = ex[:, 0].min() - 200*h, ex[:, 0].max() + 200*h
    y_min, y_max = ex[:, 1].min() - 200*h, ex[:, 1].max() + 200*h
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))

    # print output to a DataFrame
    df = pd.DataFrame(ex_km, columns=["Cluster"], index=sample_names)
    df["1st Component [X]"], df["2nd Component [Y]"] = ex[:, 0], ex[:, 1]
    df["Centroid [X]"], df["Centroid [Y]"] = [k_means.cluster_centers_[i][0] for i in df["Cluster"]], \
                                             [k_means.cluster_centers_[i][1] for i in df["Cluster"]]
    cols = df.columns.tolist()
    cols = cols[1:3] + cols[0:1] + cols[3:]
    df = df[cols]
    e_writer = pd.ExcelWriter("clustering.xls")
    df.to_excel(e_writer, "k-means", float_format='%.1f')
    e_writer.save()
    # Obtain labels for each point in mesh. Use last trained model.
    z = k_means.predict(np.c_[xx.ravel(), yy.ravel()])
    # pdb.set_trace()
    # Put the result into a color plot
    z = z.reshape(xx.shape)
    fig = plt.figure()
    ax = plt.gca()
    ax.set_xscale('linear')
    ax.set_yscale('linear')
    ax.tick_params(axis='x', colors='k', labelsize=12)
    ax.tick_params(axis='y', colors='k', labelsize=12)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.imshow(z, interpolation='nearest',
              extent=(xx.min(), xx.max(), yy.min(), yy.max()),
              cmap=plt.cm.Paired,
              aspect='auto', origin='lower')
    centroids = k_means.cluster_centers_
    ax.scatter(centroids[:, 0], centroids[:, 1], marker='x', s=169, linewidths=3, color='w', zorder=10)
    ax.plot(ex[:, 0], ex[:, 1], 'k.', markersize=10)
    for i, sample_name in enumerate(sample_names):
        ax.annotate(sample_name.split(' ')[-1], xy=(ex[i, 0], ex[i, 1]), xytext=(3, 3), textcoords='offset points',
                    fontsize=12, fontweight='bold', color='k')
    ax.set_title('Patients Exosomal miRNA - K-Means Clustering', fontweight='bold')
    ax.set_xlabel('First Component', fontweight='bold')
    ax.set_ylabel('Second Component', fontweight='bold')
    print(fig)
    # plt.xlim(x_min, x_max)
    # plt.ylim(y_min, y_max)
    # plt.xticks(())
    # plt.yticks(())
    plt.text(x_min + 25*h, y_min + 25*h,  'train time: %.2fs\ninertia: %.1f' % (t_batch, k_means.inertia_),
             fontweight='bold')
    # plt.show()
    plt.savefig("kmeans.png")
    return z
    # out = np.ndarray(shape=(5, 5))
    # for i, e in enumerate(ex):
    #     for j, x in enumerate(ex):
    #         out[i, j] = np.linalg.norm(e-x)

    # for sample_name in exp.data.index.levels[0].values:
    #     exp_arr.append(exp.data.loc[sample_name][EXP_DELTA_CT].values)
    #     sample_names.append(sample_name)

if __name__ == '__main__':
    main()