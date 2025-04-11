# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 13:35:43 2025

@author: Allison Walker
adapted from https://projects.volkamerlab.org/teachopencadd/talktorials/T005_compound_clustering.html
"""

import pandas as pd
import numpy
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import AllChem 
from rdkit.Chem import DataStructs
import argparse
from sklearn import metrics
import argparse

parser = argparse.ArgumentParser(description='Compares BGC similarity to Tanimoto similarity')
parser.add_argument('-t', '--tanimoto_file', type=str, default='tanimoto_results/NPAtlas_bm_v1.tsv', help='file containing tanimoto comparison of product structures')
parser.add_argument('-c', '--butina_threshold', type=float, default=0.2, help='Butina cluster threshold')
parser.add_argument('-b', '--butina_outfile_name', type=str, default='product_clusters/butina_clusters_pt2.csv', help='file to write butina clusters to')
parser.add_argument('-s', '--scaffolds_outfile_name', type=str, default='product_scaffolds/bm_scaffolds.csv', help='file to write scaffolds to')

args = parser.parse_args()
tanimoto_file = args.tanimoto_file
butina_outfile_name = args.butina_outfile_name
butina_threshold = args.butina_threshold
scaffolds_outfile_name = args.scaffolds_outfile_name

def cluster_fingerprints(distance_matrix, num_cmpds, cutoff=0.2):
    """Cluster fingerprints
    Parameters:
        distance matrix
        cutoff: threshold for the clustering
    """
    # Cluster the data with the Butina algorithm:
    clusters = Butina.ClusterData(distance_matrix, num_cmpds, cutoff, isDistData=True)
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters


tanimoto_matrix = pd.read_csv(tanimoto_file)
tanimoto_matrix.set_index('MiBIG_ID', inplace=True)
tanimoto_distance = 1 - tanimoto_matrix

distance_matrix = tanimoto_distance.values.tolist()

triangular_dist_matrix = []
for i in range(0, len(distance_matrix)):
    for j in range(i):
        triangular_dist_matrix.append(distance_matrix[i][j])
triangular_dist_matrix = numpy.array(triangular_dist_matrix)

# Run the clustering procedure for the dataset
clusters = cluster_fingerprints(triangular_dist_matrix, len(tanimoto_distance), cutoff=butina_threshold)

# Give a short report about the numbers of clusters and their sizes
num_clust_g1 = sum(1 for c in clusters if len(c) == 1)
num_clust_g5 = sum(1 for c in clusters if len(c) > 5)
num_clust_g25 = sum(1 for c in clusters if len(c) > 25)
num_clust_g100 = sum(1 for c in clusters if len(c) > 100)

print("total # clusters: ", len(clusters))
print("# clusters with only 1 compound: ", num_clust_g1)
print("# clusters with >5 compounds: ", num_clust_g5)
print("# clusters with >25 compounds: ", num_clust_g25)
print("# clusters with >100 compounds: ", num_clust_g100)


#write clusters to file
bgc_list = tanimoto_distance.columns.tolist()
outfile = open(butina_outfile_name,'w')

i = 0
bgc_labels = {}
for c in clusters:
    for idx in c:
        outfile.write(bgc_list[idx] + "," + str(i) + "\n")
        bgc_labels[bgc_list[idx]] = i
    i += 1
outfile.close()


label_list = numpy.zeros(len(distance_matrix))
i = 0
for bgc in bgc_list:
    label_list[i] = bgc_labels[bgc]
    i += 1
np_dist_matrix =  numpy.zeros((len(distance_matrix),len(distance_matrix)))
for i in range(0, len(distance_matrix)):
    for j in range(0, len(distance_matrix)):
        np_dist_matrix[i][j] = distance_matrix[i][j]

#calculate silhouette score of clustering
sil = metrics.silhouette_score(np_dist_matrix, numpy.array(label_list), metric='precomputed')
print("Silhouette score: " + str(sil))