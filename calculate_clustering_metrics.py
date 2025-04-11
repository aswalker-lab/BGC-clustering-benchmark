# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 10:04:28 2025

@author: Allison Walker
This script calculates supervised and unsupervised clustering metrics
"""
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from sklearn import metrics 
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Compares BGC similarity to Tanimoto similarity')
parser.add_argument('bgc_cluster_file', type=str, help='file containing BGC cluster labels, must have a format where BGC id is in the first column and cluster id is in the second column')
parser.add_argument('bgc_score_file', type=str, help='file containing BGC similarity, must have a format where columns one and two are BGC ids and column 3 is a similarity or distance measurement')
parser.add_argument('score_type', type=str, help='Score type, options are dist for distance and sim for similarity',choices=["dist","sim"])
parser.add_argument('-t', '--tanimoto_file', type=str, default='tanimoto_results/NPAtlas_bm_v1.tsv', help='file containing tanimoto comparison of product structures')
parser.add_argument('-c', '--cluster_file', type=str, default='product_clusters/BGC_butina_cluster_label_0.2.txt', help='ground truth cluster file based on structure or user-defined clusters')
parser.add_argument('-o', '--outfile_name', type=str, default='tanimoto_bgc_comparison.txt', help='file to write text results to')

args = parser.parse_args()
infile_name = args.bgc_cluster_file
score_file = args.bgc_score_file
score_type = args.score_type
tanimoto_file = args.tanimoto_file
cluster_file = args.cluster_file
outfile_name = args.outfile_name

outfile = open(outfile_name, 'w')
#load tanimoto data
tanimoto_matrix = pd.read_csv(tanimoto_file)
tanimoto_matrix.set_index('MiBIG_ID', inplace=True)

#load cluster data
labels = {}
for line in open(infile_name):
    split_line = line.split(",")
    bgc_name = split_line[0]
    cluster_label = int(split_line[1])
    labels[bgc_name] = cluster_label
    
#load bgc score data    
scores = {}
max_score = 0
for line in open(score_file):
    split_line = line.split(",")
    bgc1 = split_line[0]
    bgc2 = split_line[1]
    if bgc1 not in labels or bgc2 not in labels:
        continue
    score = float(split_line[3])
    if score > max_score:
        max_score = score
    if bgc1 not in scores:
        scores[bgc1] = {}
    if bgc2 not in scores:
        scores[bgc2] = {}
    #to make score matrix symmetric, take miminimum score if similarity and max if distance
    if bgc2 in scores[bgc1]:
        if score_type == "dist" and score < scores[bgc1][bgc2]:
            scores = scores[bgc1][bgc2]
        if score_type == "sim" and score > scores[bgc1][bgc2]:
            scores = scores[bgc1][bgc2]
    if bgc1 in scores[bgc2]:
        if score_type == "dist" and score < scores[bgc2][bgc1]:
            scores = scores[bgc2][bgc1]
        if score_type == "sim" and score > scores[bgc2][bgc1]:
            scores = scores[bgc2][bgc1]
    scores[bgc1][bgc2] = score
    scores[bgc2][bgc1] = score
    
dist_matrix = np.zeros((len(labels),len(labels)))
tanimoto_matrix_ordered = np.zeros((len(labels),len(labels)))
label_list = []
i = 0
for bgc1 in labels:
    if bgc1 not in scores:
        continue
    label_list.append(labels[bgc1])
    j =0
    for bgc2 in labels:
        tanimoto_matrix_ordered[i][j] = 1- tanimoto_matrix[bgc1[0:bgc1.find(".")]][bgc2[0:bgc2.find(".")]]
        if bgc2 == bgc1:
            dist_matrix[i][j] = 0
            j += 1
            continue
        if bgc2 not in scores[bgc1]:
            dist_matrix[i][j] = 1
            j += 1
            continue
        if score_type == "dist":
            dist_matrix[i][j] = scores[bgc1][bgc2]
        else:
            if max_score > 1:
                dist_matrix[i][j] = 1- scores[bgc1][bgc2]/max_score
            else:
                dist_matrix[i][j] = 1- scores[bgc1][bgc2]       
        j+=1
    i += 1
np.fill_diagonal(dist_matrix, 0)

#Silhouette higher better
sil = metrics.silhouette_score(dist_matrix, np.array(label_list), metric='precomputed')
print("silhouette score: " + str(sil))
outfile.write("silhouette score: " + str(sil) + "\n")
sil = metrics.silhouette_score(tanimoto_matrix_ordered, np.array(label_list), metric='precomputed')
print("tanimoto silhouette score: " + str(sil))
outfile.write("tanimoto silhouette score: " + str(sil) + "\n")

#read in true labels
true_labels = {}
for line in open(cluster_file):
    split_line = line.split(",")
    true_labels[split_line[0]] = int(split_line[1])
    
ordered_true_labels = []
ordered_pred_labels = []
for bgc in labels:
    if bgc not in true_labels:
        continue
    ordered_true_labels.append(true_labels[bgc])
    ordered_pred_labels.append(labels[bgc])

#rand: perfect is 1, in adjusted 0 is random
#The unadjusted Rand index is proportional to the number of sample pairs whose labels are the same in both labels_pred and labels_true, or are different in both.
print("rand score: " + str(metrics.rand_score(ordered_true_labels, ordered_pred_labels)))
outfile.write("rand score: " + str(metrics.rand_score(ordered_true_labels, ordered_pred_labels)) + "\n")
print("adjusted rand score: " + str(metrics.adjusted_rand_score(ordered_true_labels, ordered_pred_labels)))
outfile.write("adjusted rand score: " + str(metrics.adjusted_rand_score(ordered_true_labels, ordered_pred_labels)) + "\n")
#perfect score is 1
print("mutual information: " + str(metrics.adjusted_mutual_info_score(ordered_true_labels, ordered_pred_labels)))
outfile.write("mutual information: " + str(metrics.adjusted_mutual_info_score(ordered_true_labels, ordered_pred_labels)) + "\n")

#v score bounded to 1
print("homogeneity score: " + str(metrics.homogeneity_score(ordered_true_labels, ordered_pred_labels)))
outfile.write("homogeneity score: " + str(metrics.homogeneity_score(ordered_true_labels, ordered_pred_labels)) + "\n")
print("completeness score: " + str(metrics.completeness_score(ordered_true_labels, ordered_pred_labels)))
outfile.write("completeness score: " + str(metrics.completeness_score(ordered_true_labels, ordered_pred_labels)) + "\n")
print("v score: " + str(metrics.v_measure_score(ordered_true_labels, ordered_pred_labels)))
outfile.write("v score: " + str(metrics.v_measure_score(ordered_true_labels, ordered_pred_labels)) + "\n")
#upper bound of 1
print("fowlkes mallow score: " + str(metrics.fowlkes_mallows_score(ordered_true_labels, ordered_pred_labels)))
outfile.write("fowlkes mallow score: " + str(metrics.fowlkes_mallows_score(ordered_true_labels, ordered_pred_labels)) + "\n")
