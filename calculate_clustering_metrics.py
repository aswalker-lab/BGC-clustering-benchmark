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
from rdkit import Chem

parser = argparse.ArgumentParser(description='Compares BGC similarity to Tanimoto similarity')
parser.add_argument('bgc_cluster_file', type=str, help='file containing BGC cluster labels, must have a format where BGC id is in the first column and cluster id is in the second column')
parser.add_argument('bgc_score_file', type=str, help='file containing BGC similarity, must have a format where columns one and two are BGC ids and column 3 is a similarity or distance measurement')
parser.add_argument('score_type', type=str, help='Score type, options are dist for distance and sim for similarity',choices=["dist","sim"])
parser.add_argument('-t', '--tanimoto_file', type=str, default='tanimoto_results/NPAtlas_bm_v1.tsv', help='file containing tanimoto comparison of product structures')
parser.add_argument('-c', '--cluster_file', type=str, default='product_clusters/butina_clusters_pt2.csv', help='ground truth cluster file based on structure or user-defined clusters')
parser.add_argument('-s', '--scaffold_file', type=str, default='product_scaffolds/bm_scaffolds.csv', help='file with compound scaffolds')
parser.add_argument('-o', '--outfile_name', type=str, default='tanimoto_bgc_comparison.txt', help='file to write text results to')
parser.add_argument('-p', '--outfile_histogram', type=str, default='scaffold_histogram.png', help='file to write histogram to')

args = parser.parse_args()
infile_name = args.bgc_cluster_file
score_file = args.bgc_score_file
score_type = args.score_type
tanimoto_file = args.tanimoto_file
cluster_file = args.cluster_file
scaffold_file = args.scaffold_file
outfile_name = args.outfile_name
outfile_histogram = args.outfile_histogram

outfile = open(outfile_name, 'w')
#load tanimoto data
tanimoto_matrix = pd.read_csv(tanimoto_file)
tanimoto_matrix.set_index('MiBIG_ID', inplace=True)

#load cluster data
labels = {}
clusters = {}
for line in open(infile_name):
    split_line = line.split(",")
    bgc_name = split_line[0]
    cluster_label = int(split_line[1])
    if cluster_label not in clusters:
        clusters[cluster_label] = []
    clusters[cluster_label].append(bgc_name)
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
    score = float(split_line[2])
    if score > max_score:
        max_score = score
    if bgc1 not in scores:
        scores[bgc1] = {}
    if bgc2 not in scores:
        scores[bgc2] = {}
    #to make score matrix symmetric, take miminimum score if similarity and max if distance
    if bgc2 in scores[bgc1]:
        if score_type == "dist" and score < scores[bgc1][bgc2]:
            scores[bgc1][bgc2] = scores[bgc1][bgc2]
        if score_type == "sim" and score > scores[bgc1][bgc2]:
            scores[bgc1][bgc2] = scores[bgc1][bgc2]
    if bgc1 in scores[bgc2]:
        if score_type == "dist" and score < scores[bgc2][bgc1]:
            scores[bgc1][bgc2] = scores[bgc2][bgc1]
        if score_type == "sim" and score > scores[bgc2][bgc1]:
            scores[bgc1][bgc2] = scores[bgc2][bgc1]
    scores[bgc1][bgc2] = score
    scores[bgc2][bgc1] = score
    
new_labels = {}
for bgc1 in labels:
    if bgc1 in tanimoto_matrix:
        new_labels[bgc1]= labels[bgc1]
labels = new_labels

dist_matrix = np.zeros((len(labels),len(labels)))
tanimoto_matrix_ordered = np.zeros((len(labels),len(labels)))
label_list = []
i = 0
for bgc1 in labels:
    if bgc1 not in scores or bgc1 not in tanimoto_matrix:
        continue
    label_list.append(labels[bgc1])
    j =0
    for bgc2 in labels:
        if bgc2 not in scores or bgc2 not in tanimoto_matrix:
            continue
        tanimoto_matrix_ordered[i][j] = 1- tanimoto_matrix[bgc1][bgc2]
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

#analyze scaffolds
scaffolds = {}
for line in open(scaffold_file):
    split_line= line.split(",")
    bgc = split_line[0]
    scaffold = split_line[1].replace("\n","")
    canon_scaffold = Chem.CanonSmiles(scaffold)
    #skip empty scaffolds
    if canon_scaffold == "":
        continue
    if bgc not in scaffolds:
        scaffolds[bgc] = [canon_scaffold]
        continue
    #check if already have scaffold
    have_scaffold = False
    for smiles in scaffolds[bgc]:
        if canon_scaffold == smiles:
            have_scaffold = True
            break
    if not have_scaffold:
        scaffolds[bgc].append(canon_scaffold)
        
gcf_scaffolds = {}
gcf_max_bgc_scaffold = {}
gcf_scaffolds_normalized_counts = {}
avg_scaffolds_per_gcf = 0
for gcf in clusters:
    gcf_scaffolds[gcf] = []
    gcf_max_bgc_scaffold[gcf] = 0
    for bgc in clusters[gcf]:
        if bgc not in scaffolds:
            continue
        for scaffold in scaffolds[bgc]:
            if scaffold not in gcf_scaffolds[gcf]:
                gcf_scaffolds[gcf].append(scaffold)
        if len(scaffolds[bgc]) > gcf_max_bgc_scaffold[gcf]:
            gcf_max_bgc_scaffold[gcf] = len(scaffolds[bgc])
    if gcf_max_bgc_scaffold[gcf] > 0:
        gcf_scaffolds_normalized_counts[gcf] = len(gcf_scaffolds[gcf])/gcf_max_bgc_scaffold[gcf]
    else:
        gcf_scaffolds_normalized_counts[gcf] = 1
    avg_scaffolds_per_gcf += gcf_scaffolds_normalized_counts[gcf]

print("average normalized scaffolds per gcf: " + str(avg_scaffolds_per_gcf/len(clusters)))
outfile.write("average normalized scaffolds per gcf: " + str(avg_scaffolds_per_gcf/len(clusters)) + "\n")
#make histogram
number_list = []
fig, axs = plt.subplots(2, 1) # 2 rows, 1 column
for gcf in gcf_scaffolds_normalized_counts:
    number_list.append(gcf_scaffolds_normalized_counts[gcf])
axs[0].hist(number_list)
axs[0].set_xlabel("Normalized Scaffolds per GCF",fontsize=16)
axs[0].set_ylabel("Frequency",fontsize=16)
    
gcfs_per_scaffold = {}
for i in range(0, len(clusters)):
    gcf1 = list(clusters.keys())[i]
    for scaffold in gcf_scaffolds[gcf1]:
        if scaffold not in gcfs_per_scaffold:
            gcfs_per_scaffold[scaffold] = 0
        gcfs_per_scaffold[scaffold] += 1

number_list = []
average_gcf_counts = 0        
for scaffold in gcfs_per_scaffold:
    average_gcf_counts += gcfs_per_scaffold[scaffold]
    number_list.append(gcfs_per_scaffold[scaffold])
    
print("average gcfs per scaffold: " + str(average_gcf_counts/len(gcfs_per_scaffold)))
outfile.write("average gcfs per scaffold: " + str(average_gcf_counts/len(gcfs_per_scaffold)) + "\n")
axs[1].hist(number_list)
axs[1].set_yscale('log')
axs[1].set_xlabel("Scaffolds per GCF",fontsize=16)
axs[1].set_ylabel("Log Frequency",fontsize=16)
fig.set_constrained_layout(True)
fig.savefig(outfile_histogram)
