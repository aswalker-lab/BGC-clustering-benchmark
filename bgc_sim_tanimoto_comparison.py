# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 11:32:34 2025

@author: Allison Walker
This script compares BGC similarity to Tanimoto similarity
Biosynthetic classes of bgcs can be definied using the -c option to provide a file with 
classifications. This allows the user to define different classes than the ones used in our paper
"""

import pandas as pd
import matplotlib.pyplot as plt
import scipy
import argparse

parser = argparse.ArgumentParser(description='Compares BGC similarity to Tanimoto similarity')
parser.add_argument('bgc_file', type=str, help='file containing BGC similarity, must have a format where columns one and two are BGC ids and column 3 is a similarity or distance measurement')
parser.add_argument('-t', '--tanimoto_file', type=str, default='tanimoto_results/NPAtlas_bm_v1.tsv', help='file containing tanimoto comparison of product structures')
parser.add_argument('-c', '--class_file', type=str, default='bgc_class.csv', help='file containing BGC biosynthetic classes')
parser.add_argument('-o', '--outfile_name', type=str, default='tanimoto_bgc_comparison.txt', help='file to write text results to')
parser.add_argument('-g', '--graphfile_name', type=str, default='tanimoto_bgc_comparison.png', help='file to write graph of results to')

args = parser.parse_args()
infile_name = args.bgc_file
tanimoto_file = args.tanimoto_file
outfile_name = args.outfile_name
graphfile_name = args.graphfile_name
class_file = args.class_file

#load tanimoto data
tanimoto_matrix = pd.read_csv(tanimoto_file)
tanimoto_matrix.set_index('MiBIG_ID', inplace=True)

#load BGC classes
bgc_class_in = open(class_file, encoding='cp1252')
bgc_classes = {}
for line in bgc_class_in:
    bgc_id = line.split(",")[0]
    class_list = line.split(",")[3].replace("\n","").split(":")    
    bgc_classes[bgc_id] = class_list

#load BGC data
bgc_scores = {}
for line in open(infile_name):
    split_line = line.split(",")
    bgc1 = split_line[0]
    bgc2 = split_line[1]
    score = float(split_line[2])
    if bgc1 not in bgc_scores:
        bgc_scores[bgc1] = {}
    bgc_scores[bgc1][bgc2] = score
    
#order BGC scores and tanimoto scores
#if BGC data is missing tanimoto scores it will be excluded    
#if BGC data is missing classificaitons it will be excluded from within class correlations
bgc_score_list = []
tanimoto_list = []
hybrid_bgc_score = {}
no_hybrid_bgc_score = {}
hybrid_tanimoto_score = {}
no_hybrid_tanimoto_score = {}
for bgc1 in bgc_scores:
    if bgc1 not in tanimoto_matrix:
        continue
    for bgc2 in bgc_scores[bgc1]:
        if bgc2 not in tanimoto_matrix:
            continue
        bgc_score_list.append(bgc_scores[bgc1][bgc2])
        tanimoto_list.append(tanimoto_matrix[bgc1][bgc2])
        if bgc1 not in bgc_classes or bgc2 not in bgc_classes:
            continue
        intersection_classes = list(set(bgc_classes[bgc1]) & set(bgc_classes[bgc2]))
        #if there are no shared classes between the BGCs do not add to class specific lists
        if len(intersection_classes) < 1:
            continue
        #add to all shared class lists
        for c in intersection_classes:
            if c not in hybrid_bgc_score:
                hybrid_bgc_score[c] = []
                hybrid_tanimoto_score[c] = []
            hybrid_bgc_score[c].append(bgc_scores[bgc1][bgc2])
            hybrid_tanimoto_score[c].append(tanimoto_matrix[bgc1][bgc2])
        #only add to no hybrid list if both BGCs are not hybrids
        if len(bgc_classes[bgc1]) > 1 or len(bgc_classes[bgc2]) > 1:
            continue
        c = bgc_classes[bgc1][0]
        if c not in no_hybrid_bgc_score:
            no_hybrid_bgc_score[c] = []
            no_hybrid_tanimoto_score[c] = []
        no_hybrid_bgc_score[c].append(bgc_scores[bgc1][bgc2])
        no_hybrid_tanimoto_score[c].append(tanimoto_matrix[bgc1][bgc2])
        
#make scatter plot for all BGCs
plt.scatter(bgc_score_list, tanimoto_list)
plt.ylabel('Tanimoto Score', fontsize=18)
plt.xlabel('BGC Score', fontsize=18)
plt.ylim((0,1))
plt.xlim((0,1))
plt.savefig(graphfile_name)

outfile = open(outfile_name, 'w')

#calculate correlations for all BGCs and individual BGCs
(spearman_r, spearman_p) = scipy.stats.spearmanr(bgc_score_list, tanimoto_list) 

print("Spearman correlation for all BGCs: " + str(spearman_r) + " p-value: " + str(spearman_p)) 
outfile.write("Spearman correlation for all BGCs: " + str(spearman_r) + " p-value: " + str(spearman_p) + "\n")

print("Correlation for BGCs by class including hybrids")
outfile.write("Correlation for BGCs by class including hybrids" + "\n")
for c in hybrid_bgc_score:
    (spearman_r, spearman_p) = scipy.stats.spearmanr(hybrid_bgc_score[c], hybrid_tanimoto_score[c]) 
    print("Spearman correlation for " + c + " BGCs: " + str(spearman_r) + " p-value: " + str(spearman_p))
    outfile.write("Spearman correlation for " + c + " BGCs: " + str(spearman_r) + " p-value: " + str(spearman_p) + "\n")
    
print("Correlation for BGCs by class excluding hybrids")
outfile.write("Correlation for BGCs by class excluding hybrids")
for c in no_hybrid_bgc_score:
    (spearman_r, spearman_p) = scipy.stats.spearmanr(no_hybrid_bgc_score[c], no_hybrid_tanimoto_score[c]) 
    print("Spearman correlation for " + c + " BGCs: " + str(spearman_r) + " p-value: " + str(spearman_p))
    outfile.write("Spearman correlation for " + c + " BGCs: " + str(spearman_r) + " p-value: " + str(spearman_p) + "\n")
outfile.close()