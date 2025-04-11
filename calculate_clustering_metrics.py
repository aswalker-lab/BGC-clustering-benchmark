# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 10:04:28 2025

@author: Allison Walker
This script calculates supervised and unsupervised clustering metrics
"""

import argparse

parser = argparse.ArgumentParser(description='Compares BGC similarity to Tanimoto similarity')
parser.add_argument('bgc_cluster_file', type=str, help='file containing BGC cluster labels, must have a format where BGC id is in the first column and cluster id is in the second column')
parser.add_argument('bgc_score_file', type=str, help='file containing BGC similarity, must have a format where columns one and two are BGC ids and column 3 is a similarity or distance measurement')
parser.add_argument('score_type', type=str, help='Score type, options are dist for distance and sim for similarity',choices=["dist","sim"])
parser.add_argument('-t', '--tanimoto_file', type=str, default='tanimoto_results/NPAtlas_bm_v1.tsv', help='file containing tanimoto comparison of product structures')
parser.add_argument('-c', '--cluster_file', type=str, default='tanimoto_results/NPAtlas_bm_v1.tsv', help='ground truth cluster file based on structure or user-defined clusters')
parser.add_argument('-o', '--outfile_name', type=str, default='tanimoto_bgc_comparison.txt', help='file to write text results to')

args = parser.parse_args()
infile_name = args.bgc_cluster_file
score_file = args.bgc_score_file
score_type = args.score_type
tanimoto_file = args.tanimoto_file
cluster_file = args.cluster_file
outfile_name = args.outfile_name



