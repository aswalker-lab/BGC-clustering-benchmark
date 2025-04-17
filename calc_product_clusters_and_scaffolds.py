# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 13:35:43 2025

@author: Allison Walker
Butina clustering adapted from https://projects.volkamerlab.org/teachopencadd/talktorials/T005_compound_clustering.html
Scaffold determination adapted from https://github.com/rdkit/rdkit/discussions/6844
"""

import pandas as pd
import numpy
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import AllChem 
import argparse
from sklearn import metrics
from rdkit.Chem.Scaffolds import MurckoScaffold

parser = argparse.ArgumentParser(description='Compares BGC similarity to Tanimoto similarity')
parser.add_argument('-t', '--tanimoto_file', type=str, default='tanimoto_results/NPAtlas_bm_v1.tsv', help='file containing tanimoto comparison of product structures')
parser.add_argument('-n', '--npatlas_filename', type=str, default='source_data/NPAtlas_bm_v1.tsv', help='file containing NP structural information in the NPAtlas format')
parser.add_argument('-c', '--butina_threshold', type=float, default=0.2, help='Butina cluster threshold')
parser.add_argument('-b', '--butina_outfile_name', type=str, default='product_clusters/butina_clusters_pt2.csv', help='file to write butina clusters to')
parser.add_argument('-s', '--scaffolds_outfile_name', type=str, default='product_scaffolds/bm_scaffolds.csv', help='file to write scaffolds to')
parser.add_argument('-st', '--scaffold_type', type=str, default='real_bm', choices=['real_bm','csk','bajorath'], help='type of scaffold to calculate')

args = parser.parse_args()
tanimoto_file = args.tanimoto_file
npatlas_filename = args.npatlas_filename
butina_outfile_name = args.butina_outfile_name
butina_threshold = args.butina_threshold
scaffolds_outfile_name = args.scaffolds_outfile_name
scaffold_type = args.scaffold_type

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

def get_scaffold(mol,scaffold_type="real_bm"):
    """Get scaffolds for a molecule

    Parameters:
    mol : rdkit molecule.
    scaffold type : String indicating scaffold type
        options are real_bm, bajorath, csk
        The default is real_bm
    
    Returns the scaffold

    """
    PATT=Chem.MolFromSmarts("[$([D1]=[*])]")
    REPL=Chem.MolFromSmarts("[*]")
    Chem.RemoveStereochemistry(mol) #important for canonization of CSK!
    scaff=MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold_type == "bajorath":
        scaff=AllChem.DeleteSubstructs(scaff, PATT)
    if scaffold_type == "real_bm":
        scaff=AllChem.ReplaceSubstructs(scaff,PATT,REPL,replaceAll=True)[0]                                      
    if scaffold_type == "csk":
        scaff=MurckoScaffold.GetScaffoldForMol(scaff)
    return scaff

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

#load structures for scaffold calculation
NP_atlas_tsv = open(npatlas_filename, encoding="utf8")

smiles_dic = {}
for line in NP_atlas_tsv:
    if "npaid" in line:
        continue
    split_line = line.split("\t")
    if len(split_line) < 29:
        continue
    smiles =split_line[10]
    mibig = split_line[28]
    mibig = mibig.replace("'","")
    mibig = mibig.replace("[","").replace("]","").replace(" ","")
    if len(mibig) < 1:
        mibig_list = []
    else:
        mibig_list = mibig.split(",")

    for bgc in mibig_list:
        if bgc not in smiles_dic:
            smiles_dic[bgc] = []
        smiles_dic[bgc].append(smiles)
        
molecules = {}
scaffolds = {}

bgc_list = []
outfile = open(scaffolds_outfile_name,'w')
for bgc in smiles_dic:
    smiles = smiles_dic[bgc][0]

        
    bgc_list.append(bgc.replace("\n",""))
    if bgc not in molecules:
        molecules[bgc] = []
        scaffolds[bgc] = []
        
        
    for smiles in smiles_dic[bgc]:
        mol = Chem.MolFromSmiles(smiles)
        Chem.SanitizeMol(mol)
        molecules[bgc].append(mol)
        if scaffold_type == "real_bm":
            scaffold=get_scaffold(mol,scaffold_type ="real_bm")
        elif scaffold_type == "csk":
            scaffold=get_scaffold(mol,scaffold_type = "csk")
        else:
            scaffold=get_scaffold(mol,scaffold_type = "bajorath")
        outfile.write(bgc + "," + Chem.MolToSmiles(scaffold) + "\n")  
outfile.close()