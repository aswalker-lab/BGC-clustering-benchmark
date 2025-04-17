# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 16:18:48 2025

@author: Allison Walker
This script calculates the structural similarity based on NPAtlas input data
"""

from rdkit import Chem
from rdkit.Chem import AllChem 
from rdkit import DataStructs
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Calculates tanimoto similarity between NPs in dataset')
parser.add_argument('-f', '--filename', type=str, default='source_data/NPAtlas_bm_v1.tsv', help='file containing NP structural information in the NPAtlas format')
parser.add_argument('-o', '--outfile_name', type=str, default='tanimoto_results/NPAtlas_bm_v1.tsv', help='file to write tanimoto results to')
args = parser.parse_args()
infile_name = args.filename
outfile_name = args.outfile_name

#process NPAtlas data
NP_atlas_tsv = open(infile_name, encoding="utf8")
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


#calculate fingerprints
molecules = {}
fingerprints = {}
avg_fingerprints = {}
fingerprint_counts = {}

fpgen = AllChem.GetRDKitFPGenerator(fpSize=8192)
for bgc in smiles_dic:
    if bgc not in molecules:
        molecules[bgc] = []
        fingerprints[bgc] = []
        fingerprint_counts[bgc] =0
        
    for smiles in smiles_dic[bgc]:
        fingerprint_counts[bgc] += 1
        mol = Chem.MolFromSmiles(smiles)
        molecules[bgc].append(mol)
        fingerprints[bgc].append(fpgen.GetFingerprint(mol))
        if bgc not in avg_fingerprints:
            avg_fingerprints[bgc] = np.array(fpgen.GetFingerprint(mol))
        else:
            avg_fingerprints[bgc] += np.array(fpgen.GetFingerprint(mol))

#calculate tanimoto similarities
tanimoto = {}
for bgc1 in fingerprints:
    tanimoto[bgc1] = {}
    for bgc2 in fingerprints:
        max_tanimoto = 0
        for fp1 in fingerprints[bgc1]:
            for fp2 in fingerprints[bgc2]:
                    score = DataStructs.TanimotoSimilarity(fp1, fp2)
                    if score > max_tanimoto:
                        max_tanimoto = score
        tanimoto[bgc1][bgc2] = max_tanimoto

#write output data
outfile = open(outfile_name,'w')
outfile.write("MiBIG_ID")
for bgc1 in tanimoto:
    outfile.write(","+bgc1)
outfile.write("\n")
for bgc1 in tanimoto:
    count = 0
    outfile.write(bgc1 + ",")
    for bgc2 in tanimoto:
        if count >= len(tanimoto[bgc1]) -1:
            outfile.write(str(tanimoto[bgc1][bgc2])+"\n")
        else:
            outfile.write(str(tanimoto[bgc1][bgc2]) + ",")
        count +=1 
outfile.close()

