# BGC-clustering-benchmark
Methods for benchmarking BGC similarity comparison and clustering based on structural similarity of produts

# Required dependencies
The following packages need to be installed to run the scripts, version numbers we have used/tested are listed:\n
matplotlib 3.8.3
numpy 1.26.4
pandas 2.1.1
rdkit 2023.9.5
scikit-learn 1.4.1.post1
scipy 1.12.0

We used python version 3.11.8
# 1. Calculating Tanimoto similarity of a database
Note we have precalculated the tanimoto similiarities, if you want to use these precalculated values, proceed to step 2. If you want to do it yourself, you can use the calculate_structural_similarity.py script. To do this:
run:

python calculate_structural_similarity.py

By default, this uses the input file source_data/NPAtlas_bm_v1.tsv and output file tanimoto_results/NPAtlas_bm_v1.tsv, you can use the -f and -o options respectively to change this, e.g.:

python calculate_structural_similarity.py -f source_data/MY_DATA.tsv -o tanimoto_results/MY_DATA.tsv

Note that your input data must follow the same format as data dowloaded from NPAtlas.

# 2. Calculating similarity metrics
To calculate similarity metrics use the bgc_sim_tanimoto_comparison.py script, this has one required option: bgc_file which should contain the similarity scores for pairs of BGCs, it needs to have the format where each line has the similarity score for a pair of BGCs in comma separate form (bgc1,bgc2,score). The bgc_similarities directory has examples for existing methods and can also be used to compare performance of new methods. To run the script run:

python bgc_sim_tanimoto_comparison.py MY_BGC_DATA.csv

optional arguments are:
-t which specifies the tanimoto similarity file, as calculated in step 1
-c which specifies the file which provides the biosynthetic classes of each BGC
-o which specifies the output file name
-g which specifies the output file name to write a graph showing the correlation between the BGC similarity and Tanimoto similarity. 

-t and -c should be used if the user wants to calculate these metrics for their own custom dataset.

# 3. Calculate clusters and scaffolds for product structures
We have provided precalculated clusters using the Butina algorithm and scaffolds using the Bemis-Murcko scaffold definition. If you want to use these precalculated clusters, proceed to 4. If you want to calculate these yourself then use the calc_product_clusters_and_scaffolds.py script by running:

python calc_product_clusters_and_scaffolds.py

If you want to run this on your own dataset, use the -t option to specify Tanimoto similarity file (which can be generated as described in 1) and -n to specify the database in the NPAtlas database format. Other options for this method include -c which specifies the Butina clustering threshold, -b which specifies the output file for the Butina clusters, -s which specifies the output file for scaffolds, and -st which specifies the scaffold type which include Bemis Murcko (true_bm), Cyclic Skeleton (csk), and Bajorath Bemis Murcko (bajorath). Note these scaffold types are described here: https://github.com/rdkit/rdkit/discussions/6844

# 4. Calculate clustering metrics
To calculate clustering metrics use the calculate_clustering_metrics.py script. This script has three required arguments: bgc_cluster_file which is the file containing the BGC cluster identities as determined by the method you would like to calculate metrics for, with the format BGC_id, cluster_id on each line. Examples are provided in the bgc_clusters directory. The second required argument is bgc_score_file which has the BGC similarity scores, as described above. The final required argument is score_type which is used to specify if the score in the score file is a similarity or a distance (sim for similarity dist for distance). To run the script:

python calculate_clustering_metrics.py MY_BGC_CLUSTERS.csv MY_BGC_SCORES.csv SCORE_TYPE

There are also optional arguments if you would like to specify custom ground truth clusters, different scaffold types, or a different product database. These are:
-t specifies product similarity file as calculated in step 1. 
-c specifies the product clsuter file, as calculated in step 3.
-s specifies a file that contains the scaffolds for each BGC, as calculated in step 3. Each line should have BGC_id, scaffold smiles. If there are multiple scaffolds for a BGC due to multiple products, these should each be listed on a separate line
-o specifies the output text file
-p specifies the output file for a histogram showing the scaffold per GCF and GCF per scaffold
