# BGC-clustering-benchmark
Methods for benchmarking BGC similarity comparison and clustering based on structural similarity of produts

# 1. Calculating Tanimoto similarity of a database
Note we have precalculated the tanimoto similiarities, but if you want to do it yourself, you can use the calculate_structural_similarity.py script. To do this:
run:

python calculate_structural_similarity.py

By default, this uses the input file source_data/NPAtlas_bm_v1.tsv and output file tanimoto_results/NPAtlas_bm_v1.tsv, you can use the -f and -o options respectively to change this, e.g.:

python calculate_structural_similarity.py -f source_data/MY_DATA.tsv -o tanimoto_results/MY_DATA.tsv

Note that your input data must follow the same format as data dowloaded from NPAtlas.

# 2. Calculating similarity metrics
To calculate similarity metrics use the bgc_sim_tanimoto_comparison.py script, this has one required option: bgc_file which should contain the similarity scores for pairs of BGCs, it needs to have the format where each line has the similarity score for a pair of BGCs in comma separate form (bgc1,bgc2,score). The bgc_similarities directory has examples for existing methods and can also be used to compare performance of new methods. To run the script run:

python bgc_sim_tanimoto_comparison.py MY_BGC_DATA.csv

optional arguments are -t which specifies the tanimoto similarity file, -c which specifies the file which provides the biosynthetic classes of each BGC, -o which specifies the output file name, and -g which specifies the output file name to write a graph showing the correlation between the BGC similarity and Tanimoto similarity.
