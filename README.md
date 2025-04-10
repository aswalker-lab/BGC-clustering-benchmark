# BGC-clustering-benchmark
Methods for benchmarking BGC similarity comparison and clustering based on structural similarity of produts

#1. Calculating Tanimoto similarity of a database
Note we have precalculated the tanimoto similiarities, but if you want to do it yourself, you can use the calculate_structural_similarity.py script. To do this:
run:

python calculate_structural_similarity.py

By default, this uses the input file source_data/NPAtlas_bm_v1.tsv and output file tanimoto_results/NPAtlas_bm_v1.tsv, you can use the -f and -o options respectively to change this, e.g.:

python calculate_structural_similarity.py -f source_data/MY_DATA.tsv -o tanimoto_results/MY_DATA.tsv

Note that your input data must follow the same format as data dowloaded from NPAtlas.
