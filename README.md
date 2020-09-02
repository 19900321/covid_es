# covid_es

This is script for calculating erichment score between one drug and some disease related gene sets.

input: (1).Gene expression value of perturbagen compared with control sample without drug. (cmap here)

(2). disease gene set of upregulated one and downregulated one.

Main script:

1: ges.py for prepared the cmap data to be a matrix (drug as row and gene as columns)

2: calculate.py is the method in https://stm.sciencemag.org/content/3/96/96ra77/tab-pdf
https://www.nature.com/articles/s41421-020-0153-3 as well as generate random disease gene sets

3: cmap.py is another ES methods. https://clue.io/connectopedia/cmap_algorithms
