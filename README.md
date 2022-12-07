# Genetic-Saturation-Analysis
This is the analysis we used in [Gable and Byars et al. 2022 Systematic Biology](https://academic.oup.com/sysbio/advance-article-abstract/doi/10.1093/sysbio/syac019/6543627?redirectedFrom=fulltext).

This R script will obtain a slope of all raw and TN93-corrected pairwise distances for each alignment, based on Katie Everson's script which analyzed a single gene (https://www.kmeverson.org/blog/simple-dna-saturation-plots-in-r).  The smaller slopes are due to more saturation. The script will analyze the distribution of slopes, determine the mean for you. It will also create a list of alignments that are split by above the mean for less saturated, below the mean for more saturated. The two lists of alignment IDs can be used to make species trees from saturated and unsaturated genes. 

The R script also contains code that can analyze two files containing concatenated (1st and 2nd) and (3rd) codon positions.
