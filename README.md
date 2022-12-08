# Genetic-Saturation-Analysis
This is the analysis we used in [Gable and Byars et al. 2022 Systematic Biology](https://academic.oup.com/sysbio/advance-article-abstract/doi/10.1093/sysbio/syac019/6543627?redirectedFrom=fulltext).

This R script will obtain the slopes from regression analyses of all raw and TN93-corrected pairwise distances for a large number of DNA alignments, based on [Katie Everson's script](https://www.kmeverson.org/blog/simple-dna-saturation-plots-in-r) which analyzed a single gene, following [Philippe 1994](https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1420-9101.1994.7020247.x) and [Philippe et al. 2011](https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1420-9101.1994.7020247.x). The smaller slopes are due to more saturation (i.e., completely unsaturated genes should have a slope of 1). 

The script will analyze the distribution of slopes, determine the mean for you, perform other descriptive statistics, and plot a histogram. It will also create a list of alignments that are split by above the mean for less saturated, below the mean for more saturated. The two lists of alignment IDs can be used to make species trees from saturated and unsaturated genes, or other types of analysis. 

The R script also contains code that can analyze two files containing concatenated (1st and 2nd) and (3rd) codon positions.
