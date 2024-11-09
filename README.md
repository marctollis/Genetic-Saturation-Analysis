# Genetic-Saturation-Analysis
This is the analysis we used in [Gable and Byars et al. 2022 Systematic Biology](https://academic.oup.com/sysbio/advance-article-abstract/doi/10.1093/sysbio/syac019/6543627?redirectedFrom=fulltext).

This R script will obtain the slopes from regression analyses of all raw and TN93-corrected pairwise distances for a large number of DNA alignments, based on [Katie Everson's script](https://www.kmeverson.org/blog/simple-dna-saturation-plots-in-r) which analyzed a single gene, following [Philippe 1994](https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1420-9101.1994.7020247.x) and [Philippe et al. 2011](https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1420-9101.1994.7020247.x). Due to multiple hits, saturated genes will contain fewer phylogenetically useful sites, and this becomes apparent when we compare uncorrected distances between sequences in an alignment to corrected distances that account for things such as proper transition/tranversion ratios, etc. The smaller slopes are due to more saturation (i.e., completely unsaturated genes should have a slope of 1).

There is no objective cutoff that establishes whether a gene is "saturated", but looking at the distribution of these "saturation slopes" in a phylogenomic dataset can help you determine which genes are more versus less saturated. With this information, you can use the slopes in downstream analyses such as potential covariants to determine the effects of saturation levels on phylogenetic inference.

The script will analyze the distribution of slopes, determine the mean for you, perform other descriptive statistics, and plot a histogram. It will also create a list of alignments that are split by above the mean for less saturated, below the mean for more saturated. The two lists of alignment IDs can be used to make species trees from saturated and unsaturated genes, or other types of analysis. 

The R script also contains code that can analyze two files containing concatenated (1st and 2nd) and (3rd) codon positions.

The alignment data is provided in the repository, just download and unzip the files.

