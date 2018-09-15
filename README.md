# Project
This study is to determine how western diet consumption impact the brain of young mice using RNA-sequecing.

## Experimental design: 
Mice of 2 months of age were fed with western diet (WD) or chow diet (CD) for 2,4,6 and 8 weeks. At end of each diet consumption, the cortex region from mouse brain was harvested, snap frozen in dry ice and subjected for RNA-sequencing.  

## Data pre-processing: 
1. fastq files of RNA sequencing was processed using standard JAX in-house pipeline and raw transcript count was generated at the end of the pipeline. 
1. I used standard *edgeR*'s quasi-likelihood pipeline to determine differentially expressed (DE) genes. However, there was few DE genes detected. I have to use an alternative method (GLM mentioned below) to determine DE genes and their corresponding factors. 
1. Based on previous MDS analysis, an outlier of CD for 2 week was removed for the current study. 

## Data Analysis Using Generalized Linear Model (GLM)
1. In this analysis workflow, I used generalized linear model to determine genes affected by multi-factorial effects. To simply to model, I modelled use the duration on diet to predict the gene expression difference between WD and CD. 

## Data Analysis work flow
Raw transcript counts was converted to counts per million (CPM) using *edgeR* package
1. Quality Control 1 - perform filtering, normalization and dispersion estimation on using *edgeR* package
1. Quality Control 2 - perform form correlation and dimension reduction (MDS, PCA) to understand the sample distribution, and what factor contribute to the variance observed.
1. Perform GLM to determine the difference of gene expression between WD and CD affected by during of diet. 
1. Explore the result of GLM 
	+ analyzing coefficients and intercept, and categorize significant genes for each coefficient and intercept.
	+ hierarchical clustering samples based on top genes from GLM analysis

## GO and Kegg pathway analysis
1. Put significant genes from GLM analysis into **Go and Kegg** pathway database to predict the function of these genes or the biological pathways in which the genes are enriched

## Guidance of reviewing analysis
1. Check the code, workflow, and result in **Cortex_outlier_removed.html** or **Cortex_outlier_removed.docx** files
1. Check the full result of **GO and Kegg** in **Analysis** folder
1. Find the home-made functions in the **01_function** folder
1. Find the important analysis figures in **figure** folder
 