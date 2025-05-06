# Regenie.QRS:computationally efficient whole-genome quantile regression 
This repository contains code and instructions for performing quantile regression (QR) genome-wide association studies (GWAS) with related samples for quantitative traits, as described in our paper ["Regenie.QRS: computationally efficient whole-genome quantile regression at biobank scale"].

## Install dependent sofware and R packages

Regenie [https://rgcgithub.github.io/regenie/](https://rgcgithub.github.io/regenie/)

quantreg [https://cran.r-project.org/package=quantreg](https://cran.r-project.org/package=quantreg)

data.table [https://CRAN.R-project.org/package=data.table](https://CRAN.R-project.org/package=data.table)

dplyr [https://CRAN.R-project.org/package=dplyr](https://CRAN.R-project.org/package=dplyr)


## Example

We provide a toy dataset in [example](/example). 
The input data includes 
1) the binary genotype file in [example](/example) contains N sample = 500 individuals (column "IID") 
and N variant = 100 SNPs from the [https://github.com/cloufield/GWASTutorial.git/01_Dataset/download_sampledata.sh](https://github.com/cloufield/GWASTutorial.git/01_Dataset/download_sampledata.sh).
2) the [phenotype file](example/normalized_pheno.txt) that contains one simulated quantitative trait (column "Y1") for the N sample = 500 individuals (column "IID").
3) the [covariate file](example/covariate.txt) that contains two simulated covariates (columns "covar1" and "covar2") for the N sample = 500 individuals (column "IID").

Run the R script [example.regenie.qrs.R](example.regenie.qrs.R) to perform Regenie.QRS (single variant tests) for the phenotype and genotype data in [example](/example). Set ```is.effect.estimated = T``` to enable the estimation of quantile-specific effect size (disabled by default). Set ```data.path``` to the directory where you store your data and results, and set ```Regenie.path``` to the location where the Regenie executable is installed.


The [expected output](example/example.sumstat.tsv) is a tab-delimited text file with the QR summary statistics. 
- Column "ID": variant ID from the genotype data.
- Column "P_QR": the integrated Regenie.QRS p-value across multiple quantile levels.
- Columns from "P_QR0.1" to "P_QR0.9": quantile-specific p-value for the quantile levels 0.1, 0.2, ..., 0.9.
- Columns from "BETA_QR0.1" to "BETA_QR0.9": quantile-specific effect size for the quantile levels 0.1, 0.2, ..., 0.9.



## Suggestions for genome-wide analysis

We recommend splitting genome-wide genotype data into smaller variant subsets or small genomic regions to improve computational efficiency. Data scattering reduces the computational burden on a single machine and enables parallel execution across multiple machines, significantly enhancing scalability. 

For example, in a dataset with 382,402 individuals and 6,364,586 variants, performing Regenie.QRS GWAS for nine quantile levels took 660 CPU hours using R. To optimize this process, we split the genome into 1,073 segments, each containing ~6000 variants. Each segment required 36 minutes to complete Regenie.QRS testing (without effect size estimation) when executed on a single CPU core. By leveraging data scatter and parallel computing, we distributed the workload across 1,073 CPU cores. This reduced the wall-clock time to 36 minutes for 9 quantile levels, making the QR GWAS feasible for biobank-scale data.


Users may consider this data partitioning strategry and determine the number of genomic chunks based on data volume and available computing resources. The R code for QR GWAS requires loading a genotype matrix in memory. Loading whole-genome or chromosome-wide genotype data into R consumes significant memory and processing time. It would be much easier to load a smaller chunk of genotype data into memory and then perform QR on the small subsets. 

Users can directly import the commonly used PLINK/VCF format data into R using packages such as [genio](https://github.com/OchoaLab/genio) and [seqminer](https://github.com/zhanxw/seqminer). 


## Citation

Fan Wang, Chen Wang, Tianying Wang, Marco Masala, Edoardo Fiorillo, Marcella Devoto, Francesco Cucca, Iuliana Ionita-Laza, "Regenie.QRS: computationally efficient whole-genome quantile regression at biobank scale", biorxiv, 2025
