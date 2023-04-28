# Axonal Localisation Score from 3'end sequencing data


-   [Dependencies](#Dependencies)
-   [Repo Content](#Repo-Content)
-   [Samples description](#Samples_description)

## Dependencies
The following packages should be installed:
GenomicRanges_1.50.2
Rsamtools_2.14.0
rtracklayer_1.58.0
IRanges_2.32.0
geneplotter_1.76.0
multtest_2.54.0
mclust_6.0.0
knitr_1.42
edgeR_3.40.2
topGO_2.50.0         
SparseM_1.81         
graph_1.76.0         
plotly_4.10.1
fitdistrplus_1.1-8
GO.db_3.16.0 

Bioconductor version 3.16 (BiocManager 1.30.19), R 4.2.2 (2022-10-31)


## Repo Content
* [data](./data): folder containing the data for examples matrix of gene expression; etc. Raw sequencing data will be deposited publicly.
* [Scripts](./scripts): `R`, `Python` and `Bash` custome code

## Samples description
3’ end sequencing of RNA isolated from axons and cell bodies of sympathetic neurons exposed to either Nerve Growth factor (NGF) or Neurotrophin 3 (NT3). 
Source code related to the manuscript **The RNA Binding proteome of axonal mRNAs in sympathetic neurons**, Luisier et al. (2023).

## Overview of the data
The analysis of the read count and samples is presented [Here](https://htmlpreview.github.io/?https://github.com/RLuisier/AxonLoc/blob/main/1_overview_data.nb.html) which source code can be acessed in [./1_overview_data.nb.Rmd].












