# Axonal Localisation Score from 3'end sequencing data

![Overview of the analysis](ga.png)

This repository contains the source code and related to reproduce the figures of the manuscript entitled [The RNA binding proteome of axonal mRNAs in sympathetic neurons](https://www.biorxiv.org/content/10.1101/2022.11.23.517728v1) by R Luisier, C Andreassi, L Fournier and A Riccio.

-   [Dependencies](#Dependencies)
-   [Repo Content](#Repo-Content)
-   [Samples description](#Samples_description)

## Dependencies
### R
Bioconductor version 3.16 (BiocManager 1.30.19), R 4.2.2 (2022-10-31)

The following R packages should be installed:
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


### Python
The python requirements are listed in `python_requirements.txt`. 
To run the jupyter notebook, create a python environment with [conda](https://docs.conda.io/en/latest/) by running the command `conda create --name <env_name> --file python_requirements.txt` (tested on osx-64). You can then activate the environment by running `conda activate <env_name>`.


## Repo Content
* [data](./data): folder containing the data for examples matrix of gene expression; etc. Raw sequencing data will be deposited publicly.
* [Scripts](./scripts): `R`, `Python` and `Bash` custome code

## Samples description
3’ end sequencing of RNA isolated from axons and cell bodies of sympathetic neurons exposed to either Nerve Growth factor (NGF) or Neurotrophin 3 (NT3). 
Source code related to the manuscript **The RNA Binding proteome of axonal mRNAs in sympathetic neurons**, Luisier et al. (2023).

## Summary of the analysis and results

### 1. Overview of the data
The analysis of the read count and samples is presented [HERE](https://htmlpreview.github.io/?https://github.com/RLuisier/AxonLoc/blob/main/1_overview_data.html) that relates to Supplementary Figure 1.

### 2. Differential gene expression analysis
The analysis of the differential gene expression analysis between NT3 and NGF conditions together with transcription factor analysis is presented in [DGE analysis](https://htmlpreview.github.io/?https://github.com/RLuisier/AxonLoc/blob/main/DGE_cell_body.html). This relates to first half of Figure 1.


### 3. Differential Alternative Polyadenyation between NGF and NT3
The analysis of differential APA between NGF and NT3 and related results are presented in  [APA analysis](https://htmlpreview.github.io/?https://github.com/RLuisier/AxonLoc/blob/main/Differential_APA_CB.html), that relates to the second half of Figure 1.

### 4. Differential Alternative Polyadenyation between NGF and NT3
The RBPome study that underlies APA in developing sympathetic neurons was performed by integrating the 3'end sequencing data obtained from sympathetic neurons with publicly available CLIP-sequencing data obtained from human cell lines and lifted over the rat genome. The analysis and results are presented in [APA RBPome](https://htmlpreview.github.io/?https://github.com/RLuisier/AxonLoc/blob/main/Regulation_APA.html) and relates to Figure 2.

### 5. Axonal versus cell body compartments
The analysis of the compartment-specific mRNA pools is presented in  [compartment analysis](https://htmlpreview.github.io/?https://github.com/RLuisier/AxonLoc/blob/main/Analysis_compartment.html), that relates to the first part of the Figure 3.

### 6. Axonal Localisation Score
The modeling of the axonal localisation alongside analysis of differential localisation between NGF and NT3 condition are presented in the  [LS analysis analysis](https://htmlpreview.github.io/?https://github.com/RLuisier/AxonLoc/blob/main/4_Axonal_Localisation_Scoring.html), that relates to the second part of the Figure 3.

### 7. Predicted RBPome that shapes the axonal transcripome

### 8. Synergistic regulatory potential of RBPs in axonal localisation
The analysis of the synergistic regulatory potential of RBPs in axonal localisation in presented [synergistic analysis](https://htmlpreview.github.io/?https://github.com/RLuisier/AxonLoc/blob/main/6_RBP_regulome_localisation.html), that relates to the second part of Figure 5.












