# Axonal Localisation Score from 3'end Seq
RNA-seq data analysis from longitudinal study of patient-specific iPSC-derived cultures of differentiating motor neurons with nuclear and cytoplasmic fractionation.

-   [Dependencies](#Dependencies)
-   [Repo Content](#Repo-Content)
-   [Samples description](#Samples_description)
-   [Overview of the analysis](#Overview)
    1.   [*Import annotation file*](./0-source-annotation-files.md)
    1.   [*Preprocessing, Mapping and QC*](./1-preprocessing-mapping-qc.md)
    1.   [*Quality control analysis*](./2-QC_analysis.md)
    1.   [*Characterisation of the samples in terms of gene expression*](./2-QC_analysis.md)
    1.   [*Splicing analysis*](./4-splicing-analysis.md)
    1.   [*Intron retention focussed analysis*](./5-IR-focus-analysis.md)

## Dependencies
STAR-2.6.0a, fastqc-0.11.5, fastx_toolkit-0.0.14, samtools-1.9, bowtie-1.1.2, RSeqQC, python-3.7.3, kallisto 0.46.0, R-3.6.0. The following packages should be installed:


```R
install.packages(c('grDevices','Rsamtools','IRanges','GenomicRanges','rtracklayer','GenomicFeatures','GenomicAlignments','Segmentor3IsBack','GO.db','limma','topGO','biomaRt,'geneplotter','multtest','mclust'))
```

Finally add this to your  ~/.bashrc
```bash
export PATH=/home/r/raphaelle-luisier/Scripts/AxonLoc/scripts:$PATH
```

## Repo Content
* [annotation](./annotation): folder containing all relevant annotation files including samples annotation and qc report.
* [data](./data): folder containing the data for examples matrix of gene expression; etc. Raw sequencing data will be deposited publicly.
* [Scripts](./scripts): `R`, `Python` and `Bash` custome code

## Samples description
3'end Seq was obtained from rat sympathetic neurons...










