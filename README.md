# Optimal design for cell-type specific eQTL studies

<p align="center">

<img src="./figures/figure_s14.png" width="50%">

</p>

This repository contains the code used to analyse the single-cell
RNA-seq datasets shown in:

*Mandric, I., Schwarz, T., Majumdar, A., Hou, k., Bricscoe, L., Perez, R.,
Subramaniam, M., Hafemeister, C., Satija, R., Ye, C., Pasaniuc, B., Halperin, E. (2019) Optimal design of 
single-cell RNA sequencing experiments for cell-type-specific eQTL analysis.
<https://doi.org/10.1101/766972>*

  - [Data availability](#data-availability)
  - [Figure shortcuts](#figure-shortcuts)
  - [Analysis preliminaries](#analysis-preliminaries)
  - [1. Load and hygienize dataset](#1-load-and-hygienize-dataset)
  - [2. Knowledge-based identification of all cell
    populations](#2-knowledge-based-identification-of-all-cell-populations)

# 

## Data availability

# 

Those interested in processing the datasets independently should consider
downloading:

  - the [Smart-Seq2 dataset](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/)
    from Segerstolpe et al., 2016.  
  - the [Census of Immune Cells](https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79). We analyzed only a subset of it (Lane 7).
  - the [10X dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137029)
