# Optimal design for cell-type specific eQTL studies

<p align="center">

<img src="./figures/Figure-S15.png" width="50%">

</p>

This repository contains the code used to analyse the single-cell
RNA-seq datasets shown in:

*Mandric, I., Schwarz, T., Majumdar, A., Hou, k., Bricscoe, L., Perez, R.,
Subramaniam, M., Hafemeister, C., Satija, R., Ye, C., Pasaniuc, B., Halperin, E. (2019) Optimal design of 
single-cell RNA sequencing experiments for cell-type-specific eQTL analysis.
<https://doi.org/10.1101/766972>*

  - [Data availability](#data-availability)
  - [Figures](#figures)
  - [Figure captions](#figure-captions)
  - [ct-eQTL design calculator](#ct-eqtl-design-calculator)

## Data availability


If you want to analyze the datasets used in the paper, you can access them independently:

  - the [Smart-Seq2 dataset](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/)
    from Segerstolpe et al., 2016.  
  - the [Census of Immune Cells](https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79). We analyzed only a subset of it (Lane 7).
  - the [10X dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137029)


## Figures



|                               Figure 1                               |                              Figure 2                              |                     Figure 3                     |                Figure 4                 |
| :------------------------------------------------------------------: | :----------------------------------------------------------------: | :----------------------------------------------: | :-------------------------------------: |
|                    ![](./figures/Figure-1.png)                    |                   ![](./figures/Figure-2.png)                   |          ![](./figures/Figure-3.png)          |     ![](./figures/Figure-4.png)      |
| [Caption](#caption-figure-1), [Script](./plots/figure1.R)  | [Caption](#caption-figure-2), [Script](./plots/figure2.R) | [Caption](#caption-figure-3), [Script](./plots/figure3.R) | [Caption](#caption-figure-4) |

|            Figure 5            |            Figure S1            |                                           Figure S2                                            |                                             Figure S3                                             |
| :----------------------------: | :----------------------------: | :-------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------: |
| ![](./figures/Figure-5.png) | ![](./figures/Figure-S1.png) |                                ![](./figures/Figure-S2.png)                                 |                                  ![](./figures/Figure-S3.png)                                  |
| [Caption](#caption-figure-5)  | [Caption](#caption-figure-s1) | [Caption](#caption-figure-s2), [Script](./plots/allplots.R) | [Caption](#caption-figure-s3) |

|               Figure S4                |             Figure S5              |                 Figure S6                 |                                          Figure S7                                           |
| :------------------------------------: | :--------------------------------: | :---------------------------------------: | :------------------------------------------------------------------------------------------: |
|    ![](./figures/Figure-S4.jpg)     |  ![](./figures/Figure-S5.jpg)   |      ![](./figures/Figure-S6.jpg)      |                               ![](./figures/Figure-S7.jpg)                                |
| [Caption](#caption-figure-s4)  | [Caption](#caption-figure-s5) | [Caption](#caption-figure-s6) | [Caption](#caption-figure-s7) |

|               Figure S8                |             Figure S9              |                 Figure S10                 |                                          Figure S11                                           |
| :------------------------------------: | :--------------------------------: | :---------------------------------------: | :------------------------------------------------------------------------------------------: |
|    ![](./figures/Figure-S8.png)     |  ![](./figures/Figure-S9.png)   |      ![](./figures/Figure-S10.png)      |                               ![](./figures/Figure-S11.png)                                |
| [Caption](#caption-figure-s8)  | [Caption](#caption-figure-s9) | [Caption](#caption-figure-s10) | [Caption](#caption-figure-s11) |


|               Figure S12                |             Figure S13              |                 Figure S14                 |                                          Figure S15                                           |
| :------------------------------------: | :--------------------------------: | :---------------------------------------: | :------------------------------------------------------------------------------------------: |
|    ![](./figures/Figure-S12.png)     |  ![](./figures/Figure-S13.png)   |      ![](./figures/Figure-S14.png)      |                               ![](./figures/Figure-S15.png)                                |
| [Caption](#caption-figure-s12)  | [Caption](#caption-figure-s13) | [Caption](#caption-figure-s14) | [Caption](#caption-figure-s15) |


|               Figure S16                |             Figure S17              |                 Figure S18                 |                                          Figure S19                                           |
| :------------------------------------: | :--------------------------------: | :---------------------------------------: | :------------------------------------------------------------------------------------------: |
|    ![](./figures/Figure-S16.png)     |  ![](./figures/Figure-S17.png)   |      ![](./figures/Figure-S18.png)      |                               ![](./figures/Figure-S19.png)                                |
| [Caption](#caption-figure-s16)  | [Caption](#caption-figure-s17) | [Caption](#caption-figure-s18) | [Caption](#caption-figure-s19) |



## Figure captions

| <center> <h5 id="caption-figure-1">Caption Figure 1</h5> </center> |
| :------------------------------------: |
| The impact of read coverage on the average R2 between cell-type-specific gene expression estimates and their high-coverage values (Smart-Seq2 dataset, alpha cells).  A) Distribution of Pearson R2 computed across all the genes at different levels of read coverage, Smart-Seq2 dataset. B) Distribution of Pearson R2 at 75,000 reads per cell stratified by the expression level, Smart-Seq2 dataset. C) Distribution of Pearson R2 computed across all the genes at different levels of read coverage, 10X dataset. D) Distribution of Pearson R2 at 4,000 reads per cell stratified by the expression level, 10X dataset. |
| <center> <h5 id="caption-figure-2">Caption Figure 2</h5> </center> |
| Effective sample size computed across a grid of experimental designs with sample size N ranging from 40 to 120 individuals in steps of 8 and the number of cells per individuals M ranging from 500 to 2,750 cells per individual in steps of 250 (CD4 T cells).  A) Library preparation is assumed to be 0$ per reaction, level of multiplexing is fixed and equal to 8. B) Library preparation is set to $2000 per reaction, level of multiplexing is fixed and equal to 8. C) Library preparation is set to $2000 per reaction, greedy multiplexing. D) Library preparation is set to $2000 per reaction, greedy multiplexing, demultiplexing inaccuracy and cell-type misclassification is taken into account. |
| <center> <h5 id="caption-figure-3">Caption Figure 3</h5> </center> |
| Cell-type misclassification error rate (%) and coverage (thousands of reads per cell) for different experimental designs.  A) Cell-type misclassification error; B) Coverage. |
| <center> <h5 id="caption-figure-4">Caption Figure 4</h5> </center> |
| Experimental designs for CD4 T cells ct-eQTL with effective sample size Neff=40: A) Comparison of different experimental designs. Experimental design N=88, M=2,250, r=4,500 yields two-fold reduction in cost than the standard design. B) For a fixed sample size and number of cells per individual, increasing coverage implies increasing the effective sample size (i.e., power) only up to a point. There is little gain in power at coverages greater than 12,500 reads per cell. |
| <center> <h5 id="caption-figure-5">Caption Figure 5</h5> </center> |
| Effective sample size as a function of cell-type prevalence. Shown here is the effective sample size across the grid of experimental design when the cell-type abundance is set to different values - 5, 10, 15, 20, 25, 30%.  (CD4 T cells at fixed budget $35,000). |
| <center> <h5 id="caption-figure-6">Caption Figure 6</h5> </center> |
| Performance of ct-eQTL analysis. Shown here is recall (power estimate) as a function of coverage in the ct-eQTL analysis of CD4 T cells at fixed budget $35,000. A) Mean ct-eQTL; B) Variance ct-eQTL. |
| <center> <h5 id="caption-figure-s1">Caption Figure S1</h5> </center> |
| Pearson R2 between low-coverage estimates and the high-coverage gene expression in Smart-Seq2 dataset. A) SLC14A2 gene, 50% downsampling (375,000 reads per cell); B) SLC14A2 gene, 10% downsampling (75,000 reads per cell); C) GCG gene, 50% downsampling (375,000 reads per cell); D) GCG gene, 10% downsampling (75,000 reads per cell). |
| <center> <h5 id="caption-figure-s2">Caption Figure S2</h5> </center> |
| Stratification of genes based on the number of individuals they are expressed in (Smart-Seq2 dataset). A) Distribution of genes by the number of individuals they are expressed in, Smart-Seq2 dataset. B) Average Pearson R2 at 75,000 reads per cell (alpha cells) stratified by the number of individuals they are expressed in, Smart-Seq2 dataset (vertical bars indicate interquartile range). C) Distribution of genes by the number of individuals they are expressed in, 10X dataset (The Census of Immune Cells). D) Average Pearson R2 at 4,000 reads per cell (erythroblast cells) stratified by the number of individuals they are expressed in, 10X dataset (The Census of Immune Cells) (vertical bars indicate interquartile range). |
| <center> <h5 id="caption-figure-s3">Caption Figure S3</h5> </center> |
| Pearson R2 between low-coverage estimates and the high-coverage gene expression in a 10X dataset (subset of the Census of Immune cells). A) AGMAT gene, 50% downsampling (10,000 reads per cell); B) AGMAT gene, 10% downsampling (2,000 reads per cell); C) FYB1 gene, 50% downsampling (10,000 reads per cell); D) FYB1 gene, 10% downsampling (2,000 reads per cell). |


| <center> <h5 id="caption-figure-s4">Caption Figure S4</h5> </center> |
| Effective sample size as a function of number of individuals and number of cells per individual at budget $35,000 assuming no library preparation cost, multiplexing of 8 individuals per reaction, and known cell types. The dependence on read coverage is implicit. The maximum effective size Neff is a) 105 for B cells (N=120, M=2,750, r=14,500); b) 108 for CD14+ cells (N=120, M=2,750, r=14,500); c) 107 for CD4+ cells (N=120, M=2,750, r =14,500); d) 105 for CD8+ cells (N=120, M=2,750, r=14,500); e) 102 for dendritic cells (N=120, M=2,750, r=14,500); f) 106 for Fcgr3a cells (N=120, M=2,750, r =14,500); g) 105 for megakaryocytes (N=120, M=2,750, r=14,500); h) 106 for NK cells (N=120, M=2,750, r=14,500). | 






## ct-eQTL design calculator

The calculator for the design of single-cell RNA-Seq experiments for cell-type-specific eQTL studies is available [HERE](https://mandricigor.github.io/ct-eqtl-design/).























































