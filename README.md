# Optimal design for cell-type specific eQTL studies

<p align="center">

<img src="./figures/Figure-S14.png" width="50%">

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
  - [Scripts](#scripts)

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
| [Caption](#caption-figure-1)  | [Caption](#caption-figure-2) | [Caption](#caption-figure-3) | [Caption](#caption-figure-4) |

|            Figure 5            |            Figure 6            |                                           Figure 7                                            |                                             Figure S1                                             |
| :----------------------------: | :----------------------------: | :-------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------: |
| ![](./figures/Figure-5.png) | ![](./figures/Figure-S1.png) |                                ![](./figures/Figure-S2.png)                                 |                                  ![](./figures/Figure-S3.png)                                  |
| [Caption](#caption-figure-5)  | [Caption](#caption-figure-s1) | [Caption](#caption-figure-s2) | [Caption](#caption-figure-s3) |

|               Figure S2                |             Figure S3              |                 Figure S4                 |                                          Figure S6                                           |
| :------------------------------------: | :--------------------------------: | :---------------------------------------: | :------------------------------------------------------------------------------------------: |
|    ![](./figures/Figure-S4.jpg)     |  ![](./figures/Figure-S5.jpg)   |      ![](./figures/Figure-S6.jpg)      |                               ![](./figures/Figure-S7.jpg)                                |
| [Caption](#caption-figure-s4)  | [Caption](#caption-figure-s5) | [Caption](#caption-figure-s6) | [Caption](#caption-figure-s7) |

|               Figure S2                |             Figure S3              |                 Figure S4                 |                                          Figure S6                                           |
| :------------------------------------: | :--------------------------------: | :---------------------------------------: | :------------------------------------------------------------------------------------------: |
|    ![](./figures/Figure-S8.png)     |  ![](./figures/Figure-S9.png)   |      ![](./figures/Figure-S10.png)      |                               ![](./figures/Figure-S11.png)                                |
| [Caption](#caption-figure-s8)  | [Caption](#caption-figure-s9) | [Caption](#caption-figure-s10) | [Caption](#caption-figure-s11) |


|               Figure S2                |             Figure S3              |                 Figure S4                 |                                          Figure S6                                           |
| :------------------------------------: | :--------------------------------: | :---------------------------------------: | :------------------------------------------------------------------------------------------: |
|    ![](./figures/Figure-S12.png)     |  ![](./figures/Figure-S13.png)   |      ![](./figures/Figure-S14.png)      |                               ![](./figures/Figure-S15.png)                                |
| [Caption](#caption-figure-s12)  | [Caption](#caption-figure-s13) | [Caption](#caption-figure-s14) | [Caption](#caption-figure-s15) |


|               Figure S2                |             Figure S3              |                 Figure S4                 |                                          Figure S6                                           |
| :------------------------------------: | :--------------------------------: | :---------------------------------------: | :------------------------------------------------------------------------------------------: |
|    ![](./figures/Figure-S16.png)     |  ![](./figures/Figure-S17.png)   |      ![](./figures/Figure-S18.png)      |                               ![](./figures/Figure-S19.png)                                |
| [Caption](#caption-figure-s16)  | [Caption](#caption-figure-s17) | [Caption](#caption-figure-s18) | [Caption](#caption-figure-s19) |



## Figure captions

| <center> <h3 id="caption-figure-1">Caption Figure 1</h3> </center> |
| :------------------------------------: |
| The impact of read coverage on the average R2 between cell-type-specific gene expression estimates and their ground truth values (Smart-Seq2 dataset, alpha cells).  A) Average Pearson R2 (± 1 standard error) computed across all the genes at different levels of read coverage, Smart-Seq2 dataset. B) Pearson R2 at 75,000 reads per cell (± 1 standard error) stratified by the expression level, Smart-Seq2 dataset. C) Average Pearson R2 (± 1 standard error) computed across all the genes at different levels of read coverage, 10X dataset. D) Pearson R2 at 4,000 reads per cell (± 1 standard error) stratified by the expression level, Smart-Seq2 dataset. |

| # Caption Figure 2 |
| :------------------------------------: |
| The impact of read coverage on the accuracy of demultiplexing.  Shown here is the percentage of correctly assigned cells (to the corresponding sample) as a function of average number of reads per cell. |

| # Caption Figure 3 |
| :------------------------------------: |
| Experimental designs for CD4 T cells ct-eQTL with effective sample size Neff=40: A) Comparison of different experimental designs. Experimental design N=88, M=2,250, r=4,500 yields two-fold reduction in cost than the standard design. B) For a fixed sample size and number of cells per individual, increasing coverage implies increasing the effective sample size (i.e., power) only up to a point. There is little gain in power at coverages greater than 12,500 reads per cell. |

| # Caption Figure 4 |
| :------------------------------------: |
| Figure 4: Effective sample size as a function of cell-type prevalence. Shown here is the effective sample size across the grid of experimental design when the cell-type abundance is set to different values - 5, 10, 15, 20, 25, 30%.  (CD4 T cells at fixed budget $35,000). |

| # Caption Figure 5 |
| :------------------------------------: |
| Figure 5: Performance of ct-eQTL analysis. Shown here is recall (power estimate) as a function of coverage in the ct-eQTL analysis of CD4 T cells at fixed budget $35,000. A) Mean ct-eQTL; B) Variance ct-eQTL. |

| # Caption Figure 2 |
| :------------------------------------: |
| The impact of read coverage on the accuracy of demultiplexing.  Shown here is the percentage of correctly assigned cells (to the corresponding sample) as a function of average number of reads per cell. |

| # Caption Figure 2 |
| :------------------------------------: |
| The impact of read coverage on the accuracy of demultiplexing.  Shown here is the percentage of correctly assigned cells (to the corresponding sample) as a function of average number of reads per cell. |

| # Caption Figure 2 |
| :------------------------------------: |
| The impact of read coverage on the accuracy of demultiplexing.  Shown here is the percentage of correctly assigned cells (to the corresponding sample) as a function of average number of reads per cell. |

| # Caption Figure 2 |
| :------------------------------------: |
| The impact of read coverage on the accuracy of demultiplexing.  Shown here is the percentage of correctly assigned cells (to the corresponding sample) as a function of average number of reads per cell. |

| # Caption Figure 2 |
| :------------------------------------: |
| The impact of read coverage on the accuracy of demultiplexing.  Shown here is the percentage of correctly assigned cells (to the corresponding sample) as a function of average number of reads per cell. |

| # Caption Figure 2 |
| :------------------------------------: |
| The impact of read coverage on the accuracy of demultiplexing.  Shown here is the percentage of correctly assigned cells (to the corresponding sample) as a function of average number of reads per cell. |

| # Caption Figure 2 |
| :------------------------------------: |
| The impact of read coverage on the accuracy of demultiplexing.  Shown here is the percentage of correctly assigned cells (to the corresponding sample) as a function of average number of reads per cell. |

| # Caption Figure 2 |
| :------------------------------------: |
| The impact of read coverage on the accuracy of demultiplexing.  Shown here is the percentage of correctly assigned cells (to the corresponding sample) as a function of average number of reads per cell. |

| # Caption Figure 2 |
| :------------------------------------: |
| The impact of read coverage on the accuracy of demultiplexing.  Shown here is the percentage of correctly assigned cells (to the corresponding sample) as a function of average number of reads per cell. |



## Scripts

