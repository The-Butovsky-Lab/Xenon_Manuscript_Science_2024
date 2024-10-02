# Code for Manuscript: 
This repository contains the analysis scripts and data used for the manuscript titled "Inhaled Xenon modulates microglia and ameliorates disease in mouse models of amyloidosis and tauopathy". Each script corresponds to a specific figure in the paper and is labeled accordingly. This code was used to generate all results and visualizations included in the manuscript.

![image](https://github.com/user-attachments/assets/b417793a-9104-4085-8334-df41d133dc3e)

## Files and Figures
The files are organized by figure number. Below is a brief description of the structure:

| Figure in Manuscript | Experiment Type | Script Name | GSE |
-----------------------|-----------------|------------|-------|
Fig.1 b-c | Xenon in APPPS1 - Mice | Fig1BC_APPPS1_Mouse_XenonvsAir.Rmd | GSE269157 |
Fig.1 l-m | Xenon in Microglia in vitro - Mice | Fig1L_Microglia_InVitro_Mouse.Rmd | GSE272910|
Fig. 3 | Xenon scRNAseq Myeloid | Fig3_scRNAseq_APPPS1_Myeloid_Mouse_XenonvsAir.Rmd| GSE274764 |
Fig 4. m-n | Xenon in IPSMG H5xFAD | Fig4M-N_Xenon in IPSMG_H5xFAD_MITRG_Mouse | GSE269481 |
Fig. 4 r, u and v | 5xFAD-MITRG Astrocytes 6M | Fig4R_U-V_5xFAD-MITRG Astrocytes 6M_Mouse.Rmd | GSE271423 |
Fig. 6 a-g | Xenon in IFNgr1 cKO | Fig6A_G_Xenon_IFNgr1cKO.Rmd| GSE271471 |
Fig. 6 i-k | Xenon in APPPS1 IFNy blocker-Mice | Fig6I_K_Xenon_APPPS1_IFNyblocker_Mouse.Rmd | GSE272684 |
Fig. 8  | Tau | Fig8_Tau_Mouse.Rmd | GSE275386 | 
Supp. Fig.1 d-e | Xenon PD in MGnD paradigm | SuppFig1D_E_Xenon_PD_MGnDparadigm.Rmd | GSE272685 |
Supp. Fig.2 a-f | Xenon PK in APP.PS1 microglia -mice | SuppFig2A-F_PK_APPPS1_Microglia_Mouse.Rmd | GSE272814 |
Supp. Fig. 3 d-e | Tdt Tomato - mice | SuppFig3D_E_Tdt_Tomato_Mouse.Rmd | GSE273541 |
Supp Fig 5. a-d | Xenon-iMGL in vitro | SuppFig5A-D_Xenon_iMGLInVitro.Rmd | GSE273575 |
Supp Fig.6 a-c | 5xFAD-MITRG-iMGL 6M | SuppFig6_A-C_5xFAD_MITRG-iMGL_6M.Rmd | GSE271472 | 
Supp Fig.7 c-e | Xenon in IFNy blocker -mice | SuppFig7C-E_Xenon_IFNyblocker.Rmd | GSE272686 |
Supp Fig 9 a-g | Xenon PK in APP.PS1 Tcells -mice | SuppFig9A-G_Xenon_PK_APPPS1_Tcells.Rmd | GSE272801 |

The SuperSeries containing all sub-GSEs is  GSE277025. Also, each script contains all necessary steps to reproduce the analysis and figures presented in the manuscript. Be sure to follow the usage instructions and install the required dependencies.

## Languages and Packages
The code is primarily written in the following languages:
* R (v4.2.1): Used for data analysis and visualization. <br>
* Shell: Used for pre-processing the fastq files

### Major R Packages: <br>
* ggplot2 – for data visualization <br>
* Seurat – for single-cell RNA sequencing analysis <br>
* biomaRt – for querying gene information <br>
* fgsea – for gene set enrichment analysis <br>

## Contact
For any questions or further information, please contact:

* Madison Carpenter
* Email: mmcarpenter@bwh.harvard.edu
* GitHub: https://github.com/madison-car

## Citation 
