# **Study CD45+ tumor infiltrate of mouse tumors with single-cell RNA-seq**

This project provides the code developed in the study of Aitziber Buqué _et al._ **_"Immunosurveillance of HR+ 1 breast cancer - Prophylactic and therapeutic effects of nicotinamide"_**. It focuses on comparing control vs. nicotinamide (NAM) treated CD45+ tumor infiltrate of mouse tumors with single-cell RNA-seq.

## **Requirements**

* R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
* R librarys: DropletUtils_1.2.2, IRanges_2.16.0, R.methodsS3_1.7.1, Biobase_2.42.0, BiocGenerics_0.28.0, BiocParallel_1.16.6, cowplot_0.9.4, DelayedArray_0.8.0, dplyr_0.8.1, fgsea_1.8.0, GenomeInfoDb_1.18.2, GenomicRanges_1.34.0, ggplot2_3.1.1, ggpubr_0.2, gplots_3.0.1.1, harmony_0.1.0, kableExtra_1.1.0, magrittr_1.5, Matrix_1.2-15, matrixStats_0.54.0, pheatmap_1.0.12, R.oo_1.22.0, R.utils_2.8.0, Rcpp_1.0.1, readr_1.3.1, reshape2_1.4.3, S4Vectors_0.20.1, scater_1.10.1, Seurat_2.3.4, SingleCellExperiment_1.4.1, SingleR_0.2.2, SummarizedExperiment_1.12.0, tidyr_0.8.3

## **Data**

Raw counts data and processed Seurat object will be released after final publication.

## **Reproduce results**

#### **1-3. Data preprocess**

[1 QC.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/Rscript/QC.R)

Download cell ranger results in `~/Downloads` folder. Open terminal, enter the _Code_ root directory and run `Rscript R/Rscript/QC.R "doc/190125_scRNAseq_info.xlsx"`. This script moves all necessary files to `data` folder, run initial quality control. Output ggplot is stored in `output/{date}/g1_2_{date}.Rda`

[2 scater.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/Rscript/scater.R)

Run `Rscript R/Rscript/scater.R "doc/190125_scRNAseq_info.xlsx"` in terminal. This script uses scater pipeline to to process the single-cell RNA-seq data. The output SingleCellExperiment object is stored in `data/sce_2_{date}.Rda`

[3 Seurat_setup.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/Rscript/Seurat_setup.R)

Run `Rscript R/Rscript/Seurat_setup.R "doc/190125_scRNAseq_info.xlsx"` in terminal. This script uses Seurat 2 to read the result from the previous scater pipeline, perform normalization, scaling, dimension reduction, and unsupervised-clustering. The output Seurat object is saved in file `data/MouseTumor_2_{date}.Rda`

Above three scripts assume they are run on the same date. Otherwise the `{date}` in the file names must be updated accordingly.


#### **4-5. Identify cell types**

[4 Identify_Cell_Types_Manually.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/Rscript/Identify_Cell_Types_Manually.R)

An R shiny app is built to identify cell types. [scRNAseq-Immunosurveillance](https://weillcornellmed.shinyapps.io/scRNAseq-Immunosurveillance/)


Run `Rscript R/Rscript/Identify_Cell_Types_Manually.R "data/MouseTumor_2_{date}.Rda"`. This script use predefinde cell type markers to manually identify cell types.

[5 SingleR.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/Rscript/SingleR.R)

Run `Rscript R/Rscript/SingleR.R "data/MouseTumor_2_{date}.Rda"`. This script use SingleR package to identify cell types based reference datasets.

#### **6-9. Characterize CD45+ tumor infiltrate of mouse tumors**

[6 Figures.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/Figures.R), removing unwanted cells, relabel the clustering, generate tsne plot and bar plots. 

For more information see
[_Generate TSNE plots and compare gene expression in individual cell types (figure 8a-b)_](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/wiki/1.-Generate-TSNE-plots-and-compare-gene-expression-in-individual-cell-types-(figure-8a,-b))

[7 Differential_analysis.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/Differential_analysis.R), conducting differential analysis between control and NAM treated tumor sample, and generate heatmaps.

[8 FGESA.R](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/blob/master/R/FGESA.R), performing gene set enrichment analysis based on the results from previous step 7. 

For more information, see [_2. Differential analysis and gene set enrichment analysis (figure S7, S8)_](https://github.com/nyuhuyang/scRNAseq-Immunosurveillance/wiki/2.-Differential-analysis-and-gene-set-enrichment-analysis-(figure-S7, S8))