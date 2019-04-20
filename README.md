# scRNAseq-LorenzoBlood
CD45+ tumor infiltrate of mouse tumors (2 conditions: control versus treatment) single cell RNA-seq data anlaysis for Lorenzo Galluzzi(log3001@med.cornell.edu) and Aitziber Buque Martinez (abm2013@med.cornell.edu)

Rscript R/Rscript/QC.R "doc/190125_scRNAseq_info.xlsx"
Rscript R/Rscript/scater.R "doc/190125_scRNAseq_info.xlsx"
Rscript R/Rscript/Identify_Cell_Types_Manually.R "data/MouseTumor_2_{date}.Rda"
