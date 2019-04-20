########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
library(Seurat)
library(SingleR)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# 5.1 load data ==============
(load(file = "./data/MouseTumor_2_20190417.Rda"))

##################################
# All cell type
##################################
object %<>% SetAllIdent(id="orig.ident")
PrepareGSEA(object, k = 100)
ReportGSEA(file = "ALL_cell.h.all.NAM_versus_Control.Gsea.1555732506514",pos=T)
ReportGSEA(file = "ALL_cell.GO.NAM_versus_Control.Gsea.1555732996935",pos=F)


object %<>% SetAllIdent(id="manual")
table(object@ident)
Macrophages <- SubsetData(object, ident.use = c("Macrophages")) %>% 
        SetAllIdent(id="orig.ident")
PrepareGSEA(Macrophages, k = 100)
ReportGSEA(file = "Macrophages.Hallmark.NAM_versus_Control.Gsea.1555733693440",pos=F)
ReportGSEA(file = "Macrophages.GO.NAM_versus_Control.Gsea.1555733933549",pos=T)

Monocytes <- SubsetData(object, ident.use = c("Monocytes")) %>% 
        SetAllIdent(id="orig.ident")
PrepareGSEA(Monocytes, k = 100)
ReportGSEA(file = "Monocytes.Hallmark.NAM_versus_Control.Gsea.1555734171475",pos=F)
ReportGSEA(file = "NK_cells.GO.NAM_versus_Control.Gsea.1555736648009",pos=F)

NK_cells <- SubsetData(object, ident.use = c("NK cells")) %>% 
        SetAllIdent(id="orig.ident")
PrepareGSEA(NK_cells, k = 100)
ReportGSEA(file = "NK_cells.Hallmark.NAM_versus_Control.Gsea.1555735480255",pos=F)
ReportGSEA(file = "GO.Lymphoid.NAM_versus_Control.Gsea.1554222618473",pos=T)

T_cells <- SubsetData(object, ident.use = c("T cells")) %>% 
        SetAllIdent(id="orig.ident")
PrepareGSEA(T_cells, k = 100)
ReportGSEA(file = "T_cells.Hallmark.NAM_versus_Control.Gsea.1555736788973",pos=F)
ReportGSEA(file = "T_cells.GO.NAM_versus_Control.Gsea.1555737179894",pos=T)
