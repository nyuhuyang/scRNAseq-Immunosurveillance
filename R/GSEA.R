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
library(tidyverse)
library(readr)
library(openxlsx)
source("R/Seurat_functions.R")
path <- "output/20190422/"
if(!dir.exists(path)) dir.create(path, recursive = T)

# 5.1 load data ==============
(load(file = "data/MouseTumor_2_20190417.Rda"))

##################################
# Generate report from GSEA
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


##################################
# Generate figures from GSEA reports
##################################
data = list()
cell_types <- c("ALL_cell","Macrophages","Monocytes","NK_cells","T_cells")
# hallmark =======
for(i in seq_along(cell_types)){
        cell_type = cell_types[i]
        subfoder = paste(path,"GSEA",cell_type,"Hallmark/",sep="/")
        xls <- list.files(subfoder,pattern="gsea_report_for_.*xls")
        xls_list <- lapply(paste0(subfoder,"/",xls),function(x) {
                read_delim(x,"\t", escape_double = FALSE, trim_ws = TRUE)
        })
        df <- rbind.data.frame(xls_list[[1]],xls_list[[2]]) %>% 
                as_tibble() %>%
                arrange(desc(NES))
        colnames(df)[1] = "pathway"
        g <- ggplot(df, aes(reorder(pathway, NES), NES)) +
                geom_col(aes(fill=`NOM p-val`<0.05)) +
                coord_flip() +
                labs(x="pathway", y="Normalized Enrichment Score",
                     title=paste("Hallmark pathways NES in",cell_type,"in NAM conditions")) + 
                theme_minimal()
        jpeg(paste0(subfoder,"Hallmark_GSEA_",cell_type,".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        df$X12 = NULL
        data[[i]] = df
}
# GO =======
pathways <- read_delim("doc/GO_pathway.txt","\t", escape_double = FALSE, 
                       trim_ws = TRUE,col_names = F) %>% t %>% as.character()
for(i in seq_along(cell_types)){
        cell_type = cell_types[i]
        subfoder = paste(path,"GSEA",cell_type,"GO/",sep="/")
        xls <- list.files(subfoder,pattern="gsea_report_for_.*xls")
        xls_list <- lapply(paste0(subfoder,"/",xls),function(x) {
                read_delim(x,"\t", escape_double = FALSE, trim_ws = TRUE)
        })
        df <- rbind.data.frame(xls_list[[1]],xls_list[[2]]) %>% as.data.frame()
        df <- df[(df[,1] %in% pathways),] %>% 
                as_tibble() %>%
                arrange(desc(NES))
        colnames(df)[1] = "pathway"
        g <- ggplot(df, aes(reorder(pathway, NES), NES)) +
                geom_col(aes(fill=`NOM p-val`<0.05)) +
                coord_flip() +
                labs(x="pathway", y="Normalized Enrichment Score",
                     title=paste("GO pathways NES in",cell_type,"in NAM conditions")) + 
                theme_minimal()
        jpeg(paste0(subfoder,"GO_GSEA_",cell_type,".jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
        df$X12 = NULL
        data[[i+5]] = df
}
names(data) = paste0(rep(cell_types,time = 2),rep(c("_hallmark","_GO"),each = 5))
data = data[]
write.xlsx(data, file = paste0(path,"GSEA_results.xlsx"),
           colNames = TRUE, borders = "surrounding",colWidths = c(NA, "auto", "auto"))
