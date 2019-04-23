########################################################################
#
#  0 setup environment, install libraries if nLynchessary, load libraries
# 
# ######################################################################

invisible(sapply(c("Seurat","magrittr","tidyr","dplyr","kableExtra","gplots","pheatmap"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("./data/")) dir.create("data")

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters.
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
args <- commandArgs(trailingOnly = TRUE)
(load(file = "./data/MouseTumor_2_20190417.Rda"))
#args[1] = as.character(args[1])
#(load(file = args[1]))

#head(gde.all,10) %>% kable %>% kable_styling

# FindPairMarkers
(ident <- unique(object@meta.data$manual) %>% sort)
print(ident.1 <- paste0(ident, "_NAM"))
print(ident.2 <- paste0(ident, "_Control"))
object@meta.data$manual_orig.ident = paste0(object@meta.data$manual, "_", object@meta.data$orig.ident)
object %<>% SetAllIdent(id = "manual_orig.ident")
table(object@ident)
subfolder <- paste0(path,"NAM_vs_Control/")
gde.pair <- FindPairMarkers(object, ident.1 = ident.1, ident.2 = ident.2,
                               logfc.threshold = 0.05, 
                               return.thresh = 0.05, only.pos = FALSE, save.path = subfolder)
write.csv(gde.pair, paste0(subfolder,"pairwise_comparision.csv"))
head(gde.pair[gde.pair$cluster1.vs.cluster2 == "NK cells_NAM vs.NK cells_Control",],10) %>% 
        kable %>% kable_styling


#  DoHeatmap ============
object %<>% SetAllIdent(id = "manual")
table(object@ident)
if(!dir.exists(paste0(subfolder,"DoHeatmap/"))) dir.create(paste0(subfolder,"DoHeatmap/"), recursive = T)
for(id in ident){
        id.pair <- read.csv(paste0(subfolder,id,"_NAM vs.",id,"_Control.csv"), header = 1, row.names = 1)
        object_subset <- SubsetData(object,ident.use = id)
        object_subset %<>% ScaleData() %>% 
                SetAllIdent(id = "orig.ident")
        g <- DoHeatmap.1(object_subset, id.pair, Top_n = 50,
                         use.scaled = T, 
                         #col.low = "#FFFF00",col.mid = "#FFA500", col.high = "#FF0000",
                         ident.use = paste("Control vs. NAM in", id),
                         group.label.rot = F,cex.row = 4,remove.key =F,title.size = 12)
        jpeg(paste0(subfolder,"DoHeatmap/",id,"_NAM vs.",id,"_Control~.jpeg"), units="in", width=10, height=7,
             res=600)
        print(g)
        dev.off()

}

#  heatmap.2 ============
object %<>% SetAllIdent(id = "manual")
table(object@ident)
Top_n =50
if(!dir.exists(paste0(subfolder,"heatmap2/"))) dir.create(paste0(subfolder,"heatmap2/"), recursive = T)
for(id in ident){
        id.pair <- read.csv(paste0(subfolder,id,"_NAM vs.",id,"_Control.csv"), header = 1, row.names = 1)
        object_subset <- SubsetData(object,ident.use = id)
        object_subset %<>% ScaleData()
        colnames(id.pair)[8] = "cluster"
        top <-  id.pair %>%  group_by(cluster) %>% top_n(Top_n, avg_logFC)
        bottom <-  id.pair %>% group_by(cluster) %>% top_n(Top_n, -avg_logFC)
        all <- c(as.character(top$gene), as.character(bottom$gene))
        y = object_subset@scale.data[unique(all),]
        hc <- hclust(as.dist(1-cor(as.matrix(y), method="spearman")), method="complete")
        cc = gsub("_.*","",hc$labels)
        cc = gsub("Control","#B3DE69",cc)
        cc = gsub("NAM","#195016",cc)
        
        jpeg(paste0(subfolder,"heatmap2/",id,"_NAM vs.",id,"_Control.jpeg"), units="in", width=10, height=7,res=600)
        heatmap.2(as.matrix(y),
                  Colv = as.dendrogram(hc), Rowv= F,
                  ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",#scale="row",
                  key.xlab = "scale log nUMI",
                  cexRow = 0.5,
                  margins = c(2,5),
                  #scale = "row",
                  breaks = seq(-3,3,length.out = 300),
                  col = rev(heat.colors(299)),
                  main = paste0(id,"_NAM vs.",id,"_Control"))
        par(lend = 1)           # square line ends for the color legend
        legend(0, 0.8,       # location of the legend on the heatmap plot
               legend = c("Control", "NAM"), # category labels
               col = c("#B3DE69", "#195016"),  # color key
               lty= 1,             # line style
               lwd = 10            # line width
        )
        dev.off()
}
        
# Identify cell type markers ================
object %<>% SetAllIdent(id="orig.ident")
gde.all <- FindAllMarkers.UMI(object,logfc.threshold = 0.1, only.pos = F,return.thresh = 0.05, 
                              test.use = "MAST")
write.csv(gde.all, paste0(subfolder,"All_NAM vs. All_Control.csv"))
