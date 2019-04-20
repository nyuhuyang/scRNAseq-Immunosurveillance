########################################################################
#
#  setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(sapply(c("Seurat","magrittr","SingleR","dplyr","reshape2",
                   "kableExtra","pheatmap","tidyr"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
args <- commandArgs(trailingOnly = TRUE)
#(load(file = "./data/MouseTumor_2_20190219.Rda"))
args = sapply(args,as.character)

(load(file="data/MouseTumor_2_20190219.Rda"))
(load(file="output/singler_T_MouseTumor_2_20190219.Rda"))
(load(file = paste0("data/",args[1])))
(load(file = paste0("output/",args[2])))
##############################
# add singleR label to Seurat
###############################
singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "singler1main" = singler$singler[[1]]$SingleR.single.main$labels,
                       "singler2sub" = singler$singler[[2]]$SingleR.single$labels,
                       "singler2main" = singler$singler[[2]]$SingleR.single.main$labels,
                       "orig.ident" = object@meta.data$orig.ident,
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels))

table(rownames(singlerDF) %in% object@cell.names)

apply(singlerDF,2,function(x) length(unique(x)))
object <- AddMetaData(object = object,
                   metadata = singlerDF)
object <- SetAllIdent(object = object, id = "singler1sub")

##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_sub1~.jpeg"), units="in", width=10, height=7,
     res=600)
SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single, top.n = 50,normalize = F,
                    cells_order= 0:14,
                    clusters = object@meta.data$res.0.6,order.by.clusters=T)
dev.off()

jpeg(paste0(path,"DrawHeatmap_sub1_N~.jpeg"), units="in", width=10, height=7,
     res=600)
SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single, top.n = 50,normalize = T,
                    cells_order= as.character(0:14),
                    clusters = object@meta.data$res.0.6,order.by.clusters=T)
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

table(object@meta.data$orig.ident, object@meta.data$singler2main) %>% t %>% kable %>%
        kable_styling() 

singler2main <- table(object@meta.data$orig.ident, object@meta.data$singler2main) %>% t %>% 
        as.data.frame %>% spread(key="Var2",value = "Freq")
singler2sub <- table(object@meta.data$orig.ident, object@meta.data$singler2sub) %>% t %>% 
        as.data.frame %>% spread(key="Var2",value = "Freq")

colnames(singler2main)[1] = "Main cell types"
colnames(singler2sub)[1] = "Sub cell types"

write.csv(singler2sub,paste0(path,"singler2sub.csv"))

##############################
# process color scheme
##############################

singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors2[duplicated(singler_colors2)]
length(singler_colors1)
apply(object@meta.data[,c("singler2sub","singler2main")],
      2,function(x) length(unique(x)))
object@meta.data[,c("singler2sub")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaColor(object = object, label= "singler2sub", colors = singler_colors1[1:12])
object <- SetAllIdent(object = object, id = "singler2sub")
TSNEPlot.1(object, colors.use = ExtractMetaColor(object),no.legend = F)

SplitTSNEPlot(object,split.by = "orig.ident",do.print = T, label.size = 3, do.label = F,no.legend = T)

##############################
# draw tsne plot
##############################
p3 <- TSNEPlot.1(object = object, do.label = T, group.by = "ident",
                 do.return = TRUE, no.legend = T,
                 colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 4,force = 2)+
  ggtitle("Supervised cell type labeling by Blueprint + Encode")+
  theme(text = element_text(size=10),							
        plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_sub.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

save(object,file=paste0("data/",args[1]))
##############################
# subset Seurat
###############################
SplitTSNEPlot(CD45,do.print = T,label.size = 3)

# remove CD45- contaminants
object %<>% SetAllIdent(id="res.0.6")
CD45 <- SubsetData(object,ident.remove = 10)
CD45 %<>% SetAllIdent(id="singler2sub")
CD45 <- SubsetData(CD45,ident.remove = "Granulocytes")
TSNEPlot(CD45)

p4 <- TSNEPlot.1(object = CD45, do.label = T, group.by = "ident",
                 do.return = TRUE, no.legend = T,
                 colors.use = ExtractMetaColor(CD45),
                 pt.size = 1,label.size = 4,force = 2)+
        ggtitle("Supervised cell type labeling by Blueprint + Encode")+
        theme(text = element_text(size=10),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_sub_CD45+.jpeg"), units="in", width=10, height=7,
     res=600)
print(p4)
dev.off()

