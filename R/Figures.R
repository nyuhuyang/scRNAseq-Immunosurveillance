########################################################################
#
#  setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(sapply(c("Seurat","magrittr","SingleR","dplyr","reshape2","ggpubr",
                   "kableExtra","pheatmap","tidyr","readr"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== load data ==========================================
(load(file="data/MouseTumor_2_20190219.Rda"))
(load(file="output/singler_T_MouseTumor_2_20190219.Rda"))

#====== raw TSNE ==========================================
object %<>% SetAllIdent("res.0.6")
g <- TSNEPlot.1(object = object, do.label = F, group.by = "ident",
                        do.return = TRUE, no.legend = F, 
                        #colors.use = ExtractMetaColor(object),
                        pt.size = 1,label.size = 6 )+
        ggtitle("tSNE plot of all clusters")+
        theme(plot.title = element_text(hjust = 0.5,size = 25),
              text= element_text(size = 20))

jpeg(paste0(path,"TSNEplot-seurat.jpeg"), units="in", width=10, height=7,res=600)
print(g)
dev.off()

xy = object@dr$tsne@cell.embeddings
out = SingleR.PlotTsne(singler$singler[[2]]$SingleR.single,xy,do.label=FALSE,
                       do.letters =F,labels=singler$meta.data$clusters, label.size=5,
                       dot.size = 2,alpha=1)
jpeg(paste0(path,"TSNEplot-singer.jpeg"), units="in", width=10, height=7,res=600)
out$p+ggtitle("tSNE plot of all clusters")+
        theme(plot.title = element_text(hjust = 0.5,size = 25, face = "plain"),
              text= element_text(size = 20),
              legend.text = element_text(size=15))+
        guides(colour = guide_legend(override.aes = list(size=6)),
               shape = guide_legend(override.aes = list(size=6)))
dev.off()

#====== TSNE ==========================================
jpeg(paste0(path,"DrawHeatmap_sub1~.jpeg"), units="in", width=10, height=7,
     res=600)
SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single, top.n = 50,normalize = F,
                    clusters = object@meta.data$res.0.6,order.by.clusters=T)
dev.off()

jpeg(paste0(path,"DrawHeatmap_sub1_N~.jpeg"), units="in", width=10, height=7,
     res=600)
SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single, top.n = 50,normalize = T,
                    clusters = object@meta.data$res.0.6,order.by.clusters=T)
dev.off()

#====== cell type TSNE ==========================================
# manually select NK and T cells
object %<>% SetAllIdent("res.0.6")
p <- TSNEPlot.1(object = object, do.label = F, group.by = "ident",
                 do.return = TRUE, no.legend = F,
                 colors.use = ExtractMetaColor(object),
                 pt.size = 1,label.size = 4,force = 2)+
        ggtitle("Supervised cell type labeling by Blueprint + Encode")+
        theme(text = element_text(size = 20),							
              plot.title = element_text(hjust = 0.5,size = 25),
              legend.text = element_text(size=15))
jpeg(paste0(path,"TSNEplot-All-CellTypes.jpeg"), units="in", width=10, height=7,res=600)
print(p)
dev.off()

jpeg(paste0(path,"TSNEplot-All-CellTypes-no-Legend.jpeg"), units="in", width=10, height=7,res=600)
print(p+theme(legend.position = "none"))
dev.off()
##############################
# Adjust cell label, keep T/NK, macrophages and monocytes only
##############################
TSNEPlot.1(object, do.label = T,colors.use = ExtractMetaColor(object))
object %<>% SetAllIdent(id = "singler2sub")
object %<>% SubsetData(ident.remove = c("Granulocytes","B cells"))
object %<>% SetAllIdent("res.0.6")
TSNEPlot(object, do.label = T)
object %<>% SubsetData(ident.remove = c(10,14))

T_cells <- SubsetData(object, ident.use = c(8,9,12))
T_cells_id_copy <- FeaturePlot(object, features.plot = "Cd3g", do.identify = T)
T_cells_id = unique(c(T_cells_id,T_cells_id_copy))
length(T_cells_id)
remove(T_cells);GC()

object@meta.data$manual <- object@meta.data$singler2sub
object@meta.data[T_cells_id,"manual"] = "T cells"
object@meta.data$manual = gsub("Macrophages activated","Macrophages",object@meta.data$manual)
object@meta.data$manual = gsub("Dendritic cells","Monocytes",object@meta.data$manual)
table(object@meta.data$manual)
table(object@meta.data$singler2sub)
##############################
# process color scheme
##############################

singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors2[duplicated(singler_colors2)]
length(singler_colors2)
apply(object@meta.data[,c("singler2sub","manual")],
      2,function(x) length(unique(x)))
object@meta.data[,c("manual")] %>% table() %>% kable() %>% kable_styling()
object <- AddMetaColor(object = object, label= "manual", colors = singler_colors2)
object <- SetAllIdent(object = object, id = "manual")
TSNEPlot.1(object, colors.use = ExtractMetaColor(object),no.legend = F,do.print = T)

SplitTSNEPlot(object,split.by = "orig.ident",do.print = T, label.size = 3, do.label = F,no.legend = T)

# four cell types
p1 <- TSNEPlot.1(object = object, do.label = F, group.by = "ident",
                do.return = TRUE, no.legend = F,
                colors.use = ExtractMetaColor(object),
                pt.size = 1,label.size = 4,force = 2)+
        ggtitle("Myeloid and Lymphoid")+
        theme(text = element_text(size = 20),							
              plot.title = element_text(hjust = 0.5,size = 25),
              legend.text = element_text(size=15))
jpeg(paste0(path,"TSNEplot-Myeloid-and-Lymphoid.jpeg"), units="in", width=10, height=7,res=600)
print(p1)
dev.off()

jpeg(paste0(path,"TSNEplot-Myeloid-and-Lymphoid-no-Legend.jpeg"), units="in", width=10, height=7,res=600)
print(p1+theme(legend.position = "none"))
dev.off()

save(object, file = paste0("data/","MouseTumor_",length(sample_n),
                           "_",gsub("-","",Sys.Date()),".Rda"))

#======= cell number ===========
(load(file="data/MouseTumor_2_20190417.Rda"))
df <- table(object@meta.data$orig.ident, object@meta.data$manual) %>% t %>% 
        as.data.frame %>% spread(key="Var2",value = "Freq")

colnames(df)[1] = "Sub cell types"
rownames(df) = df[,1]
df = df[,-1]
df
p_value = c()
ColSum <- colSums(df)

for(i in 1:nrow(df)){
        conting <- rbind(df[i,],ColSum-df[i,])
        FISH <- fisher.test(conting,conf.int = T)
        p_value[i] = FISH$p.value
        #CHI = chisq.test(conting, correct = T)
        #chisq_p_value[i] = CHI$p.value             
}

df$p_value = p_value
df$p_val_adj = p.adjust(p = df$p_value, method = "bonferroni", 
                        n = nrow(df))
df %>% kable %>% kable_styling()
write.csv(df,paste0(path,"cell_numbers.csv"))

# ========= bar chart =============
# load seurat
(load(file="data/MouseTumor_2_20190417.Rda"))
object %<>% SetAllIdent(id="manual")
table(object@ident)
# load gene list
barchart_genes <- read_delim("doc/barchart_genes.txt","\t", escape_double = FALSE,
                             trim_ws = TRUE) %>% nest(-cell_types)
(cell_types <- barchart_genes$cell_types)
barchart_path <- paste0(path,"barchart/")
if(!dir.exists(barchart_path)) dir.create(barchart_path, recursive = T)

for(i in 1:length(cell_types)){
        print(cell_types[i])
        sub_object <- SubsetData(object, ident.use = cell_types[i])
        genes <- barchart_genes$data[[i]] %>% t %>% as.vector
        genes <- FilterGenes(object, genes)
        print(genes)
        sub_object@meta.data[,genes] = t(sub_object@data[genes,])
        sub_object_meta.data = sub_object@meta.data[,c("orig.ident",genes)]
        colnames(sub_object_meta.data)[1] = "Conditions"
        colnames(sub_object_meta.data) = gsub("-","",colnames(sub_object_meta.data))
        genes1 = gsub("-","",genes)
        for(k in 1:length(genes)){
                jpeg(paste0(barchart_path,paste(cell_types[i],collapse = "_"),"_",genes[k],
                            ".jpeg"), units="in", width=10, height=7,res=600)
                p <- ggbarplot(sub_object_meta.data, x = "Conditions", y = genes1[k], 
                               add = c("mean_se", "jitter"),ylab = paste(genes[k], "(UMI)"),
                               color = "Conditions", palette = c("black", "red"))+
                        stat_compare_means(label.x = 1.5, 
                                           label.y = max(sub_object_meta.data[,genes1[k]])+0.25)
                print(p)
                dev.off()
                print(paste0(k,":",length(genes)))
        }
}
