########################################################################
#
#  setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
invisible(sapply(c("magrittr","SingleR","dplyr","genefilter","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("../R/Seurat_functions.R")
source("../R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# check GBMLGG exprs data==============================
GBMLGG_RSEM <- read.delim2(paste("data","GBMLGG.mRNAseq_Preprocess.Level_3",
                                   "GBMLGG.uncv2.mRNAseq_RSEM_all.txt",
                                   sep = "/"), header = TRUE, sep= "\t",
                             stringsAsFactors = F)
gene_name <- GBMLGG_RSEM[,1] %>% as.vector
gene_name <- sub("[^[:alnum:]]", " ", gene_name)

gene_name <- gsub("  "," ",gene_name)
gene_name_df <- strsplit(gene_name, split = " ") %>% unlist %>% 
        matrix(ncol = 2, byrow = TRUE)
GBMLGG_RSEM[,1] <- gene_name_df[,1]
GBMLGG_RSEM <- GBMLGG_RSEM[(GBMLGG_RSEM[,1] != ""),]
colnames(GBMLGG_RSEM)[1] = "gene"

# remove NA rows
table(is.na(GBMLGG_RSEM))
#GBMLGG_counts <- GBMLGG_counts[!apply(GBMLGG_counts,1, function(x) all(is.na(x))),]

#remove duplicate rownames with lower rowsumns
#' @param mat input as data.frame with gene name
#' @export mat matrix with gene as rownames, no duplicated genes
RemoveDup <- function(mat){
        gene_id <- as.matrix(mat[,1])
        mat <- mat[,-1]
        if(!is.matrix(mat)) mat <- sapply(mat,as.numeric)
        rownames(mat) <- 1:nrow(mat)
        mat[is.na(mat)] = 0
        mat <- cbind(mat, "rowSums" = rowSums(mat))
        mat <- mat[order(mat[,"rowSums"],decreasing = T),]
        gene_id <- gene_id[as.numeric(rownames(mat))]
        remove_index <- duplicated(gene_id)
        mat <- mat[!remove_index,]
        rownames(mat) <- gene_id[!remove_index]
        return(mat[,-ncol(mat)])
}
GBMLGG_RSEM <- RemoveDup(GBMLGG_RSEM)
dim(GBMLGG_RSEM)
# rename colnames 
colnames(GBMLGG_RSEM) <- gsub("\\.01$|\\.11$","",colnames(GBMLGG_RSEM))
colnames(GBMLGG_RSEM) <- gsub("\\.","-",colnames(GBMLGG_RSEM))
GBMLGG_RSEM[1:5,1:5]
dup <- grep("-02$",colnames(GBMLGG_RSEM),value = T)
testMMM(GBMLGG_RSEM)

jpeg(paste0(path,"boxplot_GBMLGG_RSEM.jpeg"), units="in", width=10, height=7,res=600)
boxplot(log1p(GBMLGG_RSEM)) #slow
title(main = "boxplot for GBM and LGG")
dev.off()

# check GBMLGG annotation data==============================
subtype <- read.delim2(paste("data","GBMLGG.mRNAseq_Preprocess.Level_3",
                                   "tumormap.ucsc.edu.subtype.txt",
                                   sep = "/"), header = TRUE, sep= "\t",
                             stringsAsFactors = F)
samples_df <- strsplit(subtype$samples, split = "\\.") %>% unlist %>% 
        matrix(ncol = 2, byrow = TRUE)
subtype_df <- cbind.data.frame(samples_df,"subtype" = subtype$subtype)
colnames(subtype_df) <- c("samples","rna_methy","subtype")
table(subtype_df$rna_methy)
subtype_df[(subtype_df$samples %in% colnames(GBMLGG_RSEM)),"rna_methy"]
# add .02 duplicate
subtype_df_02<- subtype_df[(subtype_df$samples %in% gsub("-02$","",dup)),]

subtype_df_02$samples <- paste0(subtype_df_02$samples,"-02")
df_subtype <- rbind.data.frame(subtype_df,subtype_df_02)
rownames(df_subtype) <- df_subtype$samples
df_subtype <- df_subtype[(df_subtype$samples %in% colnames(GBMLGG_RSEM)),-1]
dim(df_subtype)
table(df_subtype$subtype)

df_new_subtype <- df_subtype[(df_subtype$subtype %in% c("Classical", "Mesenchymal",
                                                    "Neural", "Proneural","IDHwt")),]
dim(df_new_subtype)
GBM_RSEM <- GBMLGG_RSEM[,(colnames(GBMLGG_RSEM) %in% rownames(df_new_subtype))]
dim(GBM_RSEM)
df_new_subtype <- df_new_subtype[colnames(GBM_RSEM),]
# Create Singler Reference
ref = CreateSinglerReference(name = 'GBM_RSEM',
                             expr = log1p(GBM_RSEM), # the expression matrix
                             types = df_new_subtype$subtype, 
                             main_types = df_new_subtype$subtype)

save(ref,file='data/GeneSets/ref_GBM_RSEM.RData') # it is best to name the object and the file with the same name.



# check GBMLGG exprs counts data==============================
GBMLGG_counts <- read.delim2(paste("data","GBMLGG.mRNAseq_Preprocess.Level_3",
                                   "GBMLGG.uncv2.mRNAseq_raw_counts.txt",
                                   sep = "/"), header = TRUE, sep= "\t",
                             stringsAsFactors = F)
gene_name <- GBMLGG_counts[,1] %>% as.vector
gene_name <- sub("[^[:alnum:]]", " ", gene_name)
gene_name <- gsub("  "," ",gene_name)
gene_name_df <- strsplit(gene_name, split = " ") %>% unlist %>% 
        matrix(ncol = 2, byrow = TRUE)
GBMLGG_counts[,1] <- gene_name_df[,1]
GBMLGG_counts <- GBMLGG_counts[(GBMLGG_counts[,1] != ""),]
colnames(GBMLGG_counts)[1] = "gene"
GBMLGG_counts <- RemoveDup(GBMLGG_counts)
dim(GBMLGG_counts)
# rename colnames 
colnames(GBMLGG_counts) <- gsub("\\.01$|\\.11$","",colnames(GBMLGG_counts))
colnames(GBMLGG_counts) <- gsub("\\.","-",colnames(GBMLGG_counts))
GBMLGG_counts[1:5,1:5]
GBMLGG_counts <- GBMLGG_counts[,(colnames(GBMLGG_counts) %in% rownames(df_subtype))]
dim(GBMLGG_counts)
dup <- grep("-02$",colnames(GBMLGG_counts),value = T)

#=============Seruat========================
GBMLGG <- CreateSeuratObject(raw.data = GBMLGG_counts, project = "TCGA_GBM_LGG",
                           min.cells = 3) %>%
        NormalizeData(scale.factor = 1000) %>%
        ScaleData(display.progress = F)
VlnPlot(object = GBMLGG, features.plot = c("nGene", "nUMI"), nCol = 2)
GBMLGG <- FindVariableGenes(GBMLGG, x.low.cutoff = 0.0125, do.plot = FALSE)
length(x = GBMLGG@var.genes)

GBMLGG <- RunPCA(GBMLGG, pc.genes = GBMLGG@var.genes, pcs.compute = 50, 
                 do.print = F) %>%
        RunTSNE(reduction.use = "pca", dims.use = 1:10, do.fast = T) %>% 
        FindClusters(reduction.type = "pca", resolution = 0.6, dims.use = 1:10, 
                         force.recalc = TRUE, print.output = FALSE)

GBMLGG <- AddMetaData(object = GBMLGG, metadata = df_subtype)
head(GBMLGG@meta.data)
##############################
# process color scheme
##############################
singler_colors <- readxl::read_excel("./doc/singler.colors.xlsx")
singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors1[duplicated(singler_colors1)];singler_colors2[duplicated(singler_colors2)]
length(singler_colors1);length(singler_colors2)
length(unique(GBMLGG@meta.data[,"subtype"]))
GBMLGG@meta.data[,c("subtype")] %>% table() %>% kable() %>% kable_styling()
GBMLGG <- AddMetaColor(object = GBMLGG, label= "subtype", colors = singler_colors1[1:8])
GBMLGG <- SetAllIdent(object = GBMLGG, id = "subtype")
p3 <- TSNEPlot.1(object = GBMLGG, do.label = T, group.by = "ident", 
                 do.return = TRUE, no.legend = F, 
                 colors.use = ExtractMetaColor(GBMLGG),
                 pt.size = 2,label.size = 5,force = 2)+
        ggtitle("TCGA GBM and LGG bulk RNA-seq")+
        theme(text = element_text(size=10),							
              plot.title = element_text(hjust = 0.5,size = 18, face = "bold")) 

jpeg(paste0(path,"PlotTsne_GBMLGG.jpeg"), units="in", width=10, height=7,
     res=600)
print(p3)
dev.off()

# heatmap.2 =================
TCGA_GBMLGG_markers <- FindAllMarkers.UMI(GBMLGG,logfc.threshold = 0.1, only.pos = T,
                                            test.use = "negbinom")

y = GBMLGG@scale.data[TCGA_GBMLGG_markers$gene,]
## Column clustering (adjust here distance/linkage methods to what you need!)
hc <- hclust(dist(1-cor(as.matrix(y), method="spearman")), method="complete")
cc = GBMLGG@meta.data[hc$labels,"subtype.colors"] %>% as.vector

jpeg(paste0(path,"/Heatmap2_TCGA_GBMLGG.jpeg"), units="in", width=10, height=7,res=600)
heatmap.2(as.matrix(y),
          Colv = as.dendrogram(hc),Rowv= FALSE,
          ColSideColors = cc, trace ="none",labCol = FALSE,dendrogram = "column",
          key.xlab = "scale log RSEM",
          cexRow = 0.01,
          margins = c(5,5),
          breaks = seq(-3,3,length.out = 101),
          #scale = "row",
          col = bluered,
          main = paste("TCGA GBM and LGG bulk RNA-seq"))
par(lend = 1)           # square line ends for the color legend
legend(-0.02, 0.8,       # location of the legend on the heatmap plot
       legend = GBMLGG@meta.data[hc$labels,"subtype"] %>% as.vector %>% unique, # category labels
       col = GBMLGG@meta.data[hc$labels,"subtype.colors"] %>% as.vector %>% unique,  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
dev.off()


