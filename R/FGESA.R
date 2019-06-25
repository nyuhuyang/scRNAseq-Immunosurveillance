library(tidyverse)
library(fgsea)

# Basic function to convert mouse to human gene names
MouseGene2HumanGene <- function(x){
        
        require("biomaRt")
        human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        
        genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
        rm = duplicated(genesV2[,1])
        genesV2 = genesV2[!rm,]
        colnames(genesV2) = c("gene","GENE")
        return(genesV2)
}


res <- read_csv("output/20190419/NAM_vs_Control/All_NAM vs. All_Control.csv")
res
table(res$cluster)
Control <- res[res$cluster %in% "Control",]
NAM <- res[res$cluster %in% "NAM",]


mouse2human <- MouseGene2HumanGene(res$gene)
Control <- inner_join(Control, mouse2human, by="gene")
NAM <- inner_join(NAM, mouse2human, by="gene")

Control2 <- Control %>% 
        dplyr::select(GENE, avg_logFC) %>% 
        na.omit() %>% 
        distinct() %>% 
        group_by(GENE) %>% 
        summarize(stat=mean(avg_logFC))
NAM2 <- NAM %>% 
        dplyr::select(GENE, avg_logFC) %>% 
        na.omit() %>% 
        distinct() %>% 
        group_by(GENE) %>% 
        summarize(stat=mean(avg_logFC))
Control_ranks <- deframe(Control2)
NAM_ranks <- deframe(NAM2)


# Load the pathways into a named list
pathways.hallmark <- gmtPathways("data/msigdb/h.all.v6.2.symbols.gmt")
# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark %>% 
        head() %>% 
        lapply(head)
# Now, run the fgsea algorithm with 1000 permutations:
fgsea_Control <- fgsea(pathways=pathways.hallmark, stats=Control_ranks, nperm=1000)
fgsea_NAM <- fgsea(pathways=pathways.hallmark, stats=NAM_ranks, nperm=1000)

# Tidy the results:
fgsea_Control_Tidy <- fgsea_Control %>%
        as_tibble() %>%
        arrange(desc(NES))

fgsea_NAM_Tidy <- fgsea_NAM %>%
        as_tibble() %>%
        arrange(desc(NES))
# Show in a nice table:
fgsea_Control_Tidy %>% 
        dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
        arrange(padj) %>% 
        DT::datatable()
# Plot the normalized enrichment scores. Color the bar indicating whether or not the pathway was significant:
ggplot(fgsea_Control_Tidy, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj<0.05)) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title="Hallmark pathways NES from GSEA") + 
        theme_minimal()

ggplot(fgsea_NAM_Tidy, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj<0.05)) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title="Hallmark pathways NES from GSEA") + 
        theme_minimal()
