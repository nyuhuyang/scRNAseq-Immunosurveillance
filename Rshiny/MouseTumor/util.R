#' a supprting function for SingleFeaturePlot.1 and FeatureHeatmap.1
#' Change ggplot color scale to increase contrast gradient
#' #https://github.com/satijalab/seurat/issues/235
#' @param p ggplot object
#' @param alpha.use Define transparency of points
#' @param gradient.use Change fill and colour gradient values
#' @param scaled.expression.threshold Define lower limit of scaled gene expression level
#' @export p ggplot object
ChangeColorScale <- function(p1, alpha.use = 1,
                             gradient.use = c("yellow", "red"),
                             scaled.expression.threshold = 0.001) {
    # Order data by scaled gene expresion level
    # Compute maximum value in gene expression
    if (length(p1$data$gene)>0){                   # SingleFeaturePlot.1 
        p1$data = p1$data[order(p1$data$gene),] 
        max.scaled.exp <- max(p1$data$gene) 
    }
    
    # Define lower limit of scaled gene expression level
    if (scaled.expression.threshold == 0) {
        scaled.expression.threshold <- min(p1$data$gene)+0.001
    }

    
    # Fill points using the scaled gene expression levels
    p1$layers[[1]]$mapping$fill <- p1$layers[[1]]$mapping$colour
    
    # Define transparency of points
    p1$layers[[1]]$mapping$alpha <- alpha.use
    
    # Change fill and colour gradient values
    p1 = p1 + guides(colour = FALSE)
    p1 = p1 + scale_colour_gradientn(colours = gradient.use, guide = F,
                                     limits = c(scaled.expression.threshold,
                                                max.scaled.exp),
                                     na.value = "grey") +
        scale_fill_gradientn(colours = gradient.use,
                             name = expression(atop(expression)),
                             limits = c(scaled.expression.threshold,
                                        max.scaled.exp),
                             na.value = "grey") +
        scale_alpha_continuous(range = alpha.use, guide = F)
    
    # Return plot
    print(p1)
}


# prepare gene_exp_coordinates
PrepareExpTsne <- function(exp_dat=exp_dat, coord_dat=coord_dat, input=input){
    
    gene_exp <- matrix(exp_dat[(row.names(exp_dat) == input$text),], nrow = 1,
    dimnames = list(input$text,colnames(exp_dat)))
    gene_exp_t=t(gene_exp)
    gene_exp_t=data.frame(gene_exp_t)
    colnames(gene_exp_t)="gene"
    ##merge the gene expression by coordinates table ##
    gene_exp_coordinates=Reduce(function(x, y) merge(x, y,by = c("row.names"), all=FALSE),
                                list(coord_dat,gene_exp_t))
    gene_exp_coordinates$Sample=gene_exp_coordinates$Row.names
    rownames(gene_exp_coordinates) = gene_exp_coordinates$Row.names
    gene_exp_coordinates=gene_exp_coordinates[,c("tSNE_1","tSNE_2","Sample","gene")]
    colnames(gene_exp_coordinates) = c("x","y","Sample","gene")
    ## scaled values ##
    #gene_exp_coordinates$sdt <- rescale(as.numeric(gene_exp_coordinates$gene))  # data scaled [0, 1]
    data.cut <- as.numeric(x = as.factor(x = cut(x = as.numeric(x = gene_exp_coordinates$gene), 
                                                 breaks = 2)))
    gene_exp_coordinates$col <- as.factor(x = data.cut)
    #rbPal <- colorRampPalette(c("grey",input$Color1))
    #gene_exp_coordinates$Col <- rbPal(10)[as.numeric(cut(gene_exp_coordinates$gene,breaks = 10))]
    
    return(gene_exp_coordinates)
}


#' prepare exp and tsne file
PrepareShiny <- function(object, samples, path){
    
    cell.use <- rownames(object@meta.data)[object@meta.data$orig.ident %in%
    c("Normal",samples)]
    subset.object <- SubsetData(object, cells.use = cell.use)
    subset.object %<>% SetAllIdent(id = "orig.ident")
    subset.object@meta.data$orig.ident %>% unique %>% sort %>% print
    
    exp <- lapply(samples, function(sample) {
        single_object <- SubsetData(subset.object, ident.use = sample)
        single_object@data
    })
    names(exp) = samples
    
    #============== tsne csv===============
    tsne <- lapply(samples, function(sample) {
        single_object <- SubsetData(subset.object, ident.use = sample)
        single_object@dr$tsne@cell.embeddings
    })
    names(tsne) = samples
    save(exp,tsne, file = paste0(path, "data/",basename(path),".Rda"))
}


#' select 1/4 of cell from control
ScaleDown <- function(object, control=c("BH","DJ","MD","NZ")){
    
    normal_cells = lapply(control, function(x){
        rownames(object@meta.data)[(object@meta.data$orig.ident %in% x)]
    })
    set.seed(101)
    remove_normal_cells = lapply(normal_cells, function(x) sample(x, size = length(x)*3/4)) %>% unlist
    table(object@cell.names %in% remove_normal_cells)
    cell.use <- object@cell.names[!(object@cell.names %in% remove_normal_cells)]
    object <- SubsetData(object, cells.use = cell.use)
    object@meta.data$orig.ident = gsub(paste(control,collapse = "|"),"Normal",object@meta.data$orig.ident)
    return(object)
}


#' select input datasets
ShinnyInput <- function(exp,tsne,Dataset1){
    l <- list(exp[[Dataset1]], tsne[[Dataset1]])
    names(l) <- c("exp", "tsne")
    return(l)
}

#' prepare ggplot
Shinnyplot <- function(data = gene_exp_coordinates, input = input, 
                       cols.use = c("lightgrey","blue")){
    g <- ggplot(data = data, mapping = aes(x = x, y = y))
    g <- g + geom_point(mapping = aes(color = gene), 
                        size = as.numeric(input$dotsize), shape = 20)
    g <- g + labs(x = "tSNE_1", y = "tSNE_2")
    g = g + theme_bw()+  
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank())
    g = g + theme(axis.title.y = element_text(size = rel(1.5), angle = 90))
    g = g + theme(axis.title.x = element_text(size = rel(1.5), angle = 00))
    g = g + theme(axis.text=element_text(size=14,angle=90),
                  axis.title=element_text(size=14,face="bold"),
                  plot.title = element_text(size=20, hjust = 0.5),
                  legend.text = element_text(size = 14, colour = "black"),
                  legend.title = element_text(size = 18, colour = "black"))
    g = g + (ggtitle(paste(input$text,"in" ,input$Dataset1))) + 
        theme(plot.title = element_text(size = rel(2.4),face="bold"))
    return(g)
}
