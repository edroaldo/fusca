#' Create a weighted protein-protein interaction network based on cell-type specific correlation coefficients.
#'
#' Create a graph based on gene expression correlation for the protein-protein
#' interactions for each population. The edges are the intersection of the ones
#' in the protein-protein interaction network and the ones obtained from the
#' correlation of gene expression in the cells, the weights are the correlation
#' value. The edge lists are saved to files.
#'
#' @param ppi igraph graph; the network as calculated by the createPPI function.
#' @param expr dgCMatrix; gene expression in each cell, found in
#' CellRouter@@ndata.
#' @param samples data frame; the metadata information of the cells, found in
#' CellRouter@@sampTab.
#' @param corThr numeric; correlation threshold.
#' @param column character; character; column in the metadata table to specify
#' how the data should be aggregated (e.g. per cell type, cluster or other
#' metadata available).
#' @param dir.prefix character; prefix of the directory where the files will be
#' saved.
#'
#' @return list; the networks (igraph objects) for each cell type/cluster/other
#' metadata.
#'
#' @export
cellnetworksPPI <- function(ppi, expr, samples, corThr = 0, column,
                            dir.prefix='results/paths'){
  #
  # duvida: não usa essa corThr.
  #
  corTables <- list()
  # Run for all types in the population.
  for(celltype in as.vector(unique(samples[[column]]))){
    print(celltype)
    # Select cells in the population.
    tmp <- samples[which(samples[[column]] == celltype), ]
    tmp.expr <- expr[, rownames(tmp)]
    # Calculate variance of gene expression for all cells in the population.
    var <- apply(tmp.expr, 1, var)
    # Select the ones with variance greater than zero.
    tmp.expr <- tmp.expr[which(var > 0), ]
    # Calculate correlation.
    cor <- cor(t(as.data.frame(as.matrix((tmp.expr)))))
    print(range(cor))
    cor.m <- reshape2::melt(cor)
    rownames(cor.m) <- paste(as.vector(cor.m$Var1),
                             as.vector(cor.m$Var2), sep='_')
    # Get edgelist on the network.
    el <- igraph::get.edgelist(ppi)
    rownames(el) <- paste(el[,1], el[,2], sep='_')
    keep <- intersect(rownames(cor.m), rownames(el))
    cor.m2 <- cor.m[keep,]
    cor.m2 <- cor.m2[which(cor.m2$value > corThr), ]
    network <- igraph::graph.data.frame(cor.m2, directed=FALSE);
    igraph::E(network)$weight  <- as.vector(cor.m2$value)
    network <- igraph::simplify(network, remove.multiple = TRUE,
                                remove.loops = TRUE)
    edges <- as.data.frame(igraph::get.edgelist(network))
    edges$weight <- abs(igraph::E(network)$weight)
    edges$corr <- igraph::E(network)$weight
    celltype <- gsub(" ", ".", celltype)
    filename <- paste(dir.prefix, celltype, '.txt',sep='')
    # Input network
    write.table(edges, file=filename, sep='\t', row.names=FALSE,
                col.names = FALSE, quote=FALSE)
    corTables[[celltype]] <- network
    gc(verbose=TRUE)
  }
  return(corTables)
}



#' Create a weighted protein-protein interaction network based on cell-type specific correlation coefficients.
#'
#' Create a graph based on gene expression correlation for the protein-protein
#' interactions for each population. The edges are the intersection of the ones
#' in the protein-protein interaction network and the ones obtained from the
#' correlation of gene expression in the cells, the weights are the correlation
#' value. The edge lists are saved to files.
#' A more efficient version of cellnetworksPPI.
#'
#' @param ppi igraph graph; the network as calculated by the createPPI function.
#' @param expr dgCMatrix; gene expression in each cell, found in
#' CellRouter@@ndata.
#' @param samples data frame; the metadata information of the cells, found in
#' CellRouter@@sampTab.
#' @param corThr1 numeric; values smaller then corThr1 in the correlation
#' matrix will be replaced by 0. Therefore, there will be no connections
#' between the genes with value 0.
#' @param corThr2 numeric; keep gene-gene interactions with correlation higher
#' than the threshold.
#' @param column character; character; column in the metadata table to specify
#' how the data should be aggregated (e.g. per cell type, cluster or other
#' metadata available).
#' @param dir.prefix character; prefix of the directory where the files will be
#' saved.
#'
#' @return list; the networks (igraph objects) for each cell type/cluster/other
#' metadata.
#'
#' @export
cellnetworksPPI2 <- function(ppi, expr, samples, corThr1, corThr2, column,
                             dir.prefix='results/paths'){
  corTables <- list()
  for(celltype in as.vector(unique(samples[[column]]))){
    print(celltype)
    tmp <- samples[which(samples[[column]] == celltype),]
    if(nrow(tmp) < 4) next #it needs to contain at all 3 cells/spots to calculate a correlation.
    tmp.expr <- expr[igraph::V(ppi)$name, rownames(tmp)]
    var <- apply(tmp.expr, 1, var)
    tmp.expr <- tmp.expr[which(var > 0),]
    cor <- cor(t(as.data.frame(as.matrix((tmp.expr)))))
    print(range(cor))
    cor[cor < corThr1] <- 0
    cor_g <- igraph::graph_from_adjacency_matrix(cor, mode='undirected',
                                                 weighted = 'correlation')
    cor.m <- igraph::as_data_frame(cor_g, 'edges')
    rownames(cor.m) <- paste(as.vector(cor.m$from), as.vector(cor.m$to), sep='_')
    #overlap of edges
    el <- igraph::get.edgelist(ppi)
    rownames(el) <- paste(el[,1], el[,2], sep='_')
    keep <- intersect(rownames(cor.m), rownames(el))
    cor.m2 <- cor.m[keep,]
    cor.m2 <- cor.m2[which(cor.m2$correlation > corThr2),]
    if (nrow(cor.m2) == 0) {
        print(paste0("No correlation found for ", celltype, " cell type."))
        gc()
    }
    else {
        network <- igraph::graph.data.frame(cor.m2, directed=FALSE);
        igraph::E(network)$weight  <- as.vector(cor.m2$correlation)
        network <- igraph::simplify(network, remove.multiple = TRUE,
                                    remove.loops = TRUE)
        edges <- as.data.frame(igraph::get.edgelist(network))
        edges$weight <- abs(igraph::E(network)$weight)
        edges$corr <- igraph::E(network)$weight
        celltype <- gsub(" ", ".", celltype)
        filename <- paste(dir.prefix, celltype, '.txt',sep='')
        write.table(edges, file=filename, sep='\t', row.names=FALSE,
                    col.names = FALSE, quote=FALSE) #input network
        gc()
        corTables[[celltype]] <- network
        gc(verbose=TRUE)
    }
  }
  return(corTables)
}
