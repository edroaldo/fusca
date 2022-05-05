#' Create protein-protein interactions.
#'
#' Create a protein-protein interaction network for the reconstruction of cell
#' signaling pathways.
#'
#' @param filename character; protein interaction data.
#' @param expr dgCMatrix; gene expression in each cell, found in
#' CellRouter@@ndata.
#'
#' @return list; the network (igraph object), the ctable (data frame) with the
#' edge list, the table (data frame) with the information of the ids of gene A
#' and B and the scores, the genes in total (character vector), genesPPI
#' (character vector) with the genes in the interactions, and the remove
#' (character vector) with the genes not used in the interactions.
#'
#' @export
Omnipath_LR_network <- function(species){
  pairs <- nichenet_lr_network(only_omnipath = TRUE)
  pairs$Pair.Name <- paste(pairs$from, pairs$to, sep = "_")
  return(pairs)
}
