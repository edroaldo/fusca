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
createPPI <- function(expr, species){
  a <- expr
  if (species = "Hs") {
    organism <- 9606
  }
  else if (species = "Mm") {
    organism <- 10090
  }
  else if (species = "Rn") {
    organism <- 10116
  }
  data <- as.data.frame(import_omnipath_interactions(datasets = 'omnipath', entity_types = 'protein', organism = organism))
  cat(dim(ppi))
  cat('\nProtein interaction data loaded!')
  # The columns names are usually fixed.
  ppi <- data[, c('source_genesymbol','target_genesymbol', 'curation_effort')] #geneA and geneB, score
  idmap <- data.frame(idA=ppi$source_genesymbol, idB=ppi$target_genesymbol, score=as.numeric(ppi$curation_effort))
  allgenes <- unique(c(as.vector(idmap$idA), as.vector(idmap$idB)))
  inetwork <- igraph::graph.data.frame(idmap, directed=FALSE);
  inetwork <- igraph::simplify(inetwork, remove.multiple = TRUE,
                               remove.loops = TRUE)
  genesPPI <- intersect(rownames(expr), allgenes)
  remove <- setdiff(igraph::V(inetwork)$name, rownames(expr))
  inetwork <- igraph::delete.vertices(inetwork, remove)
  ctable <- as.data.frame(igraph::get.edgelist(inetwork))
  colnames(ctable) <- c('from', 'to')
  names <- paste(as.vector(ctable$from), as.vector(ctable$to), sep="_")
  rownames(ctable) <- names
  network <- list(network=inetwork, ctable=ctable, table=idmap, genes=allgenes,
                  genesPPI=genesPPI, removed=remove)
  return(network)
}
