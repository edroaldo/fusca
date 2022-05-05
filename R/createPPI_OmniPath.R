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
createPPI <- function(filename, expr){
  a <- expr
  data = read.csv(filename, sep='\t', header=FALSE)
  cat(dim(data))
  cat('\nProtein interaction data loaded!')
  # The columns names are usually fixed.
  ppi <- data[, c('V5','V6', 'V15')] #geneA and geneB, score
  ppi <- as.data.frame(apply(ppi, 2, function(x)gsub('\\s+', '',x)))
  geneA <- as.vector(unlist(lapply(strsplit(as.character(ppi[,'V5']),
                                            split="\\:"), "[", 2)))
  geneB <- as.vector(unlist(lapply(strsplit(as.character(ppi[,'V6']),
                                            split="\\:"), "[", 2)))
  score <- as.vector(unlist(lapply(strsplit(as.character(ppi[,'V15']),
                                            split="\\:"), "[", 2)))
  score <- as.numeric(score)
  geneA <- as.vector(unlist(lapply(strsplit(as.character(geneA),
                                            split="\\|"), "[", 1)))
  geneB <- as.vector(unlist(lapply(strsplit(as.character(geneB),
                                            split="\\|"), "[", 1)))
  idmap <- data.frame(idA=geneA, idB=geneB, score=score)
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
