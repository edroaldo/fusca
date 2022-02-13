#' Create transcription factor network.
#'
#' Create and plot the transcription factor network for the selected interaction
#' and genes. Save plot to file. Load the ggnetwork package before using this
#' function.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param apathways list; the pathways as calculated by the activeSignaling
#' function.
#' @param grn.data list; the gene regulatory network as created by the buildGRN
#' function.
#' @param interactions character; the interactions considered in the analysis.
#' @param genelist character vector; the genes considered in the analysis.
#' @param celltype character; the cell type.
#' @param column character; column corresponding to the cell type.
#'
#' @return list; the igraph transcription factor network and the plot.
#'
#' @export
createTFNetwork <- function(object, assay.type='RNA',
                            apathways, grn.data, interactions, genelist,
                            celltype, column = 'celltype5'){
  graph <- data.frame()
  for(gene in genelist){
    x <- apathways$grns[[interactions]][[gene]]
    xx <- as.data.frame(igraph::get.edgelist(x))
    graph <- rbind(xx, graph)
  }
  graph <- igraph::graph_from_data_frame(graph, directed=FALSE)
  igraph::V(graph)$label <- igraph::V(graph)$name
  geneVals <-
    apply(slot(object, 'assays')[[assay.type]]@ndata[igraph::V(graph)$label,
                                                     slot(object, 'assays')[[assay.type]]@sampTab[[column]] == celltype],
          1, mean);
  igraph::V(graph)[names(geneVals)]$color <- geneVals;
  igraph::V(graph)$type <- 'Target'
  igraph::V(graph)[intersect(igraph::V(graph)$label,
                             grn.data$tfs)]$type <- 'Regulator'
  # library(ggnetwork)
  g1 <- ggplot2::fortify(graph)
  g1$type <- factor(g1$type, levels=c('Target', 'Regulator'))
  set.seed(3)
  p <- ggplot2::ggplot(g1, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggnetwork::geom_edges(color = "grey50") +
    ggnetwork::geom_nodelabel(ggplot2::aes(color = type, label = label,
                                           fontface = "bold")) +
    ggnetwork::theme_blank()
  return(list(df = df, plot = p))
}
