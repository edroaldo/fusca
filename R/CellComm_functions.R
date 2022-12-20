#' Calculate Ligand and Receptors
#'
#' Calculate the ligand and receptor pairs from a graph. Used in the predict
#' interactions function.
#'
#' @param graph igraph; the cell graph.
#'
#' @return dataframe; the dataframe with the cell pairs.
#'
#' @export
calculateLigandReceptor <- function(graph){
  edges <- graph$pairs
  edges <- dplyr::bind_rows(edges)
  df <- data.frame(matrix(0, nrow=length(unique(edges$celltype1)),
                          ncol=length(unique(edges$celltype2))))
  rownames(df) <- unique(as.vector(edges$celltype1))
  colnames(df) <- unique(as.vector(edges$celltype2))
  aux <- list()
  for(cell1 in rownames(df)){
    for(cell2 in colnames(df)){
      c <- edges[which(edges$celltype1 == cell1 & edges$celltype2 == cell2),]
      df[cell1, cell2] <- length(unique(c$pair))
    }
  }
  df
}

#' Predict Cell Interactions
#'
#' Plot graphs of the ligands and receptors, their distances and scores.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param cluster.type character; the column where the clusters are
#' indicated. It could be either 'Cluster' or 'Subcluster'. Used to locate the
#' distance matrix in the CellRouter object.
#' @param graph igraph; the cell graph.
#' @param distance.threshold numeric; the threshold for the distance.
#'
#' @return list of ggplot2 graphs.
#'
#' @export
predictCellInteractions <- function(object, assay.type='RNA',
                                    sample.name='Sample1',
                                    cluster.type="Cluster", graph,
                                    distance.threshold=0.75){
  # library(gridExtra)
  # library(grid)
  ligand.receptor <- calculateLigandReceptor(graph = graph)
  ligand.receptor <- rescale(as.matrix(ligand.receptor), c(0,1))
  if(assay.type == "ST"){#!
    distances <- slot(object, 'assays')[[assay.type]]@image$distances[[cluster.type]][[sample.name]]
    distances <- distances[rownames(ligand.receptor), colnames(ligand.receptor)]
    distances <- rescale(distances, c(0,1))
    distances[which(distances > distance.threshold)] <- 1
    proximity.matrix <- 1 - distances
    interaction.score <- ligand.receptor * proximity.matrix
    distance <- pheatmap::pheatmap(proximity.matrix, border_color = "white",
                                 main="Cluster Centroid Spatial Distance",
                                 silent = T)
    score <- pheatmap::pheatmap(interaction.score, border_color = "white",
                              main="Interaction Score", silent = T)
    lr <- pheatmap::pheatmap(ligand.receptor, border_color = "white",
                           main="Ligand-Receptor interactions", silent = T)
    plots <- list(LR=lr$gtable, distance=distance$gtable, score=score$gtable)
  } else { #!
    lr <- pheatmap::pheatmap(ligand.receptor, border_color = "white",
                             main="Ligand-Receptor interactions", silent = T)
    plots <- lr
  }
  return(plots)
}
