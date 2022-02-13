#' Find pseudo-time curves.
#'
#' Find pseudo-time curves using the princurve library. Alternative to the
#' findPaths method. Add to CellRouter@@curves a list with the cells, their
#' distance from the origin of the curve, and the population they belong,
#' for each path curve.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param column character; column in the metadata table specifying whether
#' transitions are between clusters or other annotations, such as sorted
#' populations.
#' @param population_sequences list of character vectors; the ordered sequence
#' of the populations that will be used to calculate the curve.
#' @param find_populations boolean; if TRUE, the population sequence will be
#' found by the algorithm given a source and target population, if FALSE,
#' the population sequence will be specified by the user.
#' @param sources character vector; name of the source population.
#' @param targets character vector; name of the target population.
#' @param dims.use numeric vector; the number of dimensions to use.
#' @param reduction.type character; the dimension reduction space to be used:
#' pca, tsne, DC of diffusion components, umap, or custom.
#'
#' @return CellRouter object with the curves slot updated.
#'
#' @export
#' @docType methods
#' @rdname findCurves-methods
setGeneric("findCurves", function(object, assay.type='RNA',
                                  column = 'population',
                                  population_sequences = list(),
                                  find_populations = FALSE,
                                  sources = NULL, targets = NULL,
                                  dims.use, reduction.type)
    standardGeneric("findCurves"))
#' @rdname findCurves-methods
#' @aliases findCurves
setMethod("findCurves",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                column = 'population',
                                population_sequences = list(),
                                find_populations = FALSE,
                                sources = NULL, targets = NULL,
                                dims.use, reduction.type){
            if (find_populations == TRUE){
              population_sequences <- list()
              for (i in 1:length(sources)) {
                population_sequences[[i]] <- find_sequence(object, assay.type,
                                                           column,
                                                           sources[i],
                                                           targets[i])
              }
            }
            #gene_path = c('11','9','2','8','6')
            metadata <- slot(object, 'assays')[[assay.type]]@sampTab
            curves <- list()
            for (population_sequence in population_sequences) {
              name_path <- paste(population_sequence, collapse = '->')
              curves[[name_path]] <- find_curve(object, metadata,
                                                column, population_sequence,
                                                dims.use, reduction.type)
            }
            # Add to CellRouter.
            object@curves <- curves
            return(object)
          }
)



#' Find curve
#'
#' Find pseudo-time curve using the princurve library, given the gene path.
#' Alternative to the findPaths method.
#'
#' @param object CellRouter object.
#' @param column character; column in the metadata table specifying whether
#' transitions are between clusters or other annotations, such as sorted
#' populations.
#' @param metadata data frame; metadata table with cell annotations.
#' @param population_sequence list of character vectors; the ordered sequence
#' of the populations that will be used to calculate the curve.
#' @param dims.use numeric vector; the number of dimensions to use.
#' @param reduction.type character; the dimension reduction space to be used:
#' pca, tsne, DC of diffusion components, umap, or custom.
#'
#' @return data frame; the cells, their distance from the origin of the curve,
#' and the curve they belong.
#' @export
find_curve <- function(object, metadata, column, population_sequence, dims.use,
                       reduction.type){
  keep <- rownames(metadata[
    which(metadata[[column]] %in% population_sequence, ), ])
  reduced_dim_data <- slot(object, reduction.type)$cell.embeddings[keep,
                                                                   1:dims.use]
  # Retorna com a população na mesma ordem que esta no sampTab, crescente.
  fit <- princurve::principal_curve(reduced_dim_data,
                                    smoother = 'smooth_spline')
  # plot(fit)
  # Get the curve's first cell and its cluster.
  first_cell <- names(fit$lambda[which(fit$lambda == 0)])
  first_cluster <- metadata[first_cell, column]
  # Check if the curve is ordered correctly.
  if (first_cluster == population_sequence[1]){
    cell_dist <- fit$lambda[order(fit$lambda)]
  } else {
    cell_dist <- fit$lambda[order(fit$lambda, decreasing = TRUE)]
    cell_dist <- abs(cell_dist - max(cell_dist))
  }
  # Create the list.
  curve_data <- list(cell = names(cell_dist), dist = cell_dist,
                     clust = metadata[names(cell_dist), column])
  return(curve_data)
}



#' Find population sequence
#'
#' Find the populations between the source and target ones.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param column character; column in the metadata table specifying whether
#' transitions are between clusters or other annotations, such as sorted
#' populations.
#' @param source character; name of the source population.
#' @param target character; name of the target population.
#' @return character vector; the populations, including the source and target
#' ones.
find_sequence <- function(object, assay.type = 'RNA', column, source, target){
  metadata <- slot(object, 'assays')[[assay.type]]@sampTab
  # Select the cells in the source and target population.
  source_population <- rownames(metadata[which(metadata[[column]] == source), ])
  target_population <- rownames(metadata[which(metadata[[column]] == target), ])
  cell_graph <- object@graph$network
  # Create subgraph from the cells in the source.
  source_subgraph <- igraph::induced_subgraph(cell_graph, vids = source_population)
  # Select the cell with minimum excentricity (largest path to the other cells).
  source_vertex <- names(which.min(igraph::eccentricity(source_subgraph)))
  target_subgraph <- igraph::induced_subgraph(cell_graph, vids = target_population)
  target_vertex <- names(which.min(igraph::eccentricity(target_subgraph)))
  # Select shortest path between source and target cells.
  shortest_cell_path <- igraph::get.shortest.paths(graph = cell_graph,
                                                   from = source_vertex,
                                                   to = target_vertex)
  # Select the names of the cells in the shortest path.
  cell_path <- names(shortest_cell_path$vpath[[1]])
  # Include source to the path, target is already included.
  cell_path <- append(cell_path, source_vertex, after = 0)
  # Select the populations in the path.
  population_sequence <- unique(metadata[cell_path, column])
  return(population_sequence)
}



#' Verify gene expression.
#'
#' Plot a gene expression across a curve.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param curve_data data frame; data from the curve including cell, dist and
#' clust columns.
#' @param gene string; gene name.
verify_gene_exp <- function(object, assay.type='RNA', curve_data, gene){
  df <- data.frame(
    gene_exp = as.numeric(slot(object, 'assays')[[assay.type]]@ndata[gene, rownames(curve_data)]),
    x = curve_data$dist)
  ggplot(df, aes(x = x, y = gene_exp, group = 1)) +
    geom_point(size=5) + stat_smooth() + theme_bw() + #ggtitle(paste(gene_id, path, sep='--')) +
    scale_colour_manual(values = colors) + xlab("Trajectory") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          panel.border = element_rect(fill = NA,
                                      colour = ggplot2::alpha('black', 1),
                                      size = 1))
}
