#' Plot information about genes in selected paths.
#'
#' Plot genes along selected trajectories showing each single-cell in the
#' trajectory and curve fit, and save to file.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param paths character vector; selected trajectories.
#' @param genelist character vector; genes to show.
#'
#' @return list of ggplot2 graphs.
#'
#' @export
#' @docType methods
#' @rdname plotpaths-methods
setGeneric("plotpaths", function(object, assay.type='RNA',
                                 paths, genelist)
  standardGeneric("plotpaths"))
#' @rdname plotpaths-methods
#' @aliases plotpaths
setMethod("plotpaths",
          signature="CellRouter",
          definition=function(object, assay.type,
                              paths, genelist){
            path_distr <- object@pathsinfo
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            cellDistances <- object@graph$similarity_matrix
            for(path in paths){
              print(path)
              plots <- list()
              for(gene_id in genelist){
                if(gene_id %in% rownames(path_distr$distr[[path]])){
                  x_axis <- 1:length(path_distr$distr[[path]][gene_id,])
                  y_axis <- path_distr$distr[[path]][gene_id,]
                  df <- data.frame(cells=x_axis, Expression=as.numeric(y_axis),
                                   population=sampTab[names(y_axis), 'population'])
                  df$cells <- factor(df$cells, levels=df$cells)
                  df$population <- factor(df$population,
                                          levels=unique(df$population))
                  num_subpops <- length(unique(df$population))
                  colors <- sampTab[names(y_axis), 'colors']
                  names(colors) <- as.vector(df$population)
                  # First graph.
                  g1 <- ggplot2::ggplot(df, ggplot2::aes(x=cells, y=Expression,
                                                         group=1,
                                                         colour=population)) +
                    ggplot2::geom_point(size=5) + ggplot2::stat_smooth() +
                    ggplot2::theme_bw() +
                    ggplot2::ggtitle(paste(gene_id, path, sep='--')) +
                    ggplot2::scale_colour_manual(values = colors) +
                    ggplot2::xlab("Trajectory") +
                    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                          panel.grid.minor = ggplot2::element_blank(),
                          axis.text.x = ggplot2::element_blank(),
                          panel.background = ggplot2::element_blank(),
                          plot.background = ggplot2::element_blank(),
                          panel.border = ggplot2::element_rect(fill = NA,
                                                               colour = ggplot2::alpha('black', 1),
                                                               size=1))
                  # Save plot to list.
                  plots[[gene_id]] <- g1
                }else{
                  cat(gene_id, ' is not regulated through trajectory: ',
                      path, '\n')
                }
              }
              if(!is.null(cellDistances)){
                # Path similarity matrix.
                matrix <- as.data.frame(as.matrix(
                  cellDistances[as.vector(path_distr$path[[path]]),
                                as.vector(path_distr$path[[path]])]))
                matrix$cells <- rownames(matrix)
                matrix.m <- reshape2::melt(matrix, id.var="cells")
                matrix.m$cells <- factor(as.vector(path_distr$path[[path]]),
                                         levels=as.vector(path_distr$path[[path]]))
                # Second plot.
                g2 <- ggplot2::ggplot(matrix.m, ggplot2::aes(cells, variable)) +
                  ggplot2::geom_tile(ggplot2::aes(fill = value)) +
                  ggplot2::scale_fill_gradient2(low = "blue", mid = "white",
                                                high = "red") +
                  ggplot2::theme_bw() +
                  ggplot2::xlab("Cells") + ggplot2::ylab("Cells") +
                  ggplot2::theme(legend.position = "none",
                        axis.text.x = ggplot2::element_blank(),
                        axis.text.y = ggplot2::element_blank(),
                        axis.ticks = ggplot2::element_blank(),
                        panel.background = ggplot2::element_blank(),
                        panel.grid.major = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        plot.background = ggplot2::element_blank(),
                        panel.border = ggplot2::element_rect(fill = NA,
                                                             colour = ggplot2::alpha('black', 1),
                                                             size = 1)) +
                  ggplot2::ggtitle(path)
                # Save plot to list.
                plots[[path]] <- g2
              }
              return(plots)
            }
          }
)
