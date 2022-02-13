#' Plot genes along a trajectory.
#'
#' Plot gene expressions along a trajectory, return ggplot and save to file.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param trajectory character; selected trajectory.
#' @param geneList character vector; gene list.
#'
#' @return ggplot; the graph.
#'
#' @export
#' @docType methods
#' @rdname plottrajectory-methods
setGeneric("plottrajectory", function(object, assay.type='RNA',
                                      trajectory, geneList)
  standardGeneric("plottrajectory"))
#' @rdname plottrajectory-methods
#' @aliases plottrajectory
setMethod("plottrajectory",
          signature="CellRouter",
          definition=function(object, assay.type,
                              trajectory, geneList){
            plots <- list()
            x_axis <- 1:length(trajectory)
            for(gene_id in geneList){
              y_axis <- as.numeric(slot(object, 'assays')[[assay.type]]@ndata[gene_id, trajectory])
              lo <- loess(y_axis~x_axis)
              xl <- seq(min(x_axis), max(x_axis),
                        (max(x_axis) - min(x_axis))/1000)
              y_axis <- predict(lo, xl)
              y_axis <- scales::rescale(y_axis, newrange = c(0,1))
              df <- data.frame(cells = 1:length(y_axis),
                               Expression = as.numeric(y_axis))
              df$gene <- gene_id
              df$cells <- factor(df$cells, levels = df$cells)
              num_subpops <- length(unique(df$population))
              plots[[gene_id]] <- df
            }
            tables <- do.call(rbind, plots)
            labels <- x <- sapply(strsplit(as.vector(tables$gene),
                                           split='__', fixed=TRUE),
                                  function(x){x[1]})
            tables$gene <- labels
            tables$Expression <- scales::rescale(tables$Expression,
                                                 newrange = c(0, 1))
            # Create plot.
            g1 <- ggplot2::ggplot(tables, ggplot2::aes(x = cells,
                                                       y = Expression,
                                                       group = gene,
                                                       colour = gene)) +
              ggplot2::theme_bw() + ggplot2::geom_line(size=1) +
              ggplot2::xlab('CellRouter trajectory') +
              ggplot2::guides(col = ggplot2::guide_legend(direction = "vertical")) +
              ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             legend.position = "right",
                             panel.border = ggplot2::element_blank()) +
              ggplot2::theme(axis.line.x = ggplot2::element_line(color="black",
                                                                 size = 0.5),
                             axis.line.y = ggplot2::element_line(color="black",
                                                                 size = 0.5)) +
              ggplot2::scale_color_manual("", values = grDevices::rainbow(length(geneList)))
            return(g1)
          }
)




#' Plot genes along specified trajectories.
#'
#' Plot gene expressions along different trajectories and save to file.
#'
#' @param object CellRouter object.
#' @param trajectories character vector; selected trajectories.
#' @param geneList character vector; genes to show.
#' @param rescale boolean; whether to rescale curves simultaneously or
#' individually.
#' @param columns numeric; the number of columns in the output figure.
#'
#' @return ggplot2; plot.
#'
#' @export
#' @docType methods
#' @rdname plottrajectories-methods
setGeneric("plottrajectories", function(object, trajectories, geneList, rescale,
                                        columns = 5)
  standardGeneric("plottrajectories"))
#' @rdname plottrajectories-methods
#' @aliases plottrajectories
setMethod("plottrajectories",
          signature="CellRouter",
          definition=function(object, trajectories, geneList, rescale,
                              columns = 5){
            plotlist <- list()
            alldfs <- data.frame()
            for(t in trajectories){
              plots <- list()
              trajectory <- object@pathsinfo$path[[t]]
              x_axis <- 1:length(trajectory)
              geneList2 <- intersect(geneList, rownames(object@pathsinfo$distr[[t]]))
              for(gene_id in geneList2){
                y_axis <- as.numeric(object@pathsinfo$distr[[t]][gene_id, ])
                lo <- loess(y_axis~x_axis)
                xl <- seq(min(x_axis), max(x_axis),
                          (max(x_axis) - min(x_axis))/1000)
                y_axis <- predict(lo,xl)
                if(rescale){
                  y_axis <- scales::rescale(y_axis, newrange = c(0,1))
                }
                df <- data.frame(cells = 1:length(y_axis),
                                 Expression = as.numeric(y_axis))
                df$gene <- gene_id
                df$cells <- factor(df$cells, levels = df$cells)
                num_subpops <- length(unique(df$population))
                plots[[gene_id]] <- df
              }
              tables <- do.call(rbind, plots)
              labels <- x <- sapply(strsplit(as.vector(tables$gene),
                                             split = '__', fixed = TRUE),
                                    function(x){x[1]})
              tables$gene <- labels
              if(!rescale){
                tables$Expression <- scales::rescale(tables$Expression,
                                                     newrange = c(0, 1))
              }
              tables$trajectory <- t
              alldfs <- rbind(alldfs, tables)
            }
            g1 <- ggplot2::ggplot(alldfs, ggplot2::aes(x = cells,
                                                       y = Expression,
                                                       group = gene,
                                                       colour = gene)) +
              ggplot2::theme_bw() + ggplot2::geom_line(size=1) +
              ggplot2::xlab('CellRouter trajectory') +
              ggplot2::guides(col = ggplot2::guide_legend(direction = "vertical")) +
              ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             legend.position = "right",
                             panel.border = ggplot2::element_blank(),
                             strip.background = ggplot2::element_rect(colour = "white",
                                                                      fill = "white")) +
              ggplot2::theme(axis.line.x = ggplot2::element_line(color="black",
                                                                 size = 0.5),
                             axis.line.y = ggplot2::element_line(color="black",
                                                                 size = 0.5)) +
              ggplot2::scale_color_manual("", values = grDevices::rainbow(length(geneList))) +
              ggplot2::facet_wrap(~trajectory, ncol = columns, scales='free')
            return(g1)
          }
)
