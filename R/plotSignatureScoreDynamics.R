#' Plot signature scores along a trajectory.
#'
#' Create trajectory signature scores plot.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param trajectories character vector; selected trajectories.
#' @param geneList character vector; genes to show.
#' @param rescale boolean; whether to rescale curves simultaneously or
#' individually.
#' @param order numeric; 	integer giving the order of the neighborhood.
#'
#' @return ggplot2 graph.
#'
#' @export
#' @docType methods
#' @rdname plotSignatureScoreDynamics-methods
setGeneric("plotSignatureScoreDynamics", function(object, assay.type='RNA',
                                                  trajectories,
                                                  geneList, rescale, order=1)
  standardGeneric("plotSignatureScoreDynamics"))
#' @rdname plotSignatureScoreDynamics-methods
#' @aliases plotSignatureScoreDynamics
setMethod("plotSignatureScoreDynamics",
          signature="CellRouter",
          definition=function(object, assay.type, trajectories, geneList,
                              rescale, order=1){
            plotlist <- list()
            alldfs <- data.frame()
            networks <- object@networks
            expDat <- slot(object, 'assays')[[assay.type]]@sampTab[, geneList]
            for(t in trajectories){
              g <- networks[[t]]
              plots <- list()
              trajectory <- object@pathsinfo$path[[t]]
              x_axis <- 1:length(trajectory)
              mean <- list()
              for(cell in trajectory){
                neighs <-
                  igraph::induced.subgraph(graph=g,
                                           vids=unlist(
                                             igraph::neighborhood(graph=g,
                                                                  order=order,
                                                                  nodes=cell)))
                neigh.names <- igraph::V(neighs)$name
                # mean[[cell]] <- apply(expDat[neigh.names,], 2, mean)
                mean[[cell]] <- apply2_sp(expDat[neigh.names, ], sum)/nrow(expDat[neigh.names, ])
              }
              mean.df <- do.call(cbind, mean)
              for(gene_id in geneList){
                y_axis <- as.numeric(mean.df[gene_id,])
                lo <- loess(y_axis~x_axis)
                xl <- seq(min(x_axis), max(x_axis),
                          (max(x_axis) - min(x_axis))/1000)
                y_axis <- predict(lo,xl)
                if(rescale){
                  y_axis <- rescale(y_axis, newrange = c(0,1))
                }
                df <- data.frame(cells=1:length(y_axis),
                                 Expression=as.numeric(y_axis))
                df$gene <- gene_id
                df$cells <- factor(df$cells, levels=df$cells)
                num_subpops <- length(unique(df$population))
                plots[[gene_id]] <- df
              }
              tables <- do.call(rbind, plots)
              labels <- x <- sapply(strsplit(as.vector(tables$gene), split='__',
                                             fixed=TRUE), function(x){x[1]})
              tables$gene <- labels
              if(!rescale){
                tables$Expression <- rescale(tables$Expression,
                                             newrange = c(0,1))
              }
              tables$trajectory <- t
              alldfs <- rbind(alldfs, tables)
            }
            g1 <- ggplot2::ggplot(alldfs, ggplot2::aes(x=cells, y=Expression,
                                                       group=gene, colour=gene)) +
              ggplot2::theme_bw() + ggplot2::geom_line(size=1) +
              ggplot2::xlab('CellRouter trajectory') +
              ggplot2::guides(col=ggplot2::guide_legend(direction="vertical")) +
              ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             legend.position = "right",
                             panel.border = ggplot2::element_blank(),
                             strip.background = ggplot2::element_rect(colour="white",
                                                                      fill="white")) +
              ggplot2::theme(axis.line.x = ggplot2::element_line(color="black",
                                                                 size = 0.5),
                             axis.line.y = ggplot2::element_line(color="black",
                                                                 size = 0.5)) +
              ggplot2::scale_color_manual("", values =
                                            grDevices::rainbow(length(geneList))) +
              ggplot2::facet_wrap(~trajectory, ncol = columns, scales='free')
            return(g1)
          }
)
