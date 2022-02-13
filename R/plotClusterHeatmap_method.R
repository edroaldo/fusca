#' Plot cluster heatmap.
#'
#' Plot mean kinetic profile for each transcriptional cluster and save plots to
#' file.
#'
#' @param object CellRouter object.
#' @param show character vector; selected trajectories.
#'
#' @return list; ggplot2 plots.
#'
#' @export
#' @docType methods
#' @rdname plotClusterHeatmap-methods
setGeneric("plotClusterHeatmap", function(object, show)
    standardGeneric("plotClusterHeatmap"))
#' @rdname plotClusterHeatmap-methods
#' @aliases plotClusterHeatmap
setMethod("plotClusterHeatmap",
          signature="CellRouter",
          definition=function(object, show){
            # Select clusters to show.
            clusters <- object@clusters[show]
            plots <- list()
            for(trajectory in names(clusters)){
              exprs <- clusters[[trajectory]][['exprs']]
              clustering <- clusters[[trajectory]][['clustering']]
              df <- data.frame(matrix(0, nrow = max(clustering), ncol = 501))
              for(cluster in 1:max(clustering)){
                x <- clustering[which(clustering == cluster)]
                xx <- exprs[names(x), ]
                xxx <- apply(xx, 2, mean)
                df[cluster, ] <- scales::rescale(xxx, c(0,1))
              }
              colors <- c("navy", "yellow", "red")
              matrix <- df
              matrix$cluster <- rownames(matrix)
              matrix.m <- reshape2::melt(matrix, id.var = "cluster")
              matrix.m$cluster <- factor(rownames(df), levels=rev(rownames(df)))
              # Create plot.
              g <- ggplot2::ggplot(matrix.m, ggplot2::aes(variable, cluster)) +
                ggplot2::geom_tile(ggplot2::aes(fill = value)) +
                ggplot2::scale_fill_gradientn("", colours=colors) +
                ggplot2::theme_bw() + ggplot2::xlab("CellRouter trajectory") +
                ggplot2::ylab("") +
                ggplot2::theme(legend.position="right", #axis.text.x = element_text(size=rel(1), angle=45, hjust=1),
                      axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.3),
                                                           angle = 90),
                      panel.grid.major = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank(),
                      axis.text.x = ggplot2::element_blank(),
                      axis.ticks = ggplot2::element_blank(),
                      panel.border =
                        ggplot2::element_rect(fill = NA,
                                              colour = ggplot2::alpha('black', 1),
                                              size=1)) +
                ggplot2::ggtitle(trajectory) +
                ggplot2::scale_x_discrete(expand = c(0, 0)) +
                ggplot2::scale_y_discrete(expand = c(0, 0))
              # Add plot to list.
              plots[[trajectory]] <- g
            }
            return(plots)
          }
)
