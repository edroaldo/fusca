#' Plot kinetic trends for genes in transcriptional clusters.
#'
#' Plot the gene dynamics for each cluster in the selected trajectory and save
#' to file.
#'
#' @param object CellRouter object.
#' @param p character; selected trajectory.
#'
#' @return list; ggplot2 plots.
#'
#' @export
#' @docType methods
#' @rdname plotClusters-methods
setGeneric('plotClusters', function(object, p)
    standardGeneric('plotClusters'))
#' @rdname plotClusters-methods
#' @aliases plotClusters
setMethod('plotClusters',
          signature="CellRouter",
          definition=function(object, p){
            plots <- list()
            clustering <- object@clusters[[p]]$clustering
            for(cl in unique(clustering)){
              print(cl)
              g <- names(clustering[which(clustering == cl)])
              df <- object@dynamics[[p]][g,]
              c <- hclust(dist(df, method = 'euclidean'), method = 'ward.D')
              df <- df[c$order,]
              df <- center_with_threshold(df, 2)
              matrix <- as.data.frame(df)
              matrix$gene <- rownames(df)
              matrix.m <- reshape2::melt(matrix, id.var = "gene")
              matrix.m$gene <- factor(rownames(matrix),
                                      levels=rev(rownames(matrix)))
              g1 <- ggplot2::ggplot(matrix.m, ggplot2::aes(variable, gene)) +
                ggplot2::geom_tile(ggplot2::aes(fill = value)) +
                ggplot2::scale_fill_gradientn("zscore",
                                              colours = c("midnightblue",
                                                          "dodgerblue3",
                                                          "white",
                                                          "goldenrod1",
                                                          "darkorange2")) +
                ggplot2::theme_bw() +
                ggplot2::xlab("") + ggplot2::ylab(paste(length(g), ' genes')) +
                ggplot2::theme(legend.position="none",
                      axis.title.y = ggplot2::element_text(size = ggplot2::rel(2),
                                                           angle = 90),
                      panel.grid.major = ggplot2::element_blank(),
                      panel.grid.minor = ggplot2::element_blank(),
                      axis.text.x = ggplot2::element_blank(),
                      axis.text.y = ggplot2::element_blank(),
                      axis.ticks = ggplot2::element_blank(),
                      panel.border =
                        ggplot2::element_rect(fill = NA,
                                              colour = ggplot2::alpha('black', 1),
                                              size=1))
              plots[[cl]] <- g1
            }
            return(plots)
          }
)
