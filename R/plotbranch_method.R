#' Plot branch-specific kinetic patterns.
#'
#' Plot z-scores across specific trajectories.
#'
#' @param object CellRouter object.
#' @param direction character; plot genes up or down-regulated along
#' trajectories.
#' @param p1 character; trajectory 1.
#' @param p2 character; trajectory 2.
#'
#' @return list; ggplot2 graphs.
#'
#' @export
#' @docType methods
#' @rdname plotbranch-methods
setGeneric('plotbranch', function(object, direction = c("up", "down"), p1, p2)
  standardGeneric('plotbranch'))
#' @rdname plotbranch-methods
#' @aliases plotbranch
setMethod('plotbranch',
          signature = "CellRouter",
          definition = function(object, direction = c("up", "down"), p1, p2){
            direction <- match.arg(direction)
            # Branch of interest.
            g <- names(object@top.correlations[[direction]][[p1]])
            # Other branch.
            g <- intersect(g, rownames(cellrouter@dynamics[[p2]]))
            # For the main branch.
            df <- object@dynamics[[p1]][g, ]
            c <- hclust(dist(df, method='euclidean'), method='ward.D')
            df <- df[c$order,]
            df <- as.data.frame(t(scale(t(df))))
            plots <- list()
            matrix <- as.data.frame(df)
            matrix$gene <- rownames(df)
            matrix.m <- reshape2::melt(matrix, id.var="gene")
            matrix.m$gene <- factor(rownames(matrix), levels=rev(rownames(matrix)))
            # Create plot.
            g1 <- ggplot2::ggplot(matrix.m, ggplot2::aes(variable, gene)) +
              ggplot2::geom_tile(ggplot2::aes(fill = value)) +
              ggplot2::scale_fill_gradientn("zscore",
                                            colours = c("midnightblue",
                                                        "dodgerblue3",
                                                        "white",
                                                        "goldenrod1",
                                                        "darkorange2")) +
              ggplot2::theme_bw() +
              ggplot2::xlab("CellRouter trajectory") + ggplot2::ylab("") +
              ggplot2::theme(legend.position="right",
                             axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.3),
                                                                  angle = 90),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_blank(),
                             axis.text.y = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             panel.border = ggplot2::element_rect(fill = NA,
                                                                  colour = ggplot2::alpha('black', 1),
                                                                  size=1)) +
              ggplot2::ggtitle(p1)
            plots[[p1]] <- g1
            # For the second branch.
            # Dynamics of genes in p1 in branch p2.
            df <- object@dynamics[[p2]][g, ]
            df <- df[c$order,]
            df <- as.data.frame(t(scale(t(df))))
            matrix <- as.data.frame(df)
            matrix$gene <- rownames(df)
            matrix.m <- reshape2::melt(matrix, id.var = "gene")
            matrix.m$gene <- factor(rownames(matrix),
                                    levels = rev(rownames(matrix)))
            # Create plot.
            g2 <- ggplot2::ggplot(matrix.m, ggplot2::aes(variable, gene)) +
              ggplot2::geom_tile(ggplot2::aes(fill = value)) +
              ggplot2::scale_fill_gradientn("zscore",
                                            colours = c("midnightblue",
                                                        "dodgerblue3",
                                                        "white","goldenrod1",
                                                        "darkorange2")) +
              ggplot2::theme_bw() +
              ggplot2::xlab("CellRouter trajectory") + ggplot2::ylab("") +
              ggplot2::theme(legend.position="right",
                             axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.3),
                                                                  angle = 90),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_blank(),
                             axis.text.y = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             panel.border = ggplot2::element_rect(fill = NA,
                                                                  colour = ggplot2::alpha('black', 1),
                                                                  size=1)) +
              ggplot2::ggtitle(p2)
            plots[[p2]] <- g2
            return(plots)
          }
)
