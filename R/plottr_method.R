#' Plot derivative and kinetic patterns of predicted regulators.
#'
#' Plot dynamics and derivative plots of transcriptional regulators of cell-fate
#' transitions based on the GRN score.
#'
#' @param object CellRouter object.
#' @param p character; trajectory.
#' @param scores numeric vector; scores of transcriptional regulators.
#' @param cluster boolean; whether to cluster kinetic patterns.
#'
#' @return list; ggplot2 plots.
#'
#' @export
#' @docType methods
#' @rdname plottr-methods
setGeneric('plottr', function(object, p, scores, cluster = TRUE)
  standardGeneric('plottr'))
#' @rdname plottr-methods
#' @aliases plottr
setMethod('plottr',
          signature="CellRouter",
          definition=function(object, p, scores, cluster = TRUE){
            colors <- c('blue3', 'white','chartreuse3')
            # Derivative plot.
            matrix <- object@dynamics[[p]][names(scores),]
            time <- 1:501
            # Derivative analysis.
            matrix <- as.data.frame(t(apply(matrix, 1,
                                            function(x){diff(x)/diff(time)})))
            colnames(matrix) <- 1:500
            hc <- hclust(dist(matrix), method = 'ward.D')
            if(cluster == TRUE){
              matrix <- matrix[hc$order,]
            }
            matrix2 <- matrix
            paletteLength <- 100
            myColor <-
              grDevices::colorRampPalette(c("navy","white","red"))(paletteLength)
            myBreaks <- c(seq(min(matrix), 0,
                              length.out = ceiling(paletteLength/2) + 1),
                          seq(max(matrix)/paletteLength, max(matrix),
                              length.out = floor(paletteLength/2)))
            matrix$cluster <- rownames(matrix)
            matrix.m <- reshape2::melt(matrix, id.var = "cluster")
            matrix.m$cluster <- factor(rownames(matrix),
                                       levels = rev(rownames(matrix)))
            g2 <- ggplot2::ggplot(matrix.m, ggplot2::aes(variable, cluster)) +
              ggplot2::geom_tile(ggplot2::aes(fill = value)) +
              ggplot2::scale_fill_gradientn("Expression", colours=colors) +
              ggplot2::theme_bw() +
              ggplot2::xlab("CellRouter trajectory") + ggplot2::ylab("") +
              ggplot2::theme(legend.position = "bottom",
                             axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.3),
                                                                  angle = 90),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             panel.border = ggplot2::element_rect(fill = NA,
                                                                  colour = ggplot2::alpha('black', 1),
                                                                  size=1)) +
              ggplot2::ggtitle(p)
            # Expression plot.
            matrix <- object@dynamics[[p]][rownames(matrix),]
            matrix <- as.data.frame(t(apply(matrix, 1,
                                            function(x){scales::rescale(x, c(0, 1))})))
            matrix$cluster <- rownames(matrix)
            matrix.m <- reshape2::melt(matrix, id.var="cluster")
            matrix.m$cluster <- factor(rownames(matrix),
                                       levels = rev(rownames(matrix)))
            g1 <- ggplot2::ggplot(matrix.m, ggplot2::aes(variable, cluster)) +
              ggplot2::geom_tile(ggplot2::aes(fill = value)) +
              ggplot2::scale_fill_gradientn("Expression", colours = colors) +
              ggplot2::theme_bw() +
              ggplot2::xlab("CellRouter trajectory") + ggplot2::ylab("") +
              ggplot2::theme(legend.position="bottom",
                             axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.3),
                                                                  angle = 90),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             panel.border =
                               ggplot2::element_rect(fill = NA,
                                                     colour = ggplot2::alpha('black', 1),
                                                     size = 1)) +
              ggplot2::ggtitle(p)
            plots <- list(expression = g2, derivative = g1)
            return(plots)
          }
)
