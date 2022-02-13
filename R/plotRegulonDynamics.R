#' Plot regulon dynamics.
#'
#' Plot derivative and kinetic patterns of predicted regulators of cell-fate
#' transitions and save to file.
#'
#' @param object CellRouter object.
#' @param p character; selected trajectory.
#' @param regulator character; the selected regulator.
#' @param scores list; scores of transcriptional regulators.
#' @param cluster boolean; whether to cluster kinetic patterns.
#'
#' @return list; ggplot2 plots.
#'
#' @export
#' @docType methods
#' @rdname plotRegulonDynamics-methods
setGeneric('plotRegulonDynamics', function(object, p, regulator, scores,
                                           cluster = TRUE)
  standardGeneric('plotRegulonDynamics'))
#' @rdname plotRegulonDynamics-methods
#' @aliases plotRegulonDynamics
setMethod('plotRegulonDynamics',
          signature="CellRouter",
          definition=function(object, p, regulator, scores, cluster = TRUE){
            #
            # duvida: e esse comentario aqui?
            #
            #show only the ones with changes after regualtor changes?
            #or use derivative > 0 instead of max? >0 as a first line...
            #show line where derivative is equal zero...
            genelist <- scores[[p]][['targets']][[regulator]]
            colors <- c("navy","white","orange")
            # Derivative plot.
            matrix <- object@dynamics[[p]][genelist,]
            time <- 1:501
            # Derivative analysis.
            matrix <- as.data.frame(t(apply(matrix, 1,
                                            function(x){diff(x)/diff(time)})))
            colnames(matrix) <- 1:500
            positions <- as.numeric(colnames(matrix)[max.col(matrix,
                                                             ties.method="first")])
            names(positions) <- rownames(matrix)
            m.regulator <- positions[regulator]
            positions <- positions[positions >= m.regulator]
            position <- positions[regulator]
            matrix <- matrix[names(positions),]
            hc <- hclust(dist(matrix), method='ward.D')
            if(cluster == TRUE){
              matrix <- matrix[hc$order,]
            }
            order <- unique(c(regulator, rownames(matrix)))
            matrix2 <- matrix
            paletteLength <- 100
            myColor <- colorRampPalette(c("navy","white","red"))(paletteLength)
            myBreaks <- c(seq(min(matrix), 0,
                              length.out = ceiling(paletteLength/2) + 1),
                          seq(max(matrix)/paletteLength, max(matrix),
                              length.out = floor(paletteLength/2)))
            matrix$cluster <- rownames(matrix)
            matrix.m <- reshape2::melt(matrix, id.var="cluster")
            matrix.m$cluster <- factor(rownames(matrix), levels=rev(order))
            # First plot.
            g2 <- ggplot2::ggplot(matrix.m, ggplot2::aes(variable, cluster)) +
              ggplot2::geom_tile(ggplot2::aes(fill = value)) +
              ggplot2::geom_vline(ggplot2::aes(xintercept = position),
                                  linetype = "dotted") +
              ggplot2::scale_fill_gradient2("Derivative", low = "navy",
                                            high = "red") +
              ggplot2::theme_bw() +
              ggplot2::xlab("CellRouter trajectory") + ggplot2::ylab("") +
              ggplot2::theme(legend.position="bottom",
                             axis.title.y = element_text(size = ggplot2::rel(0.3),
                                                         angle = 90),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             panel.border = ggplot2::element_rect(fill = NA,
                                                                  colour = ggplot2::alpha('black', 1),
                                                                  size = 1)) +
              ggplot2::ggtitle(p)
            # Expression analysis.
            matrix <- object@dynamics[[p]][rownames(matrix),]
            matrix <- as.data.frame(t(apply(matrix, 1,
                                            function(x){rescale(x, c(0,1))})))
            matrix$cluster <- rownames(matrix)
            matrix.m <- reshape2::melt(matrix, id.var="cluster")
            matrix.m$cluster <- factor(rownames(matrix),
                                       levels = rev(rownames(matrix)))
            # Second plot.
            g1 <- ggplot2::ggplot(matrix.m, ggplot2::aes(variable, cluster)) +
              ggplot2::geom_tile(ggplot2::aes(fill = value)) +
              ggplot2::geom_vline(ggplot2::aes(xintercept=position),
                                  linetype = "solid") +
              ggplot2::scale_fill_gradientn("Scaled expression",
                                            colours = colors) +
              ggplot2::theme_bw() +
              ggplot2::xlab("CellRouter trajectory") + ggplot2::ylab("") +

              ggplot2::theme(legend.position="bottom",
                             axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.3),
                                                                  angle = 90),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_blank(),
                             axis.ticks = ggplot2::element_blank(),
                             panel.border = ggplot2::element_rect(fill = NA,
                                                                  colour = ggplot2::alpha('black', 1),
                                                                  size = 1)) +
              ggplot2::ggtitle(p)
            plots <- list(expression = g2, derivative = g1)
            return(plots)
          }
)
