#' Dot plot signature score.
#'
#' Present a dot plot of selected genes in all populations in
#' sampTab@@population and save to file.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genes.use character vector; genes to show.
#'
#' @return ggplot2; plot.
#'
#' @export
#' @docType methods
#' @rdname dotplotSignatureScore-methods
setGeneric("dotplotSignatureScore", function(object, assay.type='RNA',
                                             genes.use)
  standardGeneric("dotplotSignatureScore"))
#' @rdname dotplotSignatureScore-methods
#' @aliases dotplotSignatureScore
setMethod("dotplotSignatureScore",
          signature = "CellRouter",
          definition = function(object, assay.type, genes.use){
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            perc <- data.frame(matrix(0, nrow = length(genes.use), ncol = 0))
            exp <- perc
            for(i in unique(sampTab$population)){
              cells.population <- rownames(sampTab[which(sampTab$population == i), ])
              v <- apply(sampTab[cells.population, genes.use], 2, mean)
              exp <- cbind(exp, v)
            }
            colnames(exp) <- unique(sampTab$population)
            rownames(exp) <- sapply(strsplit(rownames(exp), split='__',
                                             fixed = TRUE), function(x){x[1]})
            exp$gene <- rownames(exp)
            exp$gene <- factor(exp$gene, levels = genes.use)
            exp <- reshape2::melt(exp, id.vars = 'gene')
            g <- ggplot2::ggplot(exp, ggplot2::aes(gene, variable)) +
              ggplot2::geom_point(ggplot2::aes(colour = value, size = value)) +
              ggplot2::theme_bw() + ggplot2::xlab("") + ggplot2::ylab("") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12,
                                                                 angle = 45,
                                                                 hjust = 1),
                             panel.grid.minor = ggplot2::element_blank(),
                             panel.grid.major = ggplot2::element_blank(),
                             legend.spacing.y = ggplot2::unit(0, "cm"),
                             panel.border =
                               ggplot2::element_rect(fill = NA,
                                                     colour = ggplot2::alpha('black', 1),
                                                     size = 1)) +
              ggplot2::scale_colour_gradientn("value",
                                              colours = c("midnightblue",
                                                          "dodgerblue3",
                                                          "goldenrod1",
                                                          "darkorange2")) +
              ggplot2::guides(fill = ggplot2::guide_legend(),
                              size = ggplot2::guide_legend(),
                              col = ggplot2::guide_legend(direction = "vertical",
                                                          keywidth = 0.75,
                                                          keyheight = 0.75,
                                                          override.aes = list(size = 3)))
            return(g)
          }
)
