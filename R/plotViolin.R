#' Plot violin plot.
#'
#' Plot violin plot of gene expression in the populations and save to file.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param geneList character vector; genes to show.
#' @param column character; column in the metadata table to annotate the figure.
#' @param column.color character; corresponding color column for the annotation.
#' @param cols numeric; the number of columns in the output figure.
#' @param order character vector; ordered populations. The default is ascending
#' order.
#' @param log boolean; whether to use data in logarithmic scale.
#'
#' @return ggplot2; plot.
#'
#' @export
plotViolin <- function(object, assay.type='RNA',
                       geneList, column, column.color, cols,
                       order = NULL, log = TRUE){
  plots <- list()
  sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
  expDat <- slot(object, 'assays')[[assay.type]]@ndata
  T0 <- expDat
  allgenes <- data.frame()
  for(g in geneList){
    # genes <- as.data.frame(t(T0[g,]))
    genes <- as.data.frame((T0[g,]))
    genes$gene <- g
    genes$conditions <- as.vector(sampTab[, column])
    genes.m <- reshape2::melt(genes, id.var=c('gene',"conditions"))
    allgenes <- rbind(allgenes, genes.m)
  }
  # Log to linear
  if(!log){
    allgenes$value <- 2^(allgenes$value) - 1
  }
  colors <- unique(sampTab[[column.color]])
  names(colors) <- unique(sampTab[[column]])
  if(is.null(order)){
    order <- unique(allgenes$conditions)
    order <- order[order(as.numeric(order), decreasing=FALSE)]
    allgenes$conditions <- factor(allgenes$conditions, levels=order)
  }else{
    allgenes <- allgenes[order(factor(allgenes$conditions, levels=order)),]
    allgenes$conditions <- factor(allgenes$conditions, levels=order)
  }
  p <- ggplot2::ggplot(allgenes, ggplot2::aes(x=conditions, y=value,
                                              fill=conditions)) +
    ggplot2::geom_violin(scale="width") +
    ggplot2::stat_summary(fun=mean, geom='point', size=0.5) +
    ggplot2::theme_bw() + ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme(legend.position="none") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                   axis.text.y = ggplot2::element_text(size=4),
                   strip.text.y = ggplot2::element_text(angle=180),
                   panel.border = ggplot2::element_rect(fill = NA,
                                                        colour=ggplot2::alpha('black', 0.75),
                                                        size=1),
                   strip.background = ggplot2::element_rect(colour="white",
                                                            fill="white"),
                   panel.spacing.x = ggplot2::unit(0.5, "lines"),
                   panel.spacing.y = ggplot2::unit(0, "lines")) +
    ggplot2::scale_y_continuous(position = "right") +
    ggplot2::scale_fill_manual("", values=colors) +
    ggplot2::facet_wrap(~gene, ncol = cols, strip.position = "left",
                        scales = "free_y")
  return(p)
}
