#' Plot pair dotplot.
#'
#' Plot the mean expression dotplot for the specified ligand-receptor pairs.
#'
#' @param tmp4 matrix?; the inferred cell-cell interactions and ligand-receptor
#' pairs.
#' @param interactionList character vector; the list of cell-cell interactions
#' to be plotted.
#' @param num.pairs numeric; the number of pairs that will be selected.
#'
#' @return list; the selected interactions matrix and the plot.
#'
#' @export
plotPairDotplot <- function(tmp4, interactionList, num.pairs=10){
  z <- tmp4[which(tmp4$celltypes %in% interactionList),]
  x <- split(z, z$celltypes)
  xx <- lapply(x, function(y){y <- y[order(y$mean, decreasing=TRUE), ]; if(nrow(y) < num.pairs){num.pairs=nrow(y)}; y <- y[1:num.pairs,]; y})
  xxx <- do.call(rbind, xx)
  markers <- xxx
  markers <- markers[order(markers$celltypes, -markers$mean),]
  #markers$population <- factor(markers$celltypes, levels=as.vector(unique(markers$celltypes)))
  markers$population <- factor(markers$celltypes, levels=interactionList)
  markers$gene <- factor(markers$pair, levels=as.vector(unique(markers$pair)))
  markers$pvalue2 <- -log10(markers$pvalue)
  g <- ggplot(markers, aes(gene, population)) + geom_point(aes(colour=mean, size=pvalue2)) +
    #g <- ggplot(markers, aes(gene, population)) + geom_point(size=4,aes(colour=mean)) +
    theme_bw() + xlab("") + ylab("") +
    theme(axis.text.x=element_text(size=12, angle=45, hjust=1),
          axis.text.y=element_text(size=12),
          panel.border=element_rect(fill = NA, colour=alpha('black', 1),size=1)) +
    scale_colour_gradientn("Mean(molecule 1, molecule2)",
                           colours=c("midnightblue","dodgerblue3","goldenrod1",
                                     "darkorange2")) +
    guides(col=guide_legend(direction="vertical", keywidth = 0.75,
                            keyheight = 0.75, override.aes = list(size=3))) +
    coord_flip()
  return(list(markers = markers, plot = g))
}
