#' Plot Enrichment Result
#'
#' Create enrichment results plot.
#'
#' @param enrichments list; enrichments.
#' @param annotation character; annotation of the enrichment.
#' @param num.pathways numeric; the number of pathways considered in the
#' enrichments.
#' @param order character vector; ordered populations.
#'
#' @return ggplot2 graph.
#'
#' @export
plotEnrichR <- function(enrichments, annotation = 'Reactome_2016',
                        num.pathways = 10, order){
  database <- list()
  for(r in names(enrichments)){
    x <- enrichments[[r]][[annotation]][1:num.pathways, ]
    x$timepoint <- r
    database[[r]] <- x
  }
  database <- do.call(rbind, database)
  df <- database
  df$pvalue <- -log10(df$Adjusted.P.value)
  # Only top GO terms.
  goterms <- as.vector(df$Term)
  clusters <- unique(df$timepoint)
  shared.combined <- data.frame(matrix(NA, nrow = length(clusters),
                                       ncol = length(goterms)))
  rownames(shared.combined) <- clusters
  colnames(shared.combined) <- goterms
  for(c in clusters){
    for(r in goterms){
      pval <- df[which(as.vector(df$Term) == r &
                         as.vector(df$timepoint) == c), ]
      if(nrow(pval) == 1){
        shared.combined[as.character(pval$timepoint),
                        as.character(pval$Term)] <- pval[, 'pvalue']
      }else{
        shared.combined[as.character(pval$timepoint),
                        as.character(pval$Term)] <- NA
      }
    }
  }
  shared.combined$Cluster <- rownames(shared.combined)
  shared.m <- reshape2::melt(shared.combined, id.vars = c('Cluster'))
  colnames(shared.m) <- c('Cluster', 'Description', 'pvalue')
  shared.m$Cluster <- factor(shared.m$Cluster, levels = order)
  shared.m$Description <- gsub("_.*", "", as.character(shared.m$Description))
  shared.m$Description <- factor(shared.m$Description,
                                 levels = unique(shared.m$Description))
  # Create plot.
  g <- ggplot2::ggplot(shared.m, ggplot2::aes(Cluster, Description)) +
    ggplot2::geom_point(ggplot2::aes(colour=pvalue, size=pvalue)) +
    ggplot2::theme_bw() + ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, angle=45, hjust=1),
                   panel.border = ggplot2::element_rect(fill = NA,
                                                        colour = ggplot2::alpha('black', 1),
                                                        size = 1)) +
    ggplot2::guides(color = ggplot2::guide_legend(title="-log10\n(p-value)"),
                    size = ggplot2::guide_legend(title="-log10\n(p-value)")) +
    ggplot2::scale_colour_gradientn("pvalue", colours = c("midnightblue",
                                                          "dodgerblue3",
                                                          "goldenrod1",
                                                          "darkorange2"))
  return(g)
}
