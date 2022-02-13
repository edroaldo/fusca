#' Plot pathway enrichment analysis results.
#'
#' Plot the cluster profiler and save to file.
#'
#' @param table compareClusterResult; the compareClusterResult of the
#' compareClusterResult class as obtained from the
#' clusterProfiler::compareCluster constructor.
#' @param num.pathways numeric; the maximum number of pathways considered in the
#' analysis.
#' @param number.words numeric; the minimum number of words for the cluster to
#' be considered.
#' @param order character vector; the ordered names of clusters.
#' @param flip boolean; apply coordinate flip (horizontal and vertical) to the
#' plot.
#'
#' @return ggplot2; plot.
#'
#' @export
plotClusterProfiler <- function(table, num.pathways=10, number.words=3, order,
                                flip=TRUE){
  df <- table
  df$pvalue <- -log10(df$p.adjust)
  nwords <- sapply(strsplit(df$Description, " "), length)
  df$nwords <- nwords
  df <- df[which(df$nwords <= number.words), ]

  df1 <- split(df,df$Cluster)
  df1 <- df1[lapply(df1, nrow) > 0]
  df1 <- lapply(df1, function(x){x[order(x$pvalue, decreasing=TRUE),];
    if(nrow(x) < num.pathways){
      num.pathways <- nrow(x)
    }
    x <- x[1:num.pathways,]})
  df <- do.call(rbind, df1)
  rownames(df) <- NULL

  goterms <- as.vector(df$Description) #only top GO terms
  clusters <- unique(df$Cluster)

  shared.combined <- data.frame(matrix(NA, nrow=length(clusters),
                                       ncol=length(goterms)))
  rownames(shared.combined) <- clusters
  colnames(shared.combined) <- goterms
  for(c in clusters){
    for(r in goterms){
      pval <- df[which(as.vector(df$Description) == r &
                         as.vector(df$Cluster) == c), ]
      if(nrow(pval) == 1){
        shared.combined[as.character(pval$Cluster),
                        as.character(pval$Description)] <- pval[, 'pvalue']
      }else{
        shared.combined[as.character(pval$Cluster),
                        as.character(pval$Description)] <- NA
      }
    }
  }

  shared.combined$Cluster <- rownames(shared.combined)
  shared.m <- reshape2::melt(shared.combined, id.vars = c('Cluster'))

  colnames(shared.m) <- c('Cluster', 'Description', 'pvalue')
  shared.m$Cluster <- factor(shared.m$Cluster, levels = order)
  shared.m$Description <- gsub("_.*", "", as.character(shared.m$Description))
  shared.m$Description <- factor(shared.m$Description,
                                 levels=unique(shared.m$Description))

  g <- ggplot2::ggplot(shared.m, ggplot2::aes(Cluster, Description)) +
    ggplot2::geom_point(ggplot2::aes(fill=pvalue, size=pvalue), pch=21,
                        colour='black') +
    ggplot2::theme_bw() + ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(size=12, angle=45, hjust=1),
                   axis.text.y=ggplot2::element_text(size=12, hjust=1),
          panel.border=ggplot2::element_rect(fill = NA,
                                             colour=ggplot2::alpha('black', 1),
                                             size=1)) +
    ggplot2::guides(color = ggplot2::guide_legend(title="-log10\n(p-value)"),
                    size = ggplot2::guide_legend(title="-log10\n(p-value)")) +
    ggplot2::scale_fill_gradientn("-log10(pvalue)",
                                  colours=c("midnightblue", "dodgerblue3",
                                            "goldenrod1","darkorange2"))
  if(flip){
    g <- g + ggplot2::coord_flip()
  }
  return(g)
}
