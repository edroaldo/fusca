### How refined celltypes/clusters change with age ###
# Crir coluna flsa

#' Plot spatial clusters
#'
#' Plot clusters in tissue image.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param cluster.column character; the name of the column where the clustering
#' information is stored.
#' @param colors.column character; the name of the column where the color of each
#' cluster is stored.
#' @param library.column character; the name of the column where the library
#' information is stored.
#' @param libraries character; the selected library.
#'
#' @return ggplot2; the plots for all libraries.
#'
#' @export
#' @import ggplot2
plotIntegration <- function(object, assay.type='RNA', sample.name='Sample1',
                            cluster.column, colors.column, library.column,
                            libraries){
  # Select sampTab
  sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
  # Select colors.
  colors <- na.exclude(unique(sampTab[[colors.column]]))
  unique_clusters <- na.exclude(unique(sampTab[[cluster.column]]))
  names(colors) <- unique_clusters
  # Select values.
  dr <- as.data.frame(cellrouter@umap$cell.embeddings)
  dr <- dr[rownames(sampTab),]
  colnames(dr) <- c('UMAP1', 'UMAP2')
  dr$group <- sampTab[[cluster.column]]
  dr$library <- libraries[[1]]
  # Add libraries for plotting.
  for (i in 1:length(libraries)) {
    dr[rownames(sampTab[which(sampTab[[library.column]] == libraries[[i]]),]), 'library'] <- names(libraries)[[i]]
  }
  dr$library <- factor(dr$library, levels=names(libraries))
  # Plot figures.
  p1 <- ggplot(dr, aes(x = UMAP1, y=UMAP2, colour=group))  +
    geom_point(size=0.01) +  ggtitle('') + theme_bw() +
    theme(legend.position = 'right',
          panel.border = element_rect(fill = NA, colour = "white"),
          strip.background = element_blank()) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank()) +
    scale_color_manual("", values=colors) +
    facet_wrap(~library, ncol = 3) +
    guides(col=guide_legend(direction="vertical", keywidth = 0.75,
                            keyheight = 0.85, override.aes = list(size=3)))
  return(p1)
}



