#' Plot reduced dimension space.
#'
#' Plot reduced dimension space using pca, tsne, diffusion components, UMAP or
#' custom.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param reduction.type character; the dimension reduction space to be used:
#' pca, tsne, DC of diffusion components, umap, or custom.
#' @param dims.use numeric vector; dimensions to use.
#' @param annotation character; column in the metadata table to annotate the figure.
#' @param annotation.color character; corresponding color column for the annotation.
#' @param showlabels logical; show labels in the plot.
#' @param dotsize numeric; dot size.
#' @param labelsize numeric; label size.
#' @param convex boolean; whether to connect the dots.
#' @param percentage numeric; proportion of data in the panel (between 0 and 1).
#' @param alpha numeric; transparency (between 0 and 1).
#'
#' @return ggplot2; plot.
#'
#' @export
plotReducedDimension <- function(object, assay.type='RNA',
                                 reduction.type = c('tsne', 'pca', 'DC',
                                                    'umap', 'custom'),
                                 dims.use = c(1, 2), annotation,
                                 annotation.color, showlabels,
                                 dotsize = 1,
                                 labelsize = 3, convex = FALSE,
                                 percentage = 0.80, alpha = 0.1){
  reduction.type <- match.arg(reduction.type)
  # Get the dimensions from the reduction type slot, cells as rows.
  matrix <- slot(object, reduction.type)$cell.embeddings[ , dims.use]
  scores <- as.data.frame(matrix)
  # Uses the rows from the metadata information.
  scores <- scores[rownames(slot(object, 'assays')[[assay.type]]@sampTab), ]
  colnames(scores) <- c('x', 'y')
  # Factors the information from the annotation at the metadata table.
  scores$group <- factor(as.vector(
    slot(object, 'assays')[[assay.type]]@sampTab[[annotation]]),
                         levels = unique(as.vector(
                           slot(object, 'assays')[[assay.type]]@sampTab[[annotation]])))
  # Select colors.
  colors <- unique(slot(object, 'assays')[[assay.type]]@sampTab[[annotation.color]])
  names(colors) <- unique(as.vector(slot(object, 'assays')[[assay.type]]@sampTab[[annotation]]))
  # Create plot
  g1 <- ggplot2::ggplot(scores, ggplot2::aes(x = x, y = y, colour = group)) +
    ggplot2::theme(legend.position = 'right') + ggplot2::geom_point(size = dotsize)
  if(showlabels == TRUE){
    # Gets the median of each group.
    centers <- scores %>% dplyr::group_by(group) %>%
      dplyr::summarize(x = median(x = x), y = median(x = y))
    # Adding text for the top 20 genes.
    g1 <- g1 +
      ggplot2::geom_point(data = centers, mapping = ggplot2::aes(x = x, y = y),
                          size = 0, alpha = 0) +
      ggrepel::geom_text_repel(data = centers,
                               mapping = ggplot2::aes(label = group),
                               size = labelsize, colour = 'black')
  }
  if(reduction.type == 'tsne'){
    xlab <- paste('tSNE ', dims.use[1], sep=' ')
    ylab <- paste('tSNE ', dims.use[2], sep=' ')
  }else if(reduction.type == 'pca'){
    xlab <- paste('PC', dims.use[1], sep='')
    ylab <- paste('PC', dims.use[2], sep='')
  }else if(reduction.type == 'dc'){
    xlab <- paste('DC', dims.use[1], sep='')
    ylab <- paste('DC', dims.use[2], sep='')
  }else if(reduction.type == 'umap'){
    xlab <- paste('UMAP', dims.use[1], sep='')
    ylab <- paste('UMAP', dims.use[2], sep='')
  }else{
    xlab <- 'Dim 1'
    ylab <- 'Dim 2'
  }
  # Add theme, labels...
  g1 <- g1 + ggplot2::theme_bw() +
    ggplot2::ggtitle('')  +  ggplot2::scale_color_manual("", values = colors) +
    ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, colour = "white"),
          strip.background = ggplot2::element_blank()) +
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank()) +
    ggplot2::guides(
      col = ggplot2::guide_legend(direction = "vertical", keywidth = 0.75,
                                  keyheight = 0.85, override.aes = list(size = 3)))
  # If convex.
  if(convex){
    # Prop is the proportion of data in the pannel.
    g1 <- g1 + stat_bag(prop = percentage, alpha = alpha)
  }
  return(g1)
}



#' Helper function of plotReducedDimension.
#'
#' Function for convex object.
#'
#' @inheritParams ggplot2::stat_identity
#' @param prop numeric; proportion of all the points to be included in the bag
#' (default is 0.5).
#' @param na.rm boolean; whether to remove NAs.
#' @param alpha numeric; transparency (between 0 and 1).
#' @return ggplot2 layer.
stat_bag <- function(mapping = NULL, data = NULL, geom = "polygon",
                     position = "identity", na.rm = FALSE, show.legend = NA,
                     inherit.aes = TRUE, prop = 0.5, alpha = 0.3, ...) {
  ggplot2::layer(
    stat = StatBag, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, prop = prop, alpha = alpha, ...)
  )
}


#' Helper function of stat_bag.
#'
#' Creates ggproto object. Originally from aplpack package, plotting functions
#' removed.
#'
#' @return ggproto object.
StatBag <-
  ggplot2::ggproto("Statbag",
                   ggplot2::Stat,
                   compute_group = function(data, scales, prop = 0.5) {
                     plothulls_ <-
                       function(x, y, fraction, n.hull = 1, col.hull, lty.hull,
                                lwd.hull, ensity=0, ...){
                         # Function for data peeling:
                         # x,y : data
                         # fraction.in.inner.hull : max percentage of points within the hull to be drawn
                         # n.hull : number of hulls to be plotted (if there is no fractiion argument)
                         # col.hull, lty.hull, lwd.hull : style of hull line
                         # plotting bits have been removed, BM 160321
                         # pw 130524
                         if(ncol(x) == 2){ y <- x[, 2]; x <- x[, 1] }
                         n <- length(x)
                         if(!missing(fraction)) { # Find special hull
                           n.hull <- 1
                           if(missing(col.hull)) col.hull <- 1
                           if(missing(lty.hull)) lty.hull <- 1
                           if(missing(lwd.hull)) lwd.hull <- 1
                           x.old <- x; y.old <- y
                           idx <- chull(x, y); x.hull <- x[idx]; y.hull <- y[idx]
                           for(i in 1:(length(x)/3)){
                             x <- x[-idx]; y <- y[-idx]
                             if((length(x)/n) < fraction){
                               return(cbind(x.hull, y.hull))
                             }
                             idx <- chull(x, y); x.hull <- x[idx]; y.hull <- y[idx];
                           }
                         }
                         if(missing(col.hull)) col.hull <- 1:n.hull
                         if(length(col.hull)) col.hull <- rep(col.hull, n.hull)
                         if(missing(lty.hull)) lty.hull <- 1:n.hull
                         if(length(lty.hull)) lty.hull <- rep(lty.hull, n.hull)
                         if(missing(lwd.hull)) lwd.hull <- 1
                         if(length(lwd.hull)) lwd.hull <- rep(lwd.hull, n.hull)
                         result <- NULL
                         for(i in 1:n.hull){
                           idx <- chull(x, y); x.hull <- x[idx]; y.hull <- y[idx]
                           result <- c(result, list(cbind(x.hull, y.hull)))
                           x <- x[-idx]; y <- y[-idx]
                           if(0 == length(x)) return(result)
                         }
                         result
                       } # End of definition of plothulls
                       # Prepare data to go into function below
                       the_matrix <- matrix(data = c(data$x, data$y), ncol = 2)
                       # Get data out of function as df with names
                       setNames(data.frame(plothulls_(the_matrix, fraction = prop)),
                                nm = c("x", "y"))
                       # How can we get the hull and loop vertices passed on also?
                     },
                     required_aes = c("x", "y")
)
