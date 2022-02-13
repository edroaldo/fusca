#' Proportion plot.
#'
#' Plot a graph of the cells in two populations.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param condition character; column in the metadata table specifying an
#' annotation, such as sorted populations.
#' @param population character; column in the metadata table specifying another
#' annotation.
#' @param order character vector; ordered populations. The default is ascending
#' order.
#'
#' @return list; the data frame of cells (condition) and their classifications
#' (population), and the plot.
#'
#' @export
plotProportion <- function(object, assay.type='RNA',
                           condition, population, order=NULL){
  samples <- slot(object, 'assays')[[assay.type]]@sampTab
  data2 <- data.frame(cells=samples[[condition]],
                      classification=samples[[population]])
  colors <- as.vector(unique(samples$colors))
  names(colors) <- unique(as.vector(samples$population))
  if(is.null(order)){
    data2$classification <- factor(data2$classification,
                                   levels=unique(as.vector(data2$classification)))
  }else{
    data2 <- data2[order(factor(data2$classification, levels=order)), ]
    data2$classification <- factor(data2$classification, levels=order)
  }
  g <- ggplot2::ggplot(data2, ggplot2::aes(x = classification, fill = cells)) +
    ggplot2::geom_bar(position = "fill", color='black') +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position='right',
                   legend.key.size = ggplot2::unit(0.3, "cm"),
                   legend.text = ggplot2::element_text(size=7)) +
    ggplot2::xlab("") + ggplot2::ylab("Proportion") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, angle=45,
                                                       hjust=1),
                   axis.title.y = ggplot2::element_text(size = ggplot2::rel(1),
                                                        angle = 90)) +
    ggplot2::scale_fill_brewer("", palette = 'Paired') +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0))
  return(list(df = data2, plot = g))
}



#' Proportion plot.
#'
#' Plot a graph of the cells in two populations.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param condition character; column in the metadata table specifying an
#' annotation, such as sorted populations.
#' @param population character; column in the metadata table specifying another
#' annotation.
#' @param color.column character; corresponding color column for the annotation.
#' @param order character vector; ordered populations. The default is ascending
#' order.
#'
#' @return list; the data frame of cells (condition) and their classifications
#' (population), and the plot.
#'
#' @export
plotProportion2 <- function(object, assay.type='RNA',
                            condition, population, color.column, order = NULL){
  samples <- slot(object, 'assays')[[assay.type]]@sampTab
  data2 <- data.frame(cells=samples[[condition]],
                      classification=samples[[population]])
  colors <- as.vector(unique(samples[[color.column]]))
  names(colors) <- unique(as.vector(samples[[condition]]))
  if(is.null(order)){
    data2$classification <- factor(data2$classification,
                                   levels=unique(as.vector(data2$classification)))
  }else{
    data2 <- data2[order(factor(data2$classification, levels=order)),]
    data2$classification <- factor(data2$classification, levels=order)
  }
  g <- ggplot2::ggplot(data2,aes(x = classification, fill = cells)) +
    ggplot2::geom_bar(position = "fill", color='black') +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position='right',
                   legend.key.size = ggplot2::unit(0.4, "cm"),
                   legend.text = ggplot2::element_text(size=7)) +
    ggplot2::xlab("") + ggplot2::ylab("Proportion") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, angle=45,
                                                       hjust=1),
                   axis.title.y = ggplot2::element_text(size = ggplot2::rel(1),
                                                        angle = 90)) +
    ggplot2::scale_fill_manual("", values=colors) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0))
  return(list(df = data2, plot = g))
}
