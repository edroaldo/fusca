#' Combine multiple plots into a single one.
#'
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot
#' objects). If the layout is something like matrix(c(1,2,3,3), nrow = 2,
#' byrow = TRUE), then plot 1 will go in the upper left, 2 will go in the upper
#' right, and 3 will go all the way across the bottom. Helper function of
#' plotClusterHeatmap and plotClusters.
#'
#' @param ... ggplot object; the plots.
#' @param plotlist list; the plots.
#' @param file character; filename where the plots will be saved.
#' @param cols numeric; number of columns in the generated figure.
#' @param layout matrix; A matrix specifying the layout. If present, 'cols' is
#' ignored.
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  # Make a list from the ... arguments and plotlist
  # plots <- c(list(...), plotlist) #comment this to plot several plots in list.
  plots <- plotlist
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout.
  if (is.null(layout)) {
    # Make the panel.
    # ncol: Number of columns of plots.
    # nrow: Number of rows needed, calculated from # of cols.
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page.
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout),
                                                                 ncol(layout))))
    # Make each plot, in the correct location.
    for (i in 1:numPlots) {
      # Get the i, j matrix positions of the regions that contain this subplot.
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#' Center the data using threshold.
#'
#' Helper function of plotSignaturesHeatmap, plotPathHeatmap, plotClusters...
#'
#' @param center_data data frame; data to be processed.
#' @param threshold numeric; max value present in center_data.
center_with_threshold <- function(center_data, threshold){
  # Center data (automatically ignores zeros).
  center_data <- center_data - Matrix::rowMeans(center_data, na.rm = TRUE)
  # Keep values between threshold and -threshold.
  center_data[center_data > threshold] <- threshold
  center_data[center_data < (-1 * threshold)] <- -1 * threshold
  return(center_data)
}
