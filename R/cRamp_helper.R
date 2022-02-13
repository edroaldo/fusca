#' Rescale x from 0 to 1.
#'
#' Helper function. Rescale x from 0 to 1 by subtracting the minimum value and
#' dividing by the range.
#'
#' @param x values to be converted to colors.
range01 <- function(x)(x-min(x))/diff(range(x))



#' Convert x to 4 rainbow colors.
#'
#' Helper function.
#'
#' @param x values to be converted to colors
cRamp <- function(x){
  cols <- colorRamp(rev(rainbow(4)))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}



#' Convert to values to 'white', "goldenrod1", and "darkorange2" colors.
#'
#' Helper function.
#'
#' @param x values to be converted to colors
cRamp2 <- function(x){
  cols <- colorRamp(c('white', "goldenrod1", "darkorange2"))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}



#' Convert values to colors.
#'
#' Helper function of addInfo and findClusters.
#'
#' @param x Values to be converted to columns.
#' @param num_colors Number of colors to interpolate.
cRampClust <- function(x, num_colors){
  # Creates a color pallete using the brewer.pal, then interpolates the values
  # provided to these colors.
  cols <- grDevices::colorRamp(RColorBrewer::brewer.pal(num_colors, "Paired"))(range01(x))
  # Convert the RGB values to HEX colors.
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
  #cols
}
