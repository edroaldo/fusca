#' Builds ligand-receptor network prior knowledge for fusca using Omnipath.
#'
#' @return data frame containing ligand and receptor.
#'
#' @export
Omnipath_LR_network <- function(){
  pairs <- nichenet_lr_network(ramilowski = NULL)
  pairs$Pair.Name <- paste(pairs$from, pairs$to, sep = "_")
  return(pairs)
}
