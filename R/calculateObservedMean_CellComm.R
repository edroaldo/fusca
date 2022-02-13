#' Mean expression between ligands and receptors.
#'
#' Calculate mean gene expression between ligand and receptor of each
#' interaction in the
#'
#' @param interactions data frame; the interactions as calculated by the
#' population.pairing function.
#' @param mean.expr list; the statistics of the expressed genes as calculated by
#' the computeValue method of CellRouter.
#'
#' @return data frame; the same as returned by the population.pairing function,
#' interactions with celltype1, celltype2, pair, ligand, and receptor, plus the
#' columns of mean ligand expression (mean.ligand), mean receptor expression
#' (mean.receptor), and mean pair expression (mean).
#'
#' @export
calculateObservedMean <- function(interactions, mean.expr){
  interactions2 <- interactions
  interactions2$mean.ligand <- 0
  interactions2$mean.receptor <- 0
  for(i in 1:nrow(interactions2)){
    c1 <- as.vector(interactions2[i, 'celltype1'])
    c2 <- as.vector(interactions2[i, 'celltype2'])
    l <- as.vector(interactions2[i, 'ligand'])
    r <- as.vector(interactions2[i, 'receptor'])
    interactions2[i, 'mean.ligand'] <- mean.expr[[c1]][l, 'p']
    interactions2[i, 'mean.receptor'] <- mean.expr[[c2]][r, 'p']
  }
  interactions2$mean <- as.numeric(
    interactions2$mean.ligand + as.numeric(interactions2$mean.receptor)) / 2

  return(interactions2)
}
