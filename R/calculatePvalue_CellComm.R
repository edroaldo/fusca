#' Calculate the p-value for the interactions.
#'
#' Calculate the rate of values which mean expression value for each
#' ligand-receptor pair is greater than the mean expression value for each
#' pairwise cell type.
#'
#' @param p list; the interaction data frame for the permutations calculated by
#' the clusterPermutation function.
#' @param nPerm numeric; the number of permutations.
#' @param interactions2 data frame; the interactions as calculated by the
#' calculateObservedMean function, with the mean expression for ligands,
#' receptors and pairs.
#'
#' @return matrix?; the interactions for each pair, plus the pvalue column.
#'
#' @export
calculatePvalue <- function(p, nPerm=100, interactions2){
  pc <- do.call(rbind, p)
  pc$celltypes <- paste(pc$celltype1, pc$celltype2, sep='_')
  pc$split <- paste(pc$pair, pc$celltypes, sep='_')
  pcc <- split(pc, pc$split)

  interactions2$celltypes <- paste(interactions2$celltype1,
                                   interactions2$celltype2, sep='_')
  interactions2$split <- paste(interactions2$pair,
                               interactions2$celltypes, sep='_')
  aux.int <- split(interactions2, interactions2$split)

  tmp <- list()
  tmp2 <- list()
  #pvalue <- list()
  for(pair in names(pcc)) {
    tmp2[[pair]] <- pair
  }
  
  fun_apply <- function(x) {
    pair.observed <- interactions2[which(interactions2$split == x), ]
    pair.mean.observed <- interactions2[which(interactions2$split == x), 'mean']
    num.higher <- pcc[[x]][which(pcc[[x]]$mean >= pair.mean.observed), ]
    pair.observed$pvalue <- (1+nrow(num.higher))/(1+nPerm)
    pair.observed
  }
  
  tmp2 <- lapply(tmp2, fun_apply)
  tmp3 <- do.call(rbind, tmp2)
  tmp3 <- tmp3[order(tmp3$pvalue),]

  return(tmp3)
}
