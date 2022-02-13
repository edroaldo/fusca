#' Create source target list for each interaction.
#'
#' For each interaction above the threshold, create a list of receptors,
#' containing cell surface receptor as sources and the transcriptional
#' regulators as targets; a list of ligands, sources as ligands and targets as
#' transcriptional regulators; and a data frame of ligands, receptors, and
#' interaction value.
#'
#' @param matrix data frame; the interaction matrix as calculated by the
#' interactionmatrix function.
#' @param threshold numeric; select the interactions which values are above the
#' threshold.
#' @param tfs character vector; the transcription factors.
#'
#' @return list; the receptor, the ligand, and interaction data frame for each
#' interaction.
#'
#' @export
sourcestargets <- function(matrix, threshold, tfs){
  st <- list()
  for(interaction in colnames(matrix)){
    x <- matrix[, interaction]
    names(x) <- rownames(matrix)
    x1 <- x[which(x > threshold)]
    x <- names(x[which(x > threshold)])
    ls <- unique(sapply(strsplit(x, split='_', fixed=TRUE), function(x) (x[1])))
    rs <- unique(sapply(strsplit(x, split='_', fixed=TRUE), function(x) (x[2])))
    receiver <- unique(sapply(strsplit(interaction, split='.', fixed=TRUE),
                              function(x) (x[2])))
    sender <- unique(sapply(strsplit(interaction, split='.', fixed=TRUE),
                            function(x) (x[1])))

    #if ligand is on the cell surface, also compute flows in sender cells...
    #ls <- intersect(cs$genes, ls)
    ##receptor genes based on each potential cell-cell interaction.
    #Right side interaction is the population expressing the receptors
    ##same weighted network but starting from different receptors depending on the
    ##the cell-cell interaction. The cell type at left send signals processed by
    ##receptors at the cell type at right...
    st[['receptor']][[interaction]] <- list(sources=rs, targets=tfs)
    ##membrane ligands genes based on each potential cell-cell interaction
    st[['ligand']][[interaction]] <- list(sources=ls, targets=tfs)

    ls <- sapply(strsplit(x, split='_', fixed=TRUE), function(x) (x[1]))
    rs <- sapply(strsplit(x, split='_', fixed=TRUE), function(x) (x[2]))
    st[[interaction]] <- data.frame(ligand=ls, receptor=rs, weight=x1)
  }
  return(st)
}
