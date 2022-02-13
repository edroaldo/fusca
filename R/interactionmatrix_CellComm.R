#' Compute the cell-cell interaction data frame.
#'
#' Create a data frame with putative cell-cell interactions.
#'
#' @param result data frame; the interactions for each pair, plus the
#' pvalue column, as calculated by the computePvalue function.
#'
#' @return data frame;
#'
#' @export
interactionmatrix <- function(result){
  # Slpit the result in a list with the values for the celltypes.
  system.time(tmp2 <- split(result, result[, c('celltype1', 'celltype2')]))
  # Create a table of pairs for each celltype.
  summary <- lapply(tmp2, function(x){table(as.vector(x$pair))})
  # Select summaries with length greater than 0.
  summary <- summary[lapply(summary, length) > 0]
  # Substitute marks to _.
  names(summary) <- gsub('([.])', '_', names(summary))
  # Select names for the elements in the list.
  allpairs <- unique(as.vector(unlist(lapply(summary, names))))

  matrix <- data.frame(matrix(0, nrow=length(allpairs), ncol=length(summary)))
  rownames(matrix) <- allpairs
  colnames(matrix) <- names(summary)
  for(c in colnames(matrix)){
    for(p in rownames(matrix)){
      matrix[p, c] <- summary[[c]][p]
    }
  }
  matrix[is.na(matrix)] <- 0
  #matrix <- matrix / max(matrix)

  return(matrix)
}
