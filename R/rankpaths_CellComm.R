#' Rank signaling pathways.
#'
#' Rank signaling pathways based on the number of targets, the rank, and
#' the mean flow.
#'
#' @param summary data frame; the summary of the files obtained in the
#' findpaths.simple function.
#' @param apathways list; the pathways as calculated by the activeSignaling
#' function.
#'
#' @return list; the list of ranked paths for each interaction (npath) and the
#' list with the number of targets in it (n.targets).
#'
#' @export
rankpaths <- function(summary, apathways){
  paths <- summary$paths
  npaths <- list()
  n.targets <- list()
  for(sub in names(apathways)){
    ps <- paths[[sub]]
    targets <- names(apathways[[sub]])
    ntargets <- unlist(lapply(apathways[[sub]], length))
    n.targets[[sub]] <- ntargets
    x <- ps[as.vector(ps$sink) %in% targets,]
    x$ntargets <- 0
    for(nt in names(ntargets)){
      x[which(x$sink == nt), 'ntargets'] <- ntargets[[nt]]
    }
    x$score <- as.numeric(x$ntargets) * as.numeric(x$rank) * as.numeric(x$meanflow)
    x <- x[order(x$score, decreasing = TRUE),]
    #x <- x[order(x$ntargets, decreasing = TRUE),]
    #ps[as.vector(ps$sink) %in% targets, 'ntargets']
    npaths[[sub]] <- x
  }
  list(npaths=npaths, n.targets=n.targets)
}
