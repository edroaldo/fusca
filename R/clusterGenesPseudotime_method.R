#' Cluster kinetic profiles along trajectories.
#'
#' Cluster kinetic gene expression profiles along each trajectory into the
#' specified number of clusters using the Partitioning Around Medoids method.
#'
#' @param object CellRouter object.
#' @param num.clusters numeric; number of clusters.
#'
#' @return CellRouter object with the clusters slot updated.
#'
#' @export
#' @docType methods
#' @rdname clusterGenesPseudotime-methods
setGeneric("clusterGenesPseudotime", function(object, num.clusters)
  standardGeneric("clusterGenesPseudotime"))
#' @rdname clusterGenesPseudotime-methods
#' @aliases clusterGenesPseudotime
setMethod("clusterGenesPseudotime",
          signature="CellRouter",
          definition=function(object, num.clusters){
            clusters <- list()
            trajectories <- object@dynamics
            for(trajectory in names(trajectories)){
              cat(trajectory, '\n')
              matrix <- trajectories[[trajectory]]
              x <- suppressWarnings(clustergenes(matrix, num.clusters))
              clusters[[trajectory]] <- x
            }
            object@clusters <- clusters
            return(object)
          }
)

#' Helper function of clusterGenesPseudotime.
#'
#' Cluster genes in num.clusters kinetic trends.
#'
#' @param fits matrix; gene expression by pseudotime.
#' @param num.clusters numerical; number of clusters.
clustergenes <- function(fits, num.clusters){
  # Gene expression along 500 columns of pseudotime.
  expr_matrix <- fits
  expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) == 0, ]
  expr_matrix <- t(scale(t(log10(expr_matrix))))
  expr_matrix <- expr_matrix[is.nan(rowSums(expr_matrix)) == FALSE, ]
  expr_matrix[is.na(expr_matrix)] <- 0
  # Get the variance.
  var <- apply(expr_matrix, 1, var)
  var <- var[var > 0]
  n <- as.dist((1 - cor(t(expr_matrix[names(var), ])))/2)
  # Calculate clusters.
  clusters <- cluster::pam(n, num.clusters)
  clusters$exprs <- expr_matrix
  clusters
}
