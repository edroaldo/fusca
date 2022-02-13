#' Find genes highly correlated with pseudotime.
#'
#' Find genes highly up-correlated or down-correlated genes in each path
#' according to the specified quantiles.
#'
#' @param object CellRouter object
#' @param max.quantile numeric; quantile to select positively correlated genes.
#' @param min.quantile numeric; quantile to select negatively correlated genes.
#'
#' @return CellRouter object with the top.correlations slot updated.
#'
#' @export
#' @docType methods
#' @rdname topGenes-methods
setGeneric("topGenes", function(object, max.quantile, min.quantile)
  standardGeneric("topGenes"))
#' @rdname topGenes-methods
#' @aliases topGenes
setMethod("topGenes",
          signature="CellRouter",
          definition=function(object, max.quantile, min.quantile){
            # Select correlations.
            slopes <- object@correlation
            geneTrends <- list()
            plots <- list()
            for(path in names(slopes)){
              # Select treshold according to the quantile.
              threshold.up <- quantile(slopes[[path]], max.quantile,
                                       na.rm = TRUE)
              threshold.down <- quantile(slopes[[path]], min.quantile,
                                         na.rm = TRUE)
              # Selec the genes according to the treshold.
              x <- slopes[[path]]
              genesUP <- x[which(x > threshold.up)]
              genesDOWN <- x[which(x < threshold.down)]
              geneTrends[['up']][[path]] <- genesUP[order(genesUP,
                                                          decreasing = TRUE)]
              geneTrends[['down']][[path]] <- genesDOWN[order(genesDOWN)]
            }
            object@top.correlations <- geneTrends
            return(object)
          }
)
