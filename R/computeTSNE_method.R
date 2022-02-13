#' Perform t-SNE.
#'
#' Dimensionality reduction using t-SNE.
#'
#' @param object CellRouter object.
#' @param num.pcs numeric; number of principal components used for
#' dimensionality reduction using t-SNE.
#' @param perplexity numeric; perplexity.
#' @param max_iter numeric; maximum number of iterations.
#' @param seed numeric; seed.
#'
#' @return CellRouter object with the tsne slot updated.
#'
#' @export
#' @docType methods
#' @rdname computeTSNE-methods
setGeneric("computeTSNE", function(object, num.pcs, perplexity = 30,
                                   max_iter = 2000, seed=7)
  standardGeneric("computeTSNE"))
#' @rdname computeTSNE-methods
#' @aliases computeTSNE
setMethod("computeTSNE",
          signature="CellRouter",
          definition=function(object, num.pcs, perplexity, max_iter, seed){
            if (!is.null(seed)) {
              set.seed(seed = seed)
            }
            pca <- object@pca
            # Computes t-SNE using num.pcs obtained from PCA.
            tsne.done <- Rtsne::Rtsne(pca$cell.embeddings[,1:num.pcs],
                                      perplexity = perplexity,
                                      max_iter = max_iter,
                                      check_duplicates = FALSE)
            m <- tsne.done$Y
            rownames(m) <- rownames(pca$cell.embeddings)
            colnames(m) <- c('tSNE 1', 'tSNE 2')
            object@tsne <- list(cell.embeddings = m)
            object
          }
)
