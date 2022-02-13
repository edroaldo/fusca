#' Perform UMAP.
#'
#' Dimensionality reduction using UMAP provided by the umap package.
#'
#' @param object CellRouter object.
#' @param num.pcs numeric; number of principal components used for
#' dimensionality reduction using UMAP.
#' @param spread numeric; used during automatic estimation of a/b parameters.
#' @param min_dist numeric; determines how close points appear in the final
#' layout.
#' @param n_neighbors numeric; number of nearest neighbors.
#' @param metric  character or function; determines how distances between data
#' points are computed.
#' @param set_op_mix_ratio numeric in range [0,1]; determines who the knn-graph
#' is used to create a fuzzy simplicial graph.
#' @param local_connectivity numeric; used during construction of fuzzy
#' simplicial set.
#' @param repulsion_strength numeric; weighting applied to negative samples in
#' low dimensional embedding optimization.
#' @param negative_sample_rate numeric; determines how many non-neighbor points
#' are used per point and per iteration during layout optimization.
#' @param seed numeric; seed.
#'
#' @return CellRouter object with the umap slot updated.
#'
#' @export
#' @docType methods
#' @rdname computeUMAP-methods
setGeneric("computeUMAP", function(object, num.pcs, spread=1, min_dist=0.3,
                                   n_neighbors=30, metric="cosine",
                                   set_op_mix_ratio=1,
                                   local_connectivity=1,
                                   repulsion_strength=1,
                                   negative_sample_rate=5, seed=1)
  standardGeneric("computeUMAP"))
#' @rdname computeUMAP-methods
#' @aliases computeUMAP
setMethod("computeUMAP",
          signature="CellRouter",
          definition=function(object, num.pcs, spread, min_dist,
                              n_neighbors, metric,
                              set_op_mix_ratio,
                              local_connectivity,
                              repulsion_strength,
                              negative_sample_rate, seed){
            if (!is.null(seed)) {
              set.seed(seed=seed)
            }
            pca <- object@pca
            umap.done <- umap::umap(pca$cell.embeddings[,1:num.pcs],
                                    spread=spread,
                                    min_dist=min_dist,
                                    n_neighbors=n_neighbors, metric=metric,
                                    set_op_mix_ratio=set_op_mix_ratio,
                                    local_connectivity=local_connectivity,
                                    repulsion_strength=repulsion_strength,
                                    negative_sample_rate=negative_sample_rate)
            m <- umap.done$layout
            # rownames(m) <- rownames(pca$cell.embeddings)
            colnames(m) <- c('UMAP 1', 'UMAP 2')
            object@umap <- list(cell.embeddings=m)
            object
          }
)
