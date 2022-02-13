#' Perform principal component analysis.
#'
#' Compute the correlation between cells using principal component analysis and
#' plot the standard deviation graph.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param num.pcs numeric; number of principal components to compute.
#' @param genes.use character vector; genes used in principal component
#' analysis. Default is all genes in normalized data.
#' @param seed numeric; seed.
#'
#' @return CellRouter object with the pca slot updated.
#'
#' @export
#' @docType methods
#' @rdname computePCA-methods
setGeneric("computePCA", function(object, assay.type='RNA',
                                  num.pcs, genes.use=NULL, seed=1)
  standardGeneric("computePCA"))
#' @rdname computePCA-methods
#' @aliases computePCA
setMethod("computePCA",
          signature = "CellRouter",
          definition = function(object, assay.type, num.pcs, genes.use, seed){
            if (!is.null(seed)) {
              set.seed(seed = seed)
            }
            if(is.null(genes.use)){
              genes.use <- rownames(slot(object, 'assays')[[assay.type]]@ndata)
            }
            # Fast and memory-efficient way to compute a partial SVD
            pca <- irlba::irlba(A = Matrix::t(slot(object, 'assays')[[assay.type]]@scale.data[genes.use,]),
                                nv = num.pcs)
            # nv approximate right singular vectors
            gene.loadings <- pca$v
            # nu approximate left singular vectors
            cell.embeddings <- pca$u %*% diag(pca$d)
            sdev <- pca$d/sqrt(max(1, ncol(data) - 1))
            rownames(gene.loadings) <- rownames(slot(object, 'assays')[[assay.type]]@scale.data)
            colnames(gene.loadings) <- paste0('PC', 1:num.pcs)
            rownames(cell.embeddings) <- colnames(slot(object, 'assays')[[assay.type]]@scale.data)
            colnames(cell.embeddings) <- colnames(gene.loadings)
            object@pca <- list(gene.loadings = gene.loadings,
                               cell.embeddings = cell.embeddings, sdev = sdev)
            object@rdimension <- as.data.frame(object@pca$cell.embeddings)
            # Plot graph.
            while (!is.null(dev.list()))  dev.off()
            plot(object@pca$sdev, xlab = 'PC',
                 ylab = 'Standard deviation of PC')
            return(object)
          }
)
