#' Dimensionality reduction using diffusion components.
#'
#' Performs dimensionality reduction using the DiffusionMap function from the
#' destiny package. The destiny package is presenting errors in the installation,
#' so it is not imported in the description.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genes.use character vector; genes used in diffusion components
#' dimensionality reduction. The default is all genes in normalized data.
#' @param k numeric; parameter k to be used by the DiffusionMap function from
#' destiny pakcage.
#' @param sigma character; parameter sigma to be used by the DiffusionMap
#' function from destiny pakcage.
#' @param seed numeric; random seed.
#'
#' @return CellRouter object with the dc slot updated.
#' @export
computeDC <- function(object, assay.type='RNA',
                      genes.use = NULL, k = 20, sigma = 'local',
                      seed = 1){
  # The destiny package is presenting errors in the installation, so it is not
  # in the namespace.
  # library(destiny) #anyoing error with DLLs all the time...
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  if(is.null(genes.use)){
    genes.use <- object@var.genes
  }
  diff.comp <- destiny::DiffusionMap(as.matrix(t(
    slot(object, 'assays')[[assay.type]]@scale.data[genes.use,])),
                                     k = 20, sigma = 'local')
  dc <- eigenvectors(diff.comp)
  rownames(dc) <- colnames(slot(object, 'assays')[[assay.type]]@scale.data)
  object@dc <- list(cell.embeddings = dc)
  object
}
