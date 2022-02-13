#' Normalize the data.
#'
#' Normalize the data from CellRouter.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#'
#' @return CellRouter object with the ndata slot updated.
#'
#' @export
#' @docType methods
#' @rdname Normalize-methods
setGeneric("Normalize", function(object, assay.type='RNA')
  standardGeneric("Normalize"))
#' @rdname Normalize-methods
#' @aliases Normalize
setMethod("Normalize",
          signature = "CellRouter",
          definition = function(object, assay.type){
            # ndata has cells as columns.
            # x <- object@ndata
            #x <- as.data.frame(as.matrix((object@ndata)))
            # Normalize by dividing by the total genes in the cell.
            #x <- t(t(x)/apply(x, 2, sum))
            # Applies ln(1+x)
            #x <- log1p(x * 10000)
            # object@ndata <- as(x)
            X <- slot(object, 'assays')[[assay.type]]@ndata
            x <- Matrix::t(Matrix::t(X)/apply1_sp(Matrix::t(X), sum))
            x <- log1p(x * 10000)
            slot(object, 'assays')[[assay.type]]@ndata <- x
            return(object)
          }
)
