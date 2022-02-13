#' Customized coordinates used to visualize data in a reduced dimension space.
#'
#' Attribute a matrix to \code{object@@custom} as a list.
#'
#' @param object CellRouter object.
#' @param matrix matrix; cell embeddings.
#'
#' @return CellRouter object with the custom slot updated.
#'
#' @export
#' @docType methods
#' @rdname customSpace-methods
setGeneric("customSpace", function(object, matrix)
  standardGeneric("customSpace"))
#' @rdname customSpace-methods
#' @aliases customSpace
setMethod("customSpace",
          signature = "CellRouter",
          definition = function(object, matrix){
            object@custom <- list(cell.embeddings = matrix)
            return(object)
          }
)
