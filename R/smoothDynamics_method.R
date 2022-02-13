#' Smooth dynamics along trajectories.
#'
#' Smooth the dynamics of gene expression by the path pseudotime using
#' polynomial regression. Helpful to cluster complex transcriptional patterns.
#'
#' @param object CellRouter object.
#' @param names character vector; selected trajectories. Default consists of
#' all trajectories.
#'
#' @return CellRouter object with the dynamics slot updated.
#'
#' @export
#' @docType methods
#' @rdname smoothDynamics-methods
setGeneric("smoothDynamics", function(object, names)
  standardGeneric("smoothDynamics"))
#' @rdname smoothDynamics-methods
#' @aliases smoothDynamics
setMethod("smoothDynamics",
          signature="CellRouter",
          definition=function(object, names){
            dynamics <- list()
            expDat <- object@pathsinfo$distr
            for(path in names(object@pathsinfo$distr[names])){
              print(path)
              x_axis <- as.numeric(object@pathsinfo$pseudotime[[path]])
              geneList <- rownames(expDat[[path]])
              # Creates a matrix with the genes in the path as rows and 501
              # columns.
              smoothDynamics <- data.frame(matrix(0, nrow = length(geneList),
                                                  ncol = 501))
              rownames(smoothDynamics) <- geneList
              for(gene_id in geneList){
                y_axis <- as.numeric(expDat[[path]][gene_id, ])
                # Fit a polynomial regression to the gene expression and the
                # path pseudotime.
                lo <- loess(y_axis~x_axis)
                # Divide the pseudotime length by 500 and create a sequence.
                xl <- seq(min(x_axis), max(x_axis), (max(x_axis) - min(x_axis))/500)
                # Calculate the gene expression using the polynomial regression.
                y_axis <- predict(lo, xl)
                smoothDynamics[gene_id, ] <- as.numeric(y_axis)
              }
              dynamics[[path]] <- smoothDynamics
            }
            object@dynamics <- dynamics
            return(object)
          }
)
