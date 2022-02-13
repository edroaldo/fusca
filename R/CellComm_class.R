#' The CellComm class.
#'
#' An S4 class to represent the CellRouter object and its contents.
#'
#' @details The objetct is used to analyse and find pathways of protein-protein
#' interaction connecting proteins in the surface to the proteins involved in
#' the transcriptional regulation.
#'
#' @slot expdata data frame;
#' @slot ndata data frame;
#' @slot sampTab data frame;
#' @slot rdimension data frame;
#' @slot graph list;
#' @slot signatures list;
#' @slot foldchange list;
#' @slot sources character;
#' @slot targets character;
#' @slot directory list;
#' @slot paths data frame;
#' @slot networks list;
#' @slot genes.trajectory character vector;
#' @slot pathsinfo list;
#' @slot dynamics list;
#' @slot clusters list;
#' @slot correlation list;
#' @slot top.correlations list;
#' @slot davidenrichment list;
#' @slot pathwayenrichment list;
#'
#' @importClassesFrom Matrix dgCMatrix
#'
#' @exportClass CellComm
#' @rdname CellComm
CellComm <- setClass("CellComm", slots=
                       c(expdata="data.frame", ndata="data.frame",
                         sampTab="data.frame", rdimension="data.frame",
                         graph="list", signatures="list", foldchange="list",
                         sources="character", targets="character",
                         directory="list", paths="data.frame", networks="list",
                         genes.trajectory="character", pathsinfo="list",
                         dynamics="list", clusters="list", correlation="list",
                         top.correlations="list", davidenrichment="list",
                         pathwayenrichment="list"))



#' Initialize CellComm object.
#'
#' Initialize CellComm object with no stored data.
#'
#' @param .Object object.
#'
#' @return CellComm object.
#'
#' @importMethodsFrom methods initialize
#'
#' @export
#' @rdname CellComm
#' @aliases CellComm-initialize,initializeCellComm
#' @include CellComm_class.R
setMethod("initialize",
          signature = "CellComm",
          #definition = function(.Object, expdata, annotations){
          definition = function(.Object){
            print("Initializing CellComm object")
            # If necessary, replace by the correct slot.
            #.Object@expdata <- expdata
            #.Object@ndata <- expdata
            #.Object@sampTab <- data.frame(sample_id=colnames(expdata), conditions=annotations)
            #rownames(.Object@sampTab) <- .Object@sampTab$sample_id
            #validObject(.Object)
            return(.Object)
          }
)



# Check CellComm object.
setValidity("CellComm",
            function(object){
              msg <- NULL
              if(!is.data.frame(object@expdata)){
                msg <- c(msg, "expression data must be a data.frame")
              }else if(nrow(object@expdata) < 2){
                msg <- c(msg, "expression data must have more than 2 columns")
              }else if(sum(apply(is.na(object@expdata), 1, sum) > 0)){
                msn <- c(msg, "expression data must not have NAs")
              }#else if(sum(apply(object@expdata, 1, min) < 0)){
              #  msg <- c(msg, "negative values are not allowed in expression data")
              #}
              if(is.null(msg)){
                TRUE
              }else{
                msg
              }
            })

#' Create CellComm object.
#'
#' Create the object passing no parameters.
#'
#' @export
CreateCellComm <- function(){
  object <- new(Class = "CellComm")
  return(object)
}
