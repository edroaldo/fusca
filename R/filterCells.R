#' Quality control.
#'
#' Filter out cells based on variables provided in the parameter variables.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param variables character list; column names in assay@@sampTab of the
#' variables the cells will be filtered on, such as number of detected genes or
#' mitochondrial content.
#' @param thresholds.low numeric vector; cells with values lower than the ones provided
#' here are filtered out.
#' @param thresholds.high numeric vector; cells with values higher than the ones
#' provided here are filtered out.
#'
#' @return CellRouter object with the rawdata, ndata, and sampTab slots updated.
#'
#' @export
filterCells <- function(object, assay.type='RNA',
                        variables, thresholds.low, thresholds.high){
  sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
  for(v in 1:length(variables)){
    sampTab <- sampTab[which(as.vector(sampTab[, variables[v]]) <
                               thresholds.high[v] &
                               as.vector(sampTab[, variables[v]]) >
                               thresholds.low[v]), ]
  }
  if (nrow(slot(object, 'assays')[[assay.type]]@ndata) > 0){
    slot(object, 'assays')[[assay.type]]@ndata <- slot(object, 'assays')[[assay.type]]@ndata[, rownames(sampTab)]
  }
  slot(object, 'assays')[[assay.type]]@ndata <- slot(object, 'assays')[[assay.type]]@ndata[,rownames(sampTab)]
  slot(object, 'assays')[[assay.type]]@rawdata <- slot(object, 'assays')[[assay.type]]@rawdata[,rownames(sampTab)]
  slot(object, 'assays')[[assay.type]]@sampTab <- sampTab
  object
}
