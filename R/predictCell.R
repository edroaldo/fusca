#' Predict cell cycle phase.
#'
#' Predict cell cycle phase based on scores from the metadata table and add
#' phase column with the results to it. Call the addInfo function to add the
#' information to assay@@sampTab.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param columns character vector; columns to be selected from the metadata
#' table.
#'
#' @return CellRouter object with the sampTab slot updated.
#'
#' @export
predictCellCycle <- function(object, assay.type='RNA', columns){
  cc.scores <- slot(object, 'assays')[[assay.type]]@sampTab[ , columns]
  x <- apply(cc.scores, 1, function(x){
    if(all(x < 0)){
      return('G1')
    }else{
      return(names(x)[which(x == max(x))][1])
    }
  })
  object <- addInfo(object, assay.type, x, colname = 'Phase')
  object
}

#' Predict state based on scored gene lists.
#'
#' Predict cell cycle phase based on scores from the metadata table and add
#' column with the results to it. Call the addInfo function to add the
#' information to assay@@sampTab.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param columns character vector; columns to be selected from the metadata
#' table.
#' @param col.name character; column name to be added to the metadata table
#' after state prediction.
#'
#' @return CellRouter object with the sampTab slot updated.
#'
#' @export
predictState <- function(object, assay.type='RNA', columns, col.name){
  cc.scores <- slot(object, 'assays')[[assay.type]]@sampTab[, columns]
  x <- apply(cc.scores, 1, function(x){
    if(all(x < 0)){
      return('Not assigned')
    }else{
      return(names(x)[which(x == max(x))][1])
    }
  })
  object <- addInfo(object, assay.type, x, colname = col.name)
  object
}
