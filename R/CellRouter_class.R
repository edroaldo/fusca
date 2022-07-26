# Join matrix and dgCMatrix to be used in the CellRouter object.
# The dgTMatrix matrix as added to work with ST data.
setClassUnion(name = 'AnyMatrix', members = c("matrix", "dgCMatrix"))

#' The CellRouter class.
#'
#' An S4 class to represent the CellRouter object and its contents.
#'
#' @details A multifaceted single-cell analysis platform that identifies complex
#' cell-state transition trajectories by using flow networks to explore the
#' subpopulation structure of multi-dimensional, single-cell omics data.
#'
#' @slot assays list; different data used in the analysis.
#' @slot sampTab data frame; store the information of the integrated analysis of
#' multiple assays.
#' @slot rdimension data frame; cell embeddings from PCA.
#' @slot pca list; data from principal component analysis.
#' @slot tsne list; data from t-SNE dimensionality reduction.
#' @slot umap list; data from UMAP dimensionality reduction.
#' @slot dr.custom list; data from custom dimensionality reduction.
#' @slot dc list; data from dc dimensionality reduction.
#' @slot custom list; cell embeddings for custom space.
#' @slot var.genes A character.
#' @slot graph A list.
#' @slot signatures A list.
#' @slot sources A character.
#' @slot targets A character.
#' @slot directory list; directories of files obtained in the findPaths method.
#' @slot paths data frame; paths between populations.
#' @slot networks list; graphs between populations.
#' @slot genes.trajectory character vector; genes considered in the trajectories.
#' @slot pathsinfo list; information about the paths across populations.
#' @slot dynamics list; the smoothed correlations between gene expression and
#' pseudotime in paths between populations.
#' @slot clusters list; the clusters for the selected trajectories.
#' @slot correlation list; correlations between gene expression and pseudotime
#' in paths between populations.
#' @slot top.correlations list; gene trends (up or down) for each path.
#' @slot pathwayenrichment list; pathway enrichment analysis on genes up or
#' down-regulated along trajectories.
#' @slot curves list; the curves for the specified paths.
#'
#' @importClassesFrom Matrix dgCMatrix
#'
#' @exportClass CellRouter
#' @rdname CellRouter
CellRouter <- setClass("CellRouter", slots=
                         c(assays="list", sampTab="data.frame",
                           rdimension="data.frame", pca="list",
                           tsne="list", umap="list", dr.custom="list",
                           dc="list", custom="list", var.genes="character",
                           graph="list", signatures="list", sources="character",
                           targets="character", directory="list",
                           paths="data.frame", networks="list",
                           genes.trajectory="character", pathsinfo="list",
                           dynamics="list", clusters="list", correlation="list",
                           top.correlations="list", pathwayenrichment="list",
                           curves = "list"))



#' CellRouter assay.
#'
#' Assay to create CellRouter.
#'
#' @slot rawdata AnyMatrix; raw data, cells should be represented as columns,
#' genes should be represented as rows.
#' @slot ndata AnyMatrix; original data, cells should be represented as columns,
#' genes should be represented as rows.
#' @slot scale.data AnyMatrix; scaled data.
#' @slot sampTab dataframe; data frame; metadata information with cells names as
#' row names, number of genes as nGene, number of UMI as nUMI, and cell names as
#' conditions. It also includes the clusters and subclusters of the data in the
#' assay.
#' @slot image list; image information of the ST assay.
#'
#' @exportClass FuscaAssay
FuscaAssay <- setClass(
  Class = 'FuscaAssay',
  slots = c(
    rawdata = 'AnyMatrix',
    ndata = 'AnyMatrix',
    scale.data = 'AnyMatrix',
    sampTab = 'data.frame',
    image = 'list'
  )
)



#' Create CellRouter.
#'
#' Create the object and check its parameters.
#'
#' @param rawdata AnyMatrix; raw data, cells should be represented as columns,
#' genes should be represented as rows.
#' @param assay.type character; the type of data provided to the CellRouter
#' object.
#' @param min.genes numeric; minimum of genes a cell should express for it to be
#' considered.
#' @param min.cells numeric; minimum of cells in which a gene should be
#' expressed for it to be considered.
#' @param is.expr numeric; the minimum expression value a gene should have to be
#' considered expressed.
#'
#' @export
CreateCellRouter <- function(rawdata, assay.type='RNA', min.genes, min.cells,
                             is.expr = 0){
  assay.data <- CreateAssay(rawdata, assay.type, min.genes, min.cells, is.expr);
  assay.list <- list(assay.data)
  names(assay.list) <- assay.type
  object <- new(Class = "CellRouter", assays = assay.list)
  # object@sampTab <- assay.data@sampTab
  return(object)
}


# Only works for RNA.
# Mudar para os diferentes tipos, acredito que só funcione pra o RNA agora.
CreateAssay <- function(rawdata, assay.type, min.genes, min.cells, is.expr){
  rownames(rawdata) <- make.unique(rownames(rawdata), sep = ".") #!
  colnames(rawdata) <- gsub('-','.',colnames(rawdata)) #!
  num.genes <- Matrix::colSums(rawdata > is.expr)
  num.mol <- Matrix::colSums(rawdata)
  cells.use <- names(num.genes[which(num.genes > min.genes)])
  #expdat <- rawdata[, cells.use]
  rawdata <- rawdata[, cells.use]
  genes.use <- rownames(rawdata)
  if (min.cells > 0) {
    num.cells <- Matrix::rowSums(rawdata > 0)
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    rawdata <- rawdata[genes.use, ]
  }
  nGene <- num.genes[cells.use]
  nUMI <- num.mol[cells.use]
  sampTab <- data.frame(sample_id=colnames(rawdata), nGene=nGene, nUMI=nUMI,
                        conditions=colnames(rawdata))
  rownames(sampTab) <- sampTab$sample_id
  print(dim(rawdata))
  rawdata <- rawdata[, rownames(sampTab)]
  assay <- new(Class = "FuscaAssay", rawdata = rawdata, ndata = rawdata,
               scale.data = new(Class = "matrix"), sampTab = sampTab)
  return(assay)
}

#' Add assay
#'
#' Add other assay to the CellRouter object.
#'
#' @param object CellRouter object.
#' @param rawdata AnyMatrix; raw data, cells should be represented as columns,
#' genes should be represented as rows.
#' @param assay.type character; the type of data provided to the CellRouter
#' object.
#' @param min.genes numeric; minimum of genes a cell should express for it to be
#' considered.
#' @param min.cells numeric; minimum of cells in which a gene should be
#' expressed for it to be considered.
#' @param is.expr numeric; the minimum expression value a gene should have to be
#' considered expressed.
#'
#' @export
AddAssay <- function(object, rawdata, assay.type, min.genes, min.cells, is.expr){
  assay.data <- CreateAssay(rawdata, assay.type, min.genes, min.cells, is.expr);
  object@assays[[assay.type]] <- assay.data
  return(object)
}






#' Test
#'
#' Test if the function is exported.
#'
#' @param num numeric;
#'
#' @export
testExport <- function(num){
  print('The function is exported: ', num)
}
