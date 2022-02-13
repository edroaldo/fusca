#' Score gene sets.
#'
#' Description.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param genes.list list; gene lists for which the scores will be calculated.
#' @param bins numeric; number of bins to split expression data.
#' @param genes.combine character; genes to combine. The default is all genes in
#' normalized data.
#'
#' @return CellRouter object with the sampTab slot updated.
#'
#' @export
scoreGeneSets <- function(object, assay.type='RNA', sample.name='Sample1',
                          genes.list, bins=25, genes.combine=NULL){
  if(is.null(genes.combine)){
    genes.combine=rownames(slot(object, 'assays')[[assay.type]]@ndata)
  }
  ctrl.size <- min(unlist(lapply(genes.list, length)))
  genes.list <- lapply(genes.list, function(x){intersect(x, rownames(slot(object, 'assays')[[assay.type]]@ndata))})
  cluster.length <- length(x = genes.list)
  data.avg <- Matrix::rowMeans(slot(object, 'assays')[[assay.type]]@ndata[genes.combine, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = Hmisc::cut2(
    x = data.avg,
    m = round(length(x = data.avg) / bins)
  ))
  names(data.cut) <- names(data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(genes.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(sample(data.cut[which(data.cut == data.cut[genes.use[j]])],
                     size = ctrl.size, replace = FALSE))
      )
    }
  }
  ctrl.use <- lapply(ctrl.use, unique)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(ctrl.use),
                        ncol = ncol(slot(object, 'assays')[[assay.type]]@ndata))
  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(slot(object, 'assays')[[assay.type]]@ndata[genes.use, ])
  }
  genes.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length,
                         ncol = ncol(slot(object, 'assays')[[assay.type]]@ndata))
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    genes.scores[i, ] <- Matrix::colMeans(slot(object, 'assays')[[assay.type]]@ndata[genes.use, ,
                                                       drop = FALSE])
  }
  scores <- genes.scores - ctrl.scores
  rownames(scores) <- names(genes.list)
  scores <- as.data.frame(t(scores))
  rownames(scores) <- colnames(slot(object, 'assays')[[assay.type]]@ndata)
  sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
  sampTab <- cbind(sampTab, scores[rownames(sampTab), , drop=FALSE])
  slot(object, 'assays')[[assay.type]]@sampTab <- sampTab
  # for(col.name in colnames(scores)){
  #   object <- addInfo(object, assay.type=assay.type, sample.name=sample.name,
  #                     metadata=scores, colname=col.name, metadata.column=col.name)
  # }
  object
}
