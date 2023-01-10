#' Subset
#' Subset of cells based on metadata provided in the parameter cluster
#' 
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param column character; column name in metadata (sampTab) with the target group
#' @param cluster character; specify the target group.
#'
#' @return CellRouter object with selected cells.  
#'
#' @export
subsetCells <- function(object, assay.type='RNA',column='population', cluster){
  sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
  metadata <- sampTab[grep(pattern = cluster, x = sampTab[[column]]),]
  data = slot(object, 'assays')[[assay.type]]@rawdata[,rownames(metadata), drop=FALSE]
  subset <- CreateCellRouter(data, assay.type = assay.type, min.cells = 0, min.genes=0, is.expr = 0)
  slot(subset, 'assays')[[assay.type]]@sampTab <- metadata
  slot(subset, 'assays')[[assay.type]]@ndata <- slot(object, 'assays')[[assay.type]]@ndata[,rownames(sampTab)]
  slot(subset, 'assays')[[assay.type]]@scale.data <- slot(object, 'assays')[[assay.type]]@scale.data[,rownames(sampTab)]
  slot(subset, 'rdimension') = slot(object, 'rdimension')[rownames(metadata),]
  slot(subset, 'pca')$gene.loadings <- slot(object, 'pca')$gene.loadings
  slot(subset, 'pca')$cell.embeddings <- slot(object, 'pca')$cell.embeddings[rownames(metadata),]
  slot(subset, 'pca')$sdev <- slot(object, 'pca')$sdev
  slot(subset, 'custom')$cell.embeddings <- slot(object, 'custom')$cell.embeddings[rownames(metadata),]
  slot(subset, 'umap')$cell.embeddings <- slot(object, 'umap')$cell.embeddings[rownames(metadata),]
  slot(subset, 'tsne')$cell.embeddings <- slot(object, 'tsne')$cell.embeddings[rownames(metadata),]
  if (assay.type == "ST") {
    slot(subset, 'assays')[[assay.type]]@image <- slot(object, 'assays')[[assay.type]]@image
  }
  subset@var.genes <- object@var.genes
  return (subset)
}
