#' Plot heatmap with gene signatures.
#'
#' Plot heatmap of gene expression for the specified genes in each population.
#' The plot is not shown in the RStudio Plots display, but it is saved in the
#' specified directory.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param markers tibble; genes preferentially expressed in each column.ann, as
#' obtained in the findSignatures function. For example, in clusters or sorted
#' populations.
#' @param column.ann character; column in the metadata table used to annotate
#' the kNN graph. For example, clusters, sorted cell populations.
#' @param column.color character; column corresponding to the colors.
#' @param num.cells numeric; number of cells to show in the heatmap.
#' @param threshold numeric; threshold used to center the data.
#' @param genes.show character vector;  gene names to show in the heatmap. The
#' default is the top 5 genes according to fold change for each population.
#' @param low character; color for low expression.
#' @param intermediate character; color for intermediate expression.
#' @param high character; color for high expression.
#' @param order character vector; order of population names, the default is
#' alphabetical order.
#' @param fontsize numeric; font size.
#'
#' @return pheatmap; plot.
#'
#' @export
plotSignaturesHeatmap <- function(object, assay.type='RNA',
                                  markers, column.ann, column.color,
                                  num.cells = NULL, threshold = 2,
                                  genes.show = NULL, low = 'purple',
                                  intermediate = 'black', high = 'yellow',
                                  order = NULL, fontsize = 5){
  if(is.null(num.cells)){
    print('here')
    cells.keep <- rownames(slot(object, 'assays')[[assay.type]]@sampTab)
    print(table(slot(object, 'assays')[[assay.type]]@sampTab[[column.ann]]))
  } else {
    # Split the data acording to cell type, creating a list of data from each
    # type.
    cells.use <- split(slot(object, 'assays')[[assay.type]]@sampTab,
                       slot(object, 'assays')[[assay.type]]@sampTab[[column.ann]])
    cells.use <- lapply(cells.use, function(x){
      if(nrow(x) < num.cells){
        cells.use.x <- x[sample(rownames(x), size = nrow(x)), ]
      } else {
        cells.use.xx <- x[sample(rownames(x), size = num.cells), ]
      }
    })
    cells.use.tmp <- do.call(rbind, cells.use)
    cells.keep <- as.vector(cells.use.tmp$sample_id)
  }
  matrix <- center_with_threshold(
    slot(object, 'assays')[[assay.type]]@ndata[, cells.keep], threshold)
  paletteLength <- 100
  myColor <- grDevices::colorRampPalette(c(low, intermediate, high))(paletteLength)
  myBreaks <- c(seq(min(matrix), 0, length.out = ceiling(paletteLength/2) + 1),
                seq(max(matrix)/paletteLength, max(matrix),
                    length.out = floor(paletteLength/2)))
  # Need to make sure there is no duplicated element. This is done at the
  # findSignatures function.
  markers2 <- as.data.frame(markers)
  sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
  sampTab <- sampTab[cells.keep,]
  if(column.ann == 'population') {
    markers2 <- markers2[order(as.numeric(markers2$population)),]
    rownames(markers2) <- as.vector(markers2$gene)
    sampTab <- sampTab[order(as.numeric(sampTab$population)),]
  } else if(!is.null(order)) {
    markers2 <- markers2[order(factor(markers2$population, levels = order)),]
    sampTab <- sampTab[order(factor(sampTab[[column.ann]], levels = order)),]
  } else {
    markers2 <- markers2[order(as.character(markers2$population)),]
    rownames(markers2) <- as.vector(markers2$gene)
    sampTab <- sampTab[order(as.character(sampTab[[column.ann]])),]
  }
  clusters <- as.vector(sampTab[[column.ann]])
  names(clusters) <- rownames(sampTab)
  # Data frame with cells as rows and a column with their population.
  ann_col <- data.frame(population = as.vector(clusters),
                        stringsAsFactors = FALSE)
  rownames(ann_col) <- names(clusters)
  # Data frame with genes as rows and a column with their population.
  ann_row <- data.frame(signature = as.vector(markers2$population),
                        stringsAsFactors = FALSE)
  rownames(ann_row) <- as.vector(markers2$gene)
  if(!is.null(order)){
    ann_col$population <- factor(ann_col$population, levels = order)
    ann_row$signature <- factor(ann_row$signature, levels = order)
  }
  colors <- unique(sampTab[[column.color]])
  names(colors) <- unique(as.vector(sampTab[[column.ann]]))
  color_lists <- list(population = colors, signature = colors)
  index <- getIndexes(ann_col, ann_row,
                      order.columns = unique(ann_col$population),
                      order.rows = unique(ann_row$signature))
  # Select the genes to show and attribute "" to the ohter ones.
  if(is.null(genes.show)){
    genes.show <- markers2 %>% dplyr::group_by(population) %>% dplyr::top_n(5, fc)
    genes.show <- as.vector(genes.show$gene)
    selected <- as.vector(markers2$gene)
    selected[!(selected %in% genes.show)] <- ""
  } else {
    selected <- as.vector(markers2$gene)
    selected[!(selected %in% genes.show)] <- ""
  }
  # Plot in not shown in the RStudio plots, but it is saved in the directory.
  ghm <- pheatmap::pheatmap(matrix[rownames(ann_row), rownames(ann_col)],
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     color = myColor, breaks=myBreaks, fontsize = fontsize,
                     gaps_row = index$rowsep, gaps_col = index$colsep,
                     annotation_col = ann_col, annotation_row = ann_row,
                     annotation_colors = color_lists, labels_row = selected,
                     labels_col = rep("", ncol(matrix)), silent = T)
  return(ghm)
}



#' Create heatmaps with gene signatures.
#'
#' Helper function to create heatmaps with gene signatures. Find the indexes to
#' divide the data frame cells according to their subtype and gene signatures.
#' Helper function of plotSignaturesHeatmaps.
#'
#' @param ann_col data frame; column annotations of cells and populations.
#' @param ann_row data frame; row annotations of genes and populations.
#' @param order.columns character vector; order of columns.
#' @param order.rows character vector; order of rows.
#'
#' @return list; colsep = ref_seps, rowsep = ref_seps_c.
getIndexes <- function(ann_col, ann_row, order.columns, order.rows){
  ann_col$ID <- as.vector(1:nrow(ann_col))
  # Split into cell populations.
  ref_groups <- split(ann_col, as.factor(ann_col$population))
  ref_groups <- lapply(ref_groups, function(x){x$ID})
  ref_groups <- ref_groups[order.columns]
  # Figure out where to draw lines between subtypes in the heatmap.
  ref_seps <- c()
  i_cur_idx <- 0
  order_idx <- c()
  for(ref_grp in ref_groups){
    i_cur_idx <- i_cur_idx + length(ref_grp)
    ref_seps <- c(ref_seps, i_cur_idx)
    order_idx <- c(order_idx, ref_grp)
  }
  # Figure out where to draw lines between gene signatures in the heatmap.
  ann_row$ID <- as.vector(1:nrow(ann_row))
  ref_groups <- split(ann_row, as.factor(ann_row$signature))
  ref_groups <- lapply(ref_groups, function(x){x$ID})
  ref_groups <- ref_groups[order.rows]
  ref_seps_c <- c()
  i_cur_cdx <- 0
  order_cdx <- c()
  for(ref_grp in ref_groups){
    i_cur_cdx <- i_cur_cdx + length(ref_grp)
    ref_seps_c <- c(ref_seps_c, i_cur_cdx)
    order_cdx <- c(order_cdx, ref_grp)
  }
  return(list(colsep = ref_seps, rowsep = ref_seps_c))
}
