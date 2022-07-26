#' Add metadata information to CellRouter object.
#'
#' Include metadata information in the CellRouter metadata table
#' \code{assay@@sampTab}. If you have clusters identified by other algorithm,
#' you can use this function to perform all the analysis available in
#' CellRouter using your previously identified clusters: just select the
#' corresponding columns in the table \code{assay@@sampTab}.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param metadata data frame or vector; metadata to be added.
#' @param colname character; column name to be added to \code{assay@@sampTab}
#' in case metadata is a vector.
#' @param metadata.column character; column to selected from the metadata to be
#' added and included in \code{assay@@sampTab}.
#'
#' @return CellRouter object with the sampTab slot updated.
#' @export
addInfo <- function(object, assay.type='RNA', sample.name='Sample1',
                    metadata, colname, metadata.column='population'){
  # Update to include data.frames as well.
  # Comented to add other types of metadata, not only the ones from the tutorial.
  # if (assay.type=='ST'){
  #   rownames(metadata) <- metadata[['Barcode']]
  # }
  # sampTab has cells as rows.
  sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
  # Filter spots that are only present in sampTab.
  if (assay.type=='ST'){
    metadata <- metadata[rownames(sampTab), ]
  }
  if(class(metadata) == 'data.frame'){
    # I could use merge here to avoid losing data.
    sampTab[rownames(metadata), colname] <- as.vector(metadata[[metadata.column]])
    # sampTab <- sampTab[order(sampTab[[colname]]), ]
  }else{
    sampTab[names(metadata), colname] <- metadata
  }
  # Remove NAs.
  # sampTab <- sampTab[!is.na(sampTab[[colname]]), ]
  # Convert unique values in colname to colors.
  n_colors <- length(unique(metadata[[metadata.column]]))
  # Select 8 colors to interpolate.
  colors <- cRampClust(1:length(unique(sampTab[[colname]])), 8)
  # Atribute each colname to a color.
  names(colors) <- unique(sampTab[[colname]])
  # Count the number of cells in each cluster.
  replicate_row <- as.vector(unlist(lapply(split(sampTab, sampTab[[colname]]), nrow)))
  # Replicates the colors for all cells in the cluster.
  #colors_row <- rep(colors, times=replicate_row)
  colors_row = c();for (i in sampTab[[colname]]){ colors_row = rbind(colors_row, colors[names(colors) == i])} #!
  color.column <- paste(colname, 'color', sep='_')
  sampTab[, color.column] <- colors_row
  slot(object, 'assays')[[assay.type]]@sampTab <- sampTab
  # Add information to the barcodes dataframe for plotting.
  if (assay.type=='ST'){
    bcs <- slot(object, 'assays')[['ST']]@image$bcs[[sample.name]]
    clusters <- slot(object, 'assays')[['ST']]@
      sampTab[, c('sample_id', colname, color.column)]
    bcs <- merge(bcs, clusters,
                 by.x = "barcode", by.y = "sample_id", all = TRUE)
    bcs <- bcs[!is.na(bcs$tissue), ]
    slot(object, 'assays')[['ST']]@image$bcs[[sample.name]] <- bcs
  }
  # Select only cells in sampTab.
  slot(object, 'assays')[[assay.type]]@ndata <-
    slot(object, 'assays')[[assay.type]]@ndata[, which(
      colnames(slot(object, 'assays')[[assay.type]]@ndata) %in% rownames(sampTab))]
  return(object)
}
