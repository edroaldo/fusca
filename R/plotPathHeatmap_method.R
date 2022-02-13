#' Plot path heatmap.
#'
#' Plot heatmap for specified genes in the trajectories and save plots to
#' directory.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param paths character vector; selected trajectories.
#' @param genelist character vector; genes to show.
#' @param cluster.column character; the name of the column where the clustering
#' information is stored.
#' @param color.column character; the name of the column where the color
#' information is stored.
#' @param threshold numeric; threshold to rescale gene expression.
#'
#' @return list; ggplot2 plots.
#'
#' @export
#' @docType methods
#' @rdname plotPathHeatmap-methods
setGeneric("plotPathHeatmap", function(object, assay.type='RNA',
                                       paths, genelist,
                                       cluster.column='population',
                                       color.column='colors',
                                       threshold=2)
  standardGeneric("plotPathHeatmap"))
#' @rdname plotPathHeatmap-methods
#' @aliases plotPathHeatmap
setMethod("plotPathHeatmap",
          signature="CellRouter",
          definition=function(object, assay.type,
                              paths, genelist,
                              cluster.column,
                              color.column,
                              threshold){
            corsPaths <- object@correlation
            pathsInfo <- object@pathsinfo
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            plots <- list()
            for(path in paths){
              genelist2 <- intersect(genelist,
                                     rownames(pathsInfo$distr[[path]]))
              tmpexpr <- pathsInfo$distr[[path]][genelist2,]
              tmpexpr <- center_with_threshold(tmpexpr, 1.5) #double check this...
              andf <- data.frame(sampTab[pathsInfo$path[[path]], cluster.column, ])
              rownames(andf) <- pathsInfo$path[[path]]
              colnames(andf) <- c('subpopulation')
              target_colors <- unique(sampTab[pathsInfo$path[[path]],
                                              color.column, ])
              names(target_colors) <- unique(andf$subpopulation)
              ann_colors = list(
                subpopulation = target_colors
              )
              from <- sapply(strsplit(path, split='.', fixed=TRUE),
                             function(x){x[1]})
              to <- sapply(strsplit(path, split='.', fixed=TRUE),
                           function(x){x[2]})
              title <- paste('Transition ', from, ' ', to, sep = '')
              labels <- sapply(strsplit(rownames(tmpexpr), split = '__',
                                        fixed = TRUE), function(x){x[1]})
              g <- pheatmap::pheatmap(tmpexpr, cluster_rows = FALSE,
                                 cluster_cols = FALSE, annotation_col = andf,
                                 annotation_colors=ann_colors,
                                 show_colnames = FALSE, border = FALSE,
                                 main = title, labels_row = labels, silent = T)
              plots[[path]] <- g
            }
            return(plots)
          }
)
