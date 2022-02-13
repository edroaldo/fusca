#' Plot pathway enrichment analysis heatmap.
#'
#' Plot pathway enrichment analysis heatmap and save to file.
#'
#' @param object CellRouter object.
#' @param pathway the pathway database: GO, Reactome.
#' @param numpathways numeric; number of pathways to show for each trajectory.
#' @param logTransform boolean; whether to log-transform the p-values.
#'
#' @return list; the pathway data frame and the plot.
#'
#' @export
#' @docType methods
#' @rdname pathwaycluster-methods
setGeneric("pathwaycluster", function(object, pathway, numpathways=5,
                                      logTransform)
  standardGeneric("pathwaycluster"))
#' @rdname pathwaycluster-methods
#' @aliases pathwaycluster
setMethod("pathwaycluster",
          signature="CellRouter",
          definition=function(object, pathway, numpathways=5, logTransform){
            pathsInfo <- object@pathsinfo
            #
            # duvida: esse pathway é que formato e vem da onde?
            #
            b <- pathway@compareClusterResult
            # Write csv file with the pathway.
            write.csv(b, paste(filename, '.csv', sep=''))
            b2 <- lapply(split(b, as.vector(b$Cluster)),
                         function(x){x[order(as.vector(x$p.adjust),
                                             decreasing=FALSE), ][1:numpathways, ]})
            b <- do.call(rbind, b2)
            b <- b[complete.cases(b), ]
            pathways <- data.frame(matrix(1, nrow=length(unique(b$Description)),
                                          ncol=length(unique(b$Cluster))))
            rownames(pathways) <- as.vector(unique(b$Description))
            colnames(pathways) <- as.vector(unique(b$Cluster))
            for(c in colnames(pathways)){
              for(d in rownames(pathways)){
                x <- b[which(as.vector(b$Cluster) == c &
                               as.vector(b$Description) == d),]
                if(nrow(x) == 1){
                  pathways[as.character(x$Description),
                           as.character(x$Cluster)] <- as.numeric(x$p.adjust)
                }
              }
            }
            if(logTransform){
              pathways <- -log(pathways)
            }
            pathways <- pathways[,order(colnames(pathways))]
            ann <- pathsInfo$path_data[colnames(pathways),]
            source_colors <- unique(ann$source_color)
            names(source_colors) <- unique(ann$source_population)
            target_colors <- unique(ann$target_color)
            names(target_colors) <- unique(ann$target_population)
            ann_colors = list(
              source_population = source_colors,
              target_population = target_colors
            )
            g <- pheatmap::pheatmap(pathways, border=TRUE, border_color = 'gray',
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     cellwidth=8, cellheight = 8, annotation_col = ann[,c(4,6)],
                     annotation_colors=ann_colors, clustering_method = 'ward.D',
                     silent = T)
            return(list(pathways = pathways, plot = g))
          }
)
