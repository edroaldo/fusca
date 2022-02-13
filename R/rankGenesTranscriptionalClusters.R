#' Rank genes in transcriptional clusters.
#'
#' Rank genes in each transcriptional cluster based on their correlation with
#' the mean kinetic profile.
#'
#' @param object CellRouter object.
#' @param num.genes numeric; top num.genes in each transcriptional cluster.
#'
#' @return list; cluster trends for each trajectory.
#'
#' @export
#' @docType methods
#' @rdname rankGenesTranscriptionalClusters-methods
setGeneric("rankGenesTranscriptionalClusters", function(object, num.genes)
  standardGeneric("rankGenesTranscriptionalClusters"))
#' @rdname rankGenesTranscriptionalClusters-methods
#' @aliases rankGenesTranscriptionalClusters
setMethod("rankGenesTranscriptionalClusters",
          signature="CellRouter",
          definition = function(object, num.genes){
            clusters <- object@clusters
            clusterTrends <- list()
            for(trajectory in names(clusters)){
              tmp <- list()
              exprs <- clusters[[trajectory]][['exprs']]
              clustering <- clusters[[trajectory]][['clustering']]
              for(cluster in 1:max(clustering)){
                x <- clustering[which(clustering==cluster)]
                xx <- exprs[names(x),]
                xxx <- apply(xx, 2, mean)
                cors <- apply(exprs, 1, function(x){cor(as.numeric(x),
                                                        as.numeric(xxx),
                                                        method = 'spearman')})
                cors <- cors[order(cors, decreasing=TRUE)]
                cors <- cors[1:num.genes]
                df <- data.frame(genes=names(cors), cors=cors, cluster=cluster,
                                 trend='positive')
                tmp[[paste('cl', cluster,sep='')]] <- df
              }
              clusterTrends[[trajectory]] <- do.call(rbind, tmp)
            }
            return(clusterTrends)
          }
)
