#' Plot top-ranked genes in transcriptional clusters.
#'
#' Plot top-ranked genes calculated by the rankGenesTranscriptionalClusters
#' function and saves to file.
#'
#' @param assay.type character; the type of data to use.
#' @param rankedDynamics list; ranked dynamics as calculated by the
#' rankGenesTranscriptionalClusters function.
#' @param traj.name character; trajectory of interest.
#' @param num.genes numeric; number of genes to plot.
#'
#' @return ggplot2 graph.
#'
#' @export
#' @docType methods
#' @rdname plotGenesInClusters-methods
setGeneric("plotGenesInClusters", function(assay.type='RNA', rankedDynamics,
                                           traj.name, num.genes)
  standardGeneric("plotGenesInClusters"))
#' @rdname plotGenesInClusters-methods
#' @aliases plotGenesInClusters
setMethod("plotGenesInClusters",
          signature="list",
          definition=function(assay.type, rankedDynamics, traj.name,
                              num.genes){

            trajectory <- rankedDynamics[[traj.name]]
            traj.plots <- list()
            for(cl.name in unique(as.vector(rankedDynamics[[traj.name]]$
                                            cluster))){
              x <- trajectory[which(trajectory$cluster == cl.name),]
              x <- as.vector(x[1:num.genes,'genes'])
              #
              # duvida: esse test pode colocar outra coisa? Daria pra colocar um filename e alguma coisa, se deixar assim acho que buga se nao for unix.
              #
              g <- plottrajectory(cellrouter, assay.type,
                                  cellrouter@pathsinfo$path[[traj.name]],
                                  x)
              g <- g + ggplot2::ggtitle(paste('cluster', cl.name,sep='_'))
              g <- g +
                ggplot2::guides(col = ggplot2::guide_legend(direction = "vertical",
                                                            keywidth = 0.5,
                                                            keyheight = 0.85))
              g <- g +
                ggplot2::theme(legend.text = ggplot2::element_text(colour="black",
                                                                   size=7))
              traj.plots[[paste('cluster', cl.name,sep='_')]] <- g
            }
            return(traj.plots)
          }
)
