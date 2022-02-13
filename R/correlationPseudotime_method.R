#' Compute the correlation of each gene along the trajectories.
#'
#' Computes the correlation between the gene expression and the pseudotime of
#' for each trajectory using a specified method.
#'
#' @param object CellRouter object.
#' @param type character; correlation method: linear regression (slope),
#' spearman or pearson.
#'
#' @return CellRouter object with the correlation slot updated.
#'
#' @export
#' @docType methods
#' @rdname correlationPseudotime-methods
setGeneric("correlationPseudotime", function(object, type = c('slope',
                                                              'spearman',
                                                              'pearson'))
  standardGeneric("correlationPseudotime"))
#' @rdname correlationPseudotime-methods
#' @aliases correlationPseudotime
setMethod("correlationPseudotime", signature = "CellRouter",
          definition = function(object, type = c('slope',
                                                 'spearman',
                                                 'pearson')){
            type <- match.arg(type)
            pathsInfo <- object@pathsinfo
            genelist <- object@genes.trajectory
            correlations <- list()
            #ndata <- slot(object, 'assays')[[assay.type]]@ndata
            print('computing correlation with the pseudotime')
            for(path in names(pathsInfo[['distr']])){
              cat(path, '\n')
              # Cells in the path and their pseudotime.
              x <- pathsInfo[['pseudotime']][[path]]
              genes <- intersect(genelist, rownames(pathsInfo[['distr']][[path]]))
              #genes <- intersect(genelist, pathsInfo[['distr']][[path]])
              #cells <- object@pathsinfo$path[[path]]
              #tmp.ndata <- ndata[,cells]
              for(gene in genes){
                # Select gene expression for the cells in the path.
                y <- as.numeric(pathsInfo[['distr']][[path]][gene,])
                #y <- as.numeric(tmp.ndata[gene, ])
                # Fit a model for the gene expression where the x is the
                # pseudotime and y is the expression.
                if(type == 'slope'){
                  df <- data.frame(x=x, y=y)
                  df <- lm(y ~ x, data=df)
                  correlations[[path]][[gene]] <- df$coefficients[2]
                }else if(type == 'pearson'){
                  cors <- suppressWarnings(cor.test(x, y, method="pearson"))
                  correlations[[path]][[gene]] <- as.numeric(cors$estimate)
                }else{
                  cors <- suppressWarnings(cor.test(x, y, method="spearman"))
                  correlations[[path]][[gene]] <- as.numeric(cors$estimate)
                }
                correlations[[path]] <- unlist(correlations[[path]])
              }
            }
            # Add the correlations to the CellRouter object.
            object@correlation <- correlations
            return(object)
          }
)
