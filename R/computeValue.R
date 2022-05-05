#' Compute expression statistics of each gene in each population.
#'
#' Compute maximum, median or mean expression of each gene in each population.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genelist character vector; genes to use in the analysis.
#' @param column character; column in the metadata table to group cells for
#' differential expression. For example, if 'population' is specified,
#' population-specific gene signatures will be identified.
#' @param fun character; statistical function to summary the gene expression.
#'
#' @return list; the statistics of the expressed genes in the population (p
#' column), not in the population (np), and the percentage (percent).
#'
#' @export
#' @docType methods
#' @rdname computeValue-methods
setGeneric("computeValue", function(object, assay.type='RNA',
                                    genelist, column='population', fun='max', verbose = TRUE)
  standardGeneric("computeValue"))
#' @rdname computeValue-methods
#' @aliases computeValue
setMethod("computeValue",
          signature = "CellRouter",
          definition = function(object, assay.type, genelist, column, fun, verbose){
            if(verbose == TRUE){
              print('discovering subpopulation-specific gene signatures')
            }
            expDat <- slot(object, 'assays')[[assay.type]]@ndata[genelist,]
            membs <- as.vector(slot(object, 'assays')[[assay.type]]@sampTab[[column]])
            membs_df <- as.vector(slot(object, 'assays')[[assay.type]]@sampTab[ , c('sample_id', column), drop=FALSE])
            diffs <- list()
            for(i in unique(membs)){
              if(verbose == TRUE){
                cat('cluster ', i, '\n')
              }
              if(sum(membs == i) == 0) next
              m_indexes <- membs_df[which(membs_df[[column]] != i), 'sample_id']
              n_indexes <- membs_df[which(membs_df[[column]] == i), 'sample_id']
              if(fun == 'max'){
                m <- if(sum(membs != i) > 1) apply(expDat[, m_indexes], 1, max) else expDat[, m_indexes]
                n <- if(sum(membs == i) > 1) apply(expDat[, n_indexes], 1, max) else expDat[, n_indexes]
              }else if(fun == 'median'){
                m <- if(sum(membs != i) > 1) apply(expDat[, m_indexes], 1, median) else expDat[, m_indexes]
                n <- if(sum(membs == i) > 1) apply(expDat[, n_indexes], 1, median) else expDat[, n_indexes]
              }else{
                m <- if(sum(membs != i) > 1) apply(expDat[, m_indexes], 1, mean) else expDat[, m_indexes]
                n <- if(sum(membs == i) > 1) apply(expDat[, n_indexes], 1, mean) else expDat[, n_indexes]
              }
              d <- data.frame(np=m, p=n) #log scale
              diffs[[i]] <- d
            }
            return (diffs)
          }
)

detectGenes <- function(expr, min_expr = 0.1){
  detected <- list()
  genes_detected <- do.call(rbind, apply(expr, 1, function(x){return (data.frame(num_cells_detected=sum(x > min_expr)))}))
  detected[['genes_detected']] <- genes_detected
  detected
}


#' Compute expression statistics of each gene in each population.
#'
#' Compute maximum, median or mean expression of each gene in each population.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genelist character vector; genes to use in the analysis.
#' @param column character; column in the metadata table to group cells for
#' differential expression. For example, if 'population' is specified,
#' population-specific gene signatures will be identified.
#' @param subcluster.column character; the name of the column where the
#' subclustering information will be stored.
#' @param clusters character; selected clusters.
#' @param fun character; statistical function to summary the gene expression.
#'
#' @return list; the statistics of the expressed genes in the population (p
#' column), not in the population (np), and the percentage (percent).
#'
#' @export
#' @docType methods
#' @rdname computeValueSubclusters-methods
setGeneric("computeValueSubclusters", function(object, assay.type='RNA',
                                               genelist,
                                               column='population',
                                               subcluster.column='Subpopulation',
                                               clusters,
                                               fun='max',
                                               verbose = TRUE)
  standardGeneric("computeValueSubclusters"))
#' @rdname computeValueSubclusters-methods
#' @aliases computeValueSubclusters
setMethod("computeValueSubclusters",
          signature = "CellRouter",
          definition = function(object, assay.type, genelist, column,
                                subcluster.column, clusters, fun, verbose){
            if(verbose == TRUE){
              print('discovering subpopulation-specific gene signatures')
            }
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab[
              slot(object, 'assays')[[assay.type]]@sampTab[[column]] %in% clusters,]
            expDat <- slot(object, 'assays')[[assay.type]]@ndata[genelist, rownames(sampTab)]
            #membs <- as.vector(object@sampTab$population)
            membs <- as.vector(sampTab[[subcluster.column]])
            diffs <- list()
            for(i in unique(membs)){
              if(verbose == TRUE){
                cat('cluster ', i, '\n')
              }
              if(sum(membs == i) == 0) next
              if(fun == 'max'){
                m <- if(sum(membs != i) > 1) apply(expDat[, membs != i], 1, max) else expDat[, membs != i]
                n <- if(sum(membs == i) > 1) apply(expDat[, membs == i], 1, max) else expDat[, membs == i]
              }else if(fun == 'median'){
                m <- if(sum(membs != i) > 1) apply(expDat[, membs != i], 1, median) else expDat[, membs != i]
                n <- if(sum(membs == i) > 1) apply(expDat[, membs == i], 1, median) else expDat[, membs == i]
              }else{
                m <- if(sum(membs != i) > 1) apply(expDat[, membs != i], 1, mean) else expDat[, membs != i]
                n <- if(sum(membs == i) > 1) apply(expDat[, membs == i], 1, mean) else expDat[, membs == i]
              }
              d <- data.frame(np=m, p=n) #log scale
              diffs[[i]] <- d
            }
            return (diffs)
          }
)
