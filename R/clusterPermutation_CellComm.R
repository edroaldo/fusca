#' Performs cluster/cell type permutations to estimate an empirical P-value.
#'
#' Sort cluster labels and calculate interaction with mean expression for each
#' permutation. Save interaction to file every 1000 permutations.
#'
#' @param cellrouter CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genelist character vector; genes to use in the analysis.
#' @param nPerm numeric; the number of permutations.
#' @param nCores numeric; the number of cores.
#' @param interactions the interactions as calculated by the
#' population.pairing function or the calculateObservedMean function.
#' @param cluster.label character; column in the metadata table
#' (CellRouter@@sampTab) corresponding to the population.
#'
#' @return list; the interaction data frame for each permutation, including the
#' columns of mean expression of ligands, receptors and pairs.
#'
#' @export
clusterPermutation <- function(cellrouter, assay.type='RNA', genelist,
                               nPerm, interactions, cluster.label, nCores=1){
  cl <- parallel::makeCluster(nCores) 
  doParallel::registerDoParallel(cl) 
  pcellrouter <- cellrouter
  mean.pairs <- list()
  mean.pairs <- foreach (j = 1:nPerm) %dopar% { #for(j in 1:nPerm){
    cat(j, '______________________________\n')
    pclusters <- as.vector(slot(pcellrouter, 'assays')[[assay.type]]@sampTab[[cluster.label]])
    pclusters <- pclusters[sample(length(pclusters))]
    slot(pcellrouter, 'assays')[[assay.type]]@sampTab[[cluster.label]] <- pclusters
    mean.expr <- computeValue(pcellrouter, assay.type, genelist, cluster.label,
                              fun='mean')
    # interactions <- population.pairing(mean.expr = mean.expr, ligands=ligands,
    #                                    receptors=receptors, threshold = 0.25, pairs.m)
    interactions2 <- interactions
    interactions2$mean.ligand <- 0
    interactions2$mean.receptor <- 0
    for(i in 1:nrow(interactions2)){
      c1 <- as.vector(interactions2[i, 'celltype1'])
      c2 <- as.vector(interactions2[i, 'celltype2'])
      l <- as.vector(interactions2[i, 'ligand'])
      r <- as.vector(interactions2[i, 'receptor'])
      interactions2[i, 'mean.ligand'] <- mean.expr[[c1]][l, 'p']
      interactions2[i, 'mean.receptor'] <- mean.expr[[c2]][r, 'p']
    }
    interactions2$mean <- as.numeric(interactions2$mean.ligand +
                                       as.numeric(interactions2$mean.receptor)) / 2
    #rownames(interactions2) <- as.vector(interactions2$pair)
    mean.pairs[[paste('perm',j,sep='')]] <- interactions2
    #if((j %% 1000) == 0){
    #  save(interactions2, file=paste('permutations/permutation_number',j,sep='_'))
    #}
  }
  name <- foreach (j = 1:nPerm, .combine='c')  %do% {paste('perm',j,sep='')} 
  names(mean.pairs) = name; 
  parallel::stopCluster(cl) 
  return(mean.pairs)
}


#' Performs cluster/cell type permutations to estimate an empirical P-value.
#'
#' Sort cluster labels and calculate interaction with mean expression for each
#' permutation. Save interaction to file every 1000 permutations.
#'
#' @param cellrouter CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genelist character vector; genes to use in the analysis.
#' @param nPerm numeric; the number of permutations.
#' @param nCores numeric; the number of cores.
#' @param interactions the interactions as calculated by the
#' population.pairing function or the calculateObservedMean function.
#' @param cluster.label character; column in the metadata table
#' corresponding to the population.
#' @param subcluster.column character; the name of the column where the
#' subclustering information will be stored.
#' @param clusters character; selected clusters.
#'
#' @return list; the interaction data frame for each permutation, including the
#' columns of mean expression of ligands, receptors and pairs.
#'
#' @export
clusterPermutationSubcluster <- function(cellrouter, assay.type='RNA', genelist,
                                         nPerm, interactions, cluster.label,
                                         subcluster.column='Subpopulation',
                                         clusters, nCores=1){
  cl <- parallel::makeCluster(nCores) 
  doParallel::registerDoParallel(cl) 
  pcellrouter <- cellrouter
  sampTab <- slot(pcellrouter, 'assays')[[assay.type]]@sampTab[
    slot(pcellrouter, 'assays')[[assay.type]]@sampTab[[cluster.label]] %in% clusters, ]
  expDat <- slot(pcellrouter, 'assays')[[assay.type]]@ndata[genelist, rownames(sampTab)]
  slot(pcellrouter, 'assays')[[assay.type]]@sampTab <- sampTab
  slot(pcellrouter, 'assays')[[assay.type]]@ndata <- expDat
  mean.pairs <- list()
  mean.pairs <- foreach (j = 1:nPerm) %dopar% { #for(j in 1:nPerm){
    cat(j, '______________________________\n')
    pclusters <- as.vector(slot(pcellrouter, 'assays')[[assay.type]]@sampTab[[subcluster.column]])
    pclusters <- pclusters[sample(length(pclusters))]
    slot(pcellrouter, 'assays')[[assay.type]]@sampTab[[subcluster.column]] <- pclusters
    mean.expr <- computeValueSubclusters(object = pcellrouter,
                                         assay.type = assay.type,
                                         genelist = genelist,
                                         column = cluster.label,
                                         subcluster.column = subcluster.column,
                                         clusters = clusters, fun='mean')
    # interactions <- population.pairing(mean.expr = mean.expr, ligands=ligands,
    #                                    receptors=receptors, threshold = 0.25, pairs.m)
    interactions2 <- interactions
    interactions2$mean.ligand <- 0
    interactions2$mean.receptor <- 0
    for(i in 1:nrow(interactions2)){
      c1 <- as.vector(interactions2[i, 'celltype1'])
      c2 <- as.vector(interactions2[i, 'celltype2'])
      l <- as.vector(interactions2[i, 'ligand'])
      r <- as.vector(interactions2[i, 'receptor'])
      interactions2[i, 'mean.ligand'] <- mean.expr[[c1]][l, 'p']
      interactions2[i, 'mean.receptor'] <- mean.expr[[c2]][r, 'p']
    }
    interactions2$mean <- as.numeric(interactions2$mean.ligand +
                                       as.numeric(interactions2$mean.receptor)) / 2
    #rownames(interactions2) <- as.vector(interactions2$pair)
    mean.pairs[[paste('perm',j,sep='')]] <- interactions2
    #if((j %% 1000) == 0){
    #  save(interactions2, file=paste('permutations/permutation_number',j,sep='_'))
    #}
  }
  name <- foreach (j = 1:nPerm, .combine='c')  %do% {paste('perm',j,sep='')} 
  names(mean.pairs) = name; 
  parallel::stopCluster(cl) 
  return(mean.pairs)
}
