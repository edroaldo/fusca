#' Pathway enrichment analysis
#'
#' Pathway enrichment analysis on genes up or down-regulated along trajectories.
#'
#' @param object CellRouter object.
#' @param names character vector; selected trajectories.
#' @param cc character vector; list of cell cycle genes or other gene lists to
#' be removed from the genes up or down-regulated in each trajectory.
#' @param species character; species: Hs for Homo Sapiens or Mm for Mus Musculus.
#' @param annotation character; organism-specific annotations: 'mouse' or
#' 'org.Mm.eg.db', or 'human' or 'org.Hs.eg.db'.
#' @param ids table containing mappings between gene identifiers.
#'
#' @return CellRouter object with the pathwayenrichment slot updated.
#'
#' @export
#' @docType methods
#' @rdname pathwayenrichment-methods
setGeneric("pathwayenrichment", function(object, names, cc=NULL,
                                         species=c('Hs', 'Mm'),
                                         annotation, ids)
  standardGeneric("pathwayenrichment"))
#' @rdname pathwayenrichment-methods
#' @aliases pathwayenrichment
setMethod("pathwayenrichment",
          signature="CellRouter",
          definition=function(object, names, cc=NULL, species=c('Hs', 'Mm'),
                              annotation, ids){
            species <- match.arg(species)
            # Up-regulated.
            upregulated <- list()
            geneNames <- lapply(object@top.correlations$up, names)
            if(!is.null(cc)){
              cat('removing genes...\n')
              # Remove a selected gene set.
              geneNames <- lapply(geneNames, function(x){setdiff(x, cc)})
            }
            geneList <- lapply(geneNames, function(x){
              convertIDs(ids, x, from='external_gene_name', to="entrezgene")})
            geneList <- lapply(geneList, names)
            geneList <- lapply(geneList, function(x){x[!is.na(x)]})
            print('pathway enrichment for up-regulated genes')
            #
            # duvida: nao achei essa funcao
            #
            ck1 <- compareCluster(geneCluster = geneList[names],
                                  fun = "enrichPathway",
                                  organism = species,
                                  pvalueCutoff = 0.05,
                                  readable=T)
            ck3 <- compareCluster(geneCluster = geneList[names],
                                  fun = "enrichGO",
                                  ont='BP',
                                  OrgDb=annotation,
                                  pvalueCutoff = 0.05,
                                  readable=T)
            upregulated[['REACTOME']] <- ck1
            upregulated[['GOBP']] <- ck3
            # Down-regulated
            downregulated <- list()
            geneNames <- lapply(object@top.correlations$down, names)
            geneList <- lapply(geneNames, function(x){
              convertIDs(ids, x, from='external_gene_name', to="entrezgene")})
            geneList <- lapply(geneList, names)
            print('pathway enrichment for down-regulated genes')
            ck1 <- compareCluster(geneCluster = geneList[names],
                                  fun = "enrichPathway",
                                  organism = species,
                                  pvalueCutoff = 0.05,
                                  readable=T)
            ck3 <- compareCluster(geneCluster = geneList[names],
                                  fun = "enrichGO",
                                  ont='BP',
                                  OrgDb=annotation,
                                  pvalueCutoff = 0.05,
                                  readable=T)
            downregulated[['GOBP']] <- ck3
            # Writes pathway enrichment to CellRouter.
            object@pathwayenrichment <- list("UP"=upregulated,
                                             "DOWN"=downregulated)
            return(object)
          })

#' Helper function of convertIDs.
#'
#' Convert gene identifiers.
#'
#' @param ids list; table containing mappings between gene identifiers.
#' @param ens character; gene name.
#' @param from character; gene identifier.
#' @param to character; gene identifier.
#'
#' @return list; gene list.
#'
#' @export
convertID <- function(ids, ens, from, to){
  symbol <-  ids[which(ids[,from]==ens), to]
  if(length(symbol) == 0){
    ens
  }else{
    symbol[1]
  }
}



#' Helper function of pathwayenrichment.
#'
#' Convert gene lists from one identifier to another.
#'
#' @param ids list; table containing mappings between gene identifiers.
#' @param geneList character vector; gene list to be converted.
#' @param from character; gene identifier.
#' @param to character; gene identifier.
#'
#' @return list; gene list.
#'
#' @export
convertIDs <- function(ids, geneList, from, to){
  genes <- vector()
  for(gene in geneList){
    g <- convertID(ids, gene, from, to)
    genes <- append(genes, g)
  }
  names(geneList) <- genes
  genes <- genes[!is.na(genes)]
  geneList
}
