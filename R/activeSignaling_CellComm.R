#' Create gene regulator subnetworks based on gene signature value.
#'
#' Create subnetworks based on the gene signature value for each node in the
#' summary and for each target in the node.
#'
#' @param data dgCMatrix; gene expression in each cell, found in
#' CellRouter@@ndata.
#' @param sampTab data frame; the metadata information of the cells, found in
#' CellRouter@@sampTab.
#' @param markers data frame; the gene signatures as calculated by the
#' findSignatures function of CellRouter.
#' @param grn data frame; the gene regulatory network table (GRN_table) as
#' calculated by the buildGRN method of CellRouter.
#' @param summary data frame; the summary of the files obtained in the
#' findpaths.simple function.
#' @param fc numeric; keep only vertex which signature value is smaller than the
#' fold changes.
#'
#' @return list; the pathways and the gene regulatory networks (gnrs). The
#' pathways contain the lists of genes that are the markers and the
#' regulators of the target populations for each interaction. The gnrs
#' contain the regultory sub network for each node in the summary and for each
#' target in the node.
#'
#' @export
activeSignaling <- function(data, sampTab, markers,  grn, summary, fc = 5){
  bla <- list()
  nodes <- summary$nodes
  grns <- list()
  for(sub in names(nodes)){
    #x <- data[,rownames(sampTab[sampTab$population == sub,])]
    tmp <- nodes[[sub]]
    targets <- rownames(tmp[which(tmp$category == 'TARGET'),])
    seed <- rownames(tmp[which(tmp$category == 'SOURCE'),])
    celltype <- strsplit(sub, split='_')[[1]][2]
    print(celltype)
    #cors <- vector()
    #names <- vector()
    cat('SUB: ', sub,' SEED: ', seed, '\n')
    cat('SUB: ', sub,' TARGETS: ', targets, '\n')
    overlaps <- list()
    for(g in targets){
      xx <- grn[which(grn$reg == g),] #targets of regulator g
      xxx <- as.vector(xx$target)

      queryGenes <- as.vector(xx$target)
      nGenesQuery <- length(queryGenes);
      universe <- rownames(data)
      #annGenes <- intersect(universe, rownames(cellchat@signatures[[sub]]))
      #genesBoth <- intersect(queryGenes, rownames(cellchat@signatures[[sub]]))
      signature <- markers[which(markers$population == celltype),]
      annGenes <- intersect(universe, rownames(signature))
      genesBoth <- intersect(queryGenes, rownames(signature))

      nBoth <- length(genesBoth)
      if(nBoth > 0){
        genesQueryOnly <- setdiff(queryGenes, genesBoth);
        genesAnnOnly   <- setdiff(annGenes,   genesBoth);
        genesNeither   <- setdiff(universe,   c(queryGenes, annGenes));
        nQueryOnly     <- length(genesQueryOnly);
        nAnnOnly       <- length(genesAnnOnly  );
        nNeither       <- length(genesNeither  );

        mt <- matrix(c(nBoth,nQueryOnly,nAnnOnly, nNeither ), nrow=2, byrow=T);
        #print(mt)
        tResult <- fisher.test(mt);
        pvals <- tResult$p.value;
        if(pvals < 0.05){
          cat(g, ': ', nBoth, '\t', dim(xx), '\t ', pvals,'\n')
          overlaps[[g]] <- genesBoth
        }
        if(length(overlaps) > 0){
          bla[[sub]] <- overlaps
        }
        ##create regulator subnetwork
        xxtmp <- xx
        xxtmp <- xxtmp[xxtmp$target %in% genesBoth,]
        #print(xxtmp)
        colnames(xxtmp) <- c('target', 'reg', 'zscore', 'weight')
        gxx <- igraph::graph.data.frame(xxtmp[ ,c('reg', 'target', 'weight')],
                                        directed = TRUE)
        igraph::V(gxx)$type <- 'target'
        #V(gxx)[g]$type <- 'regulator'
        igraph::V(gxx)[c(g, intersect(genesBoth, as.vector(xxtmp$reg)))]$
          type <- 'regulator'
        #keep <- intersect(rownames(cellchat@signatures[[sub]]), V(gxx)$name)
        #V(gxx)[keep]$size <- log2(cellchat@signatures[[sub]][keep,'fc'])
        keep <- intersect(rownames(signature), igraph::V(gxx)$name)
        igraph::V(gxx)[keep]$size <- signature[keep,'fc']

        igraph::V(gxx)$label <- igraph::V(gxx)$name
        igraph::V(gxx)[igraph::V(gxx)$size < fc &
                         igraph::V(gxx)$type != 'regulator']$label <- ""

        grns[[sub]][[g]] <- gxx

      }
    }
    cat('______________________________\n')
  }
  return(list(pathways=bla, grns=grns))
}
