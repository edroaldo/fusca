#' Predict a gene regulatory network.
#'
#' Predict a gene regulatory network using BioConductor annotations.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param species character; species: Hs for Homo Sapiens or Mm for Mus Musculus.
#' @param genes.use character vector; genes to include in the gene regulatory
#' network. The default is to use all genes.
#' @param zscore numeric; zscore threshold to identify putative regulatory
#' interactions.
#' @param blocksize numeric; size of the blocks in which genes will be scaled.
#' @param filename character; filename where the GRN data will be saved.
#'
#' @return list; the GNR with the gene regulatory network graph, the GRN_table
#' with the gene regulatory network table, and the the tfs with the transcription
#' factors.
#'
#' @export
#' @docType methods
#' @rdname buildGRN-methods
setGeneric("buildGRN", function(object, assay.type='RNA',
                                species, genes.use = NULL, zscore = 5,
                                blocksize = 500, filename='GRN.R')
  standardGeneric("buildGRN"))
#' @rdname buildGRN-methods
#' @aliases buildGRN
setMethod("buildGRN",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                species = c('Hs', 'Mm'), genes.use,
                                zscore = 5, blocksize, filename='GRN.R'){
            species <- match.arg(species)
            if(is.null(genes.use)){
              genes.use <- rownames(slot(object, 'assays')[[assay.type]]@ndata)
            }
            # Trancription factors.
            tfs <- find_tfs(species = species)
            # Gene regulated network.
            grn <- globalGRN(slot(object, 'assays')[[assay.type]]@ndata[genes.use, , drop=FALSE],
                             tfs, zscore, blocksize)
            colnames(grn)[1:2]<-c("TG", "TF");
            ggrn <- ig_tabToIgraph(grn, simplify = TRUE)
            x <- list(GRN = ggrn, GRN_table = grn, tfs = tfs)
            save(x, file = filename)
            return(x)
          }
)



#' Predict a gene regulatory network.
#'
#' Predict a gene regulatory network allowing a customized list of transcription
#' factors or genes (tfs).
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param species character; species: Hs for Homo Sapiens or Mm for Mus Musculus.
#' @param tfs character vector; transcription factors to use for GRN reconstruction.
#' @param genes.use character vector; genes to include in the gene regulatory
#' network. The default is to use all genes.
#' @param zscore numeric; zscore threshold to identify putative regulatory
#' interactions.
#' @param blocksize numeric; size of the blocks in which genes will be scaled.
#' @param filename character; filename where the GRN data will be saved.
#'
#' @return list; the GNR with the gene regulatory network graph, the GRN_table
#' with the gene regulatory network table, and the the tfs with the transcription
#' factors.
#'
#' @export
#' @docType methods
#' @rdname buildGRN2-methods
setGeneric("buildGRN2", function(object, assay.type='RNA',
                                 species = c('Hs', 'Mm'), tfs = NULL,
                                 genes.use = NULL, zscore = 5, blocksize = 500,
                                 filename = 'GRN.R')
  standardGeneric("buildGRN2"))
#' @rdname buildGRN2-methods
#' @aliases buildGRN2
setMethod("buildGRN2",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                species, tfs, genes.use, zscore = 5,
                                blocksize, filename = 'GRN.R'){
            if(is.null(genes.use)){
              genes.use <- rownames(slot(object, 'assays')[[assay.type]]@ndata)
            }
            if(is.null(tfs)){
              tfs <- find_tfs(species = species)
            }
            grn <- globalGRN(slot(object, 'assays')[[assay.type]]@ndata[genes.use,],
                             tfs, zscore, blocksize)
            colnames(grn)[1:2] <- c("TG", "TF");
            ggrn <- ig_tabToIgraph(grn, simplify = TRUE)
            x <- list(GRN = ggrn, GRN_table = grn, tfs = tfs)
            save(x, file = filename)
            return(x)
          }
)



#' Find transcription factors.
#'
#' Identify transcriptional regulators based on gene ontology annotations.
#' Helper function of buildGRN.
#'
#' @param species character; species: Hs for human or Mm for mouse.
find_tfs <- function(species = 'Hs'){
  # Hs is species abbreviation.
  cat("Loading gene annotations...\n")
  if(species=='Hs'){
    egSymbols <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL);
    goegs <- as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS);
  }
  else{
    egSymbols <- as.list(org.Mm.eg.db::org.Mm.egSYMBOL);
    goegs <- as.list(org.Mm.eg.db::org.Mm.egGO2ALLEGS);
  }
  goterms <- as.list(GO.db::GOTERM);
  goids <- names(goegs);
  onts <- lapply(goids, AnnotationDbi::Ontology);
  bps <- onts[onts == 'BP'];
  goids <- names(unlist(bps));
  cat("Matching gene symbols and annotations...")
  gobpList<-list();
  for(goid in goids){
    egs <- goegs[[goid]];
    goterm <- AnnotationDbi::Term(goterms[[goid]]);
    genes <- sort(unique(as.vector(unlist(egSymbols[egs]))));
    gobpList[[goterm]] <- genes;
  }
  regNames <- names(gobpList)[grep("regulation of transcription", names(gobpList))];
  trs <- unique(unlist(gobpList[regNames]));
  mfs <- onts[onts == 'MF'];
  goidsMF <- names(unlist(mfs));
  gomfList <- list();
  for(goid in goidsMF){
    egs <- goegs[[goid]];
    goterm <- AnnotationDbi::Term(goterms[[goid]]);
    genes <- sort(unique(as.vector(unlist(egSymbols[egs]))));
    gomfList[[goterm]] <- genes;
  }
  dbs <- gomfList[['DNA binding']];
  sort(intersect(trs, dbs));
}

#' Create global gene regulatory network.
#'
#' Helper function of buildGRN.
#'
#' @param expr dgCMatrix; gene expression in each cell, found in
#' CellRouter@@ndata.
#' @param tfs character vector; transcription factors.
#' @param threshold numeric; zscore threshold to identify regulatory
#' interactions.
#' @param blocksize numeric; size of the blocks in which genes will be scaled.
globalGRN <- function(expr, tfs, threshold, blocksize){
  # expr <- as.data.frame(as.matrix(expr))
  tfs <- intersect(tfs, rownames(expr))
  # Modified to calculate correlation os sparse matrix.
  #corrAll <- abs(sparse.cor(Matrix::t(expr)))
  corrAll <- abs(sparse.cor2(Matrix::t(expr)))
  # corrAll <- abs(cor(t(expr)))
  zscores <- mat_zscores(corrAll, blocksize)
  zscores <- zscores[, tfs]
  grnTable <- extractRegulators(zscores, corrAll, rownames(expr), threshold)
  grnTable
}



#' Create an igraph object for the gene regulatory network.
#'
#' Helper function of buildGRN.
#'
#' @param grnTab matrix; table containing gene regulatory relationships. Table
#' of TF, TF, maybe zscores, maybe correlations.
#' @param simplify boolean; simplify the graph. For example, remove loops.
#' @param directed boolean; create a directed or undirected graph.
#' @param weights boolean; add gene correlations as weights of the graph.
#' @return iGraph object.
ig_tabToIgraph <- function(grnTab, simplify = FALSE, directed = FALSE,
                           weights = TRUE){
  # Simplify failed when iranges is loaded.
  tmpAns <- as.matrix(grnTab[, c("TF", "TG")]);
  regs <- as.vector(unique(grnTab[, "TF"]));
  targs <- setdiff(as.vector(grnTab[,"TG"]), regs);
  myRegs <- rep("Regulator", length = length(regs));
  myTargs <- rep("Target", length = length(targs));
  types <- c(myRegs, myTargs);
  verticies <- data.frame(name = c(regs, targs), label = c(regs, targs),
                          type = types);
  # Greate graph.
  iG <- igraph::graph.data.frame(tmpAns, directed = directed, v = verticies);
  if(weights){
    igraph::E(iG)$weight <- as.numeric(as.character(grnTab$corr));
  }
  if(simplify){
    iG <- igraph::simplify(iG);
  }
  igraph::V(iG)$nEnts <- 1;
  return(iG)
}



#' GRN reconstruction.
#'
#' GRN reconstruction code: from CellNet (Cahan et al., Cell 2014).
#'
#' @param corrMat matrix; correlation matrix.
#' @param blocksize numeric; size of the blocks in which genes will be scaled.
#'
#' @return matrix;
mat_zscores <- function(corrMat, blocksize){
  print('mat_zscores')
  z_row <- scale.sparse(Matrix::t(corrMat), blocksize = blocksize) ** 2;
  z_col <- scale.sparse(corrMat, blocksize=blocksize) ** 2;
  cat("\n", dim(z_col), "\n");
  ans <- sqrt(z_row + z_col);
  ans;
}




#' Extract transcriptional regulators.
#'
#' Helper function globalGNR.
#'
#' @param zscores matrix; zscores.
#' @param corrMatrix matrix; correlation matrix.
#' @param genes character vector; genes.
#' @param threshold numeric; threshold.
extractRegulators <- function(zscores, corrMatrix, genes, threshold){
  targets <- vector()
  regulators <- vector()
  zscoresX <- vector()
  correlations <- vector()
  targets <- rep('', 1e6)
  regulators <- rep('', 1e6)
  zscoresX <- rep('', 1e6)
  correlations <- rep('', 1e6)
  i <- 1
  j <- 1
  k = 1
  zscores <- as.matrix(zscores)
  corrMatrix <- as.matrix(corrMatrix)
  for(target in genes){
    x <- zscores[target,]
    regs <- names(which(x > threshold))
    if(length(regs) > 0){
      zzs <- x[regs]
      ncount <- length(zzs)
      corrs <- as.numeric(corrMatrix[target, regs])
      j <- i + ncount - 1
      targets[i:j] <- rep(target, ncount)
      regulators[i:j] <- regs
      zscoresX[i:j] <- zzs
      correlations[i:j] <- corrs
      i <- j + 1
    }
    k <- k + 1
    if((k %% 1000) == 0){
      cat(k, ' of ', length(genes), ' processed\n')
    }
  }
  targets <- targets[1:j]
  regulators <- regulators[1:j]
  zscoresX <- as.numeric(zscoresX[1:j])
  correlations <- as.numeric(correlations[1:j])
  data.frame(target=targets, reg=regulators, zscore=zscoresX, corr=correlations)
}
