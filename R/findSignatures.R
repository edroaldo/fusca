#' Identify gene signatures.
#'
#' Identify gene signatures using the wilcoxon test or template-matching
#' methods.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param column character; specify the groups to compare.
#' @param test.use character; differential expression test to use. Alternative
#' is based on template-matching. Possible values: wilcox or template.
#' @param pos.only logical; use only upregulated genes.
#' @param min.pct numeric; minimum percentage.
#' @param fc.threshold numerical; fold change threshold.
#' @param fc.tm logical; whether to include fold change values in the
#' template-matching differential expression analysis.
#' @param nCores numeric; the number of cores.
#' @param max.cells numeric; maximum number of cells per cluster (for large clusters)
#'
#' @return data frame; the population, fold change (fc), mean p, mean np,
#' p-value, p-value adjusted, and gene.
#'
#' @export
findSignatures <- function(object, assay.type='RNA',
                           column='population', test.use='wilcox',
                           pos.only=TRUE, min.pct=0.1, fc.threshold=0.25,
                           fc.tm=FALSE, nCores=1, max.cells=Inf){ #!
  if(test.use == 'wilcox'){
    cat('Calculating fold changes...', '\n')
    object <- computeFC(object, assay.type, column=column, min.pct=min.pct,
                        pos.only=pos.only, fc.threshold=fc.threshold, nCores=nCores, max.cells=max.cells) #!
    if (length(object@signatures) == 0){
      print('No population statisfies the selected fc.threshold.')
      markers <- NULL
    } else {
      markers <- findmarkers(object)
    }
  }else if(test.use == 'template'){
    cat('Calculating template-matchings...', '\n')
    signatures <- ranked_findSpecGenes(slot(object, 'assays')[[assay.type]]@ndata,
                                       slot(object, 'assays')[[assay.type]]@sampTab,
                                       qtile=0.99,
                                       remove=TRUE, dLevel=column)
    mylist <- signatures
    for(i in 1:length(mylist) ){ mylist[[i]] <- cbind(mylist[[i]], population=rep(names(mylist[i]), nrow(mylist[[i]])) ) }
    markers <- as.data.frame(do.call(rbind, mylist))
    rownames(markers) <- as.vector(markers$gene)
    # Add fold changes to template method.
    if(fc.tm){
      object <- computeFC(object, assay.type, column, pos.only, fc.threshold, nCores) #!
      for(signature in names(object@signatures)){
        markers[rownames(signatures[[signature]]), 'log2FC'] <- object@signatures[[signature]][rownames(signatures[[signature]]), 'fc']
        markers[rownames(signatures[[signature]]), 'log2FC_pval'] <- object@signatures[[signature]][rownames(signatures[[signature]]), 'pv']
        markers[rownames(signatures[[signature]]), 'log2FC_p.adj'] <- object@signatures[[signature]][rownames(signatures[[signature]]), 'p.adj']
      }
    }
  }
  markers
}



#' Compute fold changes and find gene signatures.
#'
#' Compute fold changes and find gene signatures in each group of cells.
#'
#' @param object CellRouter object
#' @param assay.type character; the type of data to use.
#' @param column character; column in the metadata table to group cells for
#' differential expression. For example, if 'population' is specified,
#' population-specific gene signatures will be identified.
#' @param pos.only logical; only uses genes upregulated.
#' @param min.pct numeric; minimum percentage.
#' @param fc.threshold numerical; fold change threshold.
#' @param nCores numeric; the number of cores. 
#' @param max.cells numeric; maximum number of cells per cluster (for large clusters)
#'
#' @return CellRouter object with the signatures slot updated. This slot
#' contains a list with a dataframe for each group containing information about
#' mean.np, mean.p, fc (fold changes), pv (p-value), p.adj (p-value adjusted.).
#'
#' @export
#' @docType methods
#' @rdname computeFC-methods
setGeneric("computeFC", function(object, assay.type='RNA',
                                 column='population', pos.only=TRUE,
                                 min.pct=0.1, fc.threshold=0.25, nCores=1, max.cells=Inf) #!
  standardGeneric("computeFC"))
#' @rdname computeFC-methods
#' @aliases computeFC
setMethod("computeFC",
          signature="CellRouter",
          definition=function(object, assay.type,
                                column, pos.only, min.pct, fc.threshold, nCores, max.cells){ #!
            print('Identify cluster-specific gene signatures')
            expDat <- slot(object, 'assays')[[assay.type]]@ndata
            membs <- as.vector(slot(object, 'assays')[[assay.type]]@sampTab[[column]])
            membs_df <- as.data.frame(slot(object, 'assays')[[assay.type]]@sampTab[ , c('sample_id', column), drop=FALSE])
            samptab = slot(object, 'assays')[[assay.type]]@sampTab 
            if (max.cells < Inf) {
              set.seed(42); membs_df = data.frame(); subSamp = data.frame();
              for (i in unique(membs)){ #if(sum(membs == i) == 0) next
                submembs_df = data.frame();
                if (length(samptab[[column]][samptab[[column]] == i]) > max.cells) {
                  submembs_df <- sample_n(samptab[samptab[[column]] == i,], max.cells)
                } else {submembs_df <- samptab[samptab[[column]] == i,]}
                subSamp = rbind(subSamp, submembs_df)
              }
              membs_df = as.data.frame(subSamp[ , c('sample_id', column), drop=FALSE])
              membs <- as.vector(subSamp[[column]])
            }
            diffs <- list()
            cl <- parallel::makeCluster(nCores, outfile = "") #!
            doParallel::registerDoParallel(cl) #!
            diffs <- foreach (i = unique(membs), .packages = c("fusca", "Matrix")) %dopar% { # for(i in unique(membs)){ #!
              cat('cluster ', i, '\n')
              if(sum(membs == i) == 0) next
              n_indexes <- membs_df[which(membs_df[[column]] == i), 'sample_id']
              m_indexes <- membs_df[which(membs_df[[column]] != i), 'sample_id']
              n <- log(1 + Matrix::rowMeans(expm1(expDat[, n_indexes, drop=FALSE])), base=2)
              m <- log(1 + Matrix::rowMeans(expm1(expDat[, m_indexes, drop=FALSE])), base=2)
              names(n) <- rownames(expDat)
              names(m) <- rownames(expDat)
              pct.1 <- round(Matrix::rowSums(x=expDat[, n_indexes, drop=FALSE] > 0) / length(x=n_indexes), digits=3)
              pct.2 <- round(Matrix::rowSums(x=expDat[, m_indexes, drop=FALSE] > 0) / length(x=m_indexes), digits=3)
              d <- data.frame(mean.np=m, mean.p=n, fc=n-m, pct.1=pct.1, pct.2=pct.2)
              pct.min <- pmax(d$pct.1, d$pct.2)
              names(pct.min) <- rownames(d)
              genes <- names(which(pct.min >= min.pct))
              if(length(genes) == 0){
                print("No genes passed min.pct threshold")
              }else{
                d <- d[genes, ]
              }
              if(pos.only){
                genes.use <- rownames(d[which(d$fc >= fc.threshold), ])
              }else{
                genes.use <- rownames(d[which(abs(d$fc) >= fc.threshold), ])
              }
              if(length(genes.use) != 0){
                m <- m[genes.use]
                n <- n[genes.use]
                #print(length(genes.use))
                coldata <- slot(object, 'assays')[[assay.type]]@sampTab
                coldata[n_indexes, "group"] <- "Group1"
                coldata[m_indexes, "group"] <- "Group2"
                coldata$group <- factor(x=coldata$group)
                countdata.test <- expDat[genes.use, rownames(x=coldata),
                                         drop=FALSE]
                p_val <- sapply(X=1:nrow(x=countdata.test), FUN=function(x) {
                  return(wilcox.test(countdata.test[x, ] ~ coldata$group)$p.value)
                })
                d2 <- data.frame(mean.np=m, mean.p=n, fc=n-m, pv=p_val,
                                 p.adj=p.adjust(p_val, method='bonferroni'))
                #d2 <- d2[which(d2$pv < 0.01),]
                #print(nrow(d2))
                diffs[[i]] <- d2
              }else{
                print('No genes passed fold change threshold')
              }
            }
            names(diffs) = unique(membs); #!
            object@signatures <- diffs
            parallel::stopCluster(cl) #!
            return(object)
          }
)



#' Helper function of findSignatures.
#'
#' Find markers for each cell population. Create a table of genes,
#' fc_subpopulation, subpopulation with max expression.
#'
#' @param object CellRouter object.
#'
#' @return data frame;
setGeneric("findmarkers", function(object) standardGeneric("findmarkers"))
setMethod("findmarkers",
          signature="CellRouter",
          definition=function(object){
            print('finding subpopulation markers')
            # Select most representative gentes from each group of cells.
            genes <- unique(as.vector(unlist(lapply(object@signatures,
                                                    rownames))))
            df <- data.frame(matrix(0, nrow=length(genes),
                                    ncol=length(object@signatures)))
            rownames(df) <- genes
            colnames(df) <- names(object@signatures)
            # Attribute fold changes of genes in each group to the data frame.
            for(gene in rownames(df)){
              for(pop in colnames(df)){
                if(gene %in% rownames(object@signatures[[pop]])){
                  df[gene, pop] <- object@signatures[[pop]][gene, 'fc']
                }else{
                  df[gene, pop] <- 0
                }
              }
            }
            # Get the cell group which the genes have the greatest fold changes.
            x <- apply(df, 1, function(x){names(x)[which(x == max(x))][1]})
            # Get the values of maximum fold changes for each gene.
            xx <- apply(df, 1, function(x){max(x)})
            df$population <- x
            df$fc <- xx
            # Use just the population and fold change columns.
            df <- df[, c('population', 'fc')]
            # Get other values for the most expressed genens in each population.
            for(pop in names(object@signatures)){
              dfx <- df[which(df$population == pop),]
              df[rownames(dfx), 'mean.p'] <- object@signatures[[pop]][rownames(dfx), 'mean.p']
              df[rownames(dfx), 'mean.np'] <- object@signatures[[pop]][rownames(dfx), 'mean.np']
              df[rownames(dfx), 'pval'] <- object@signatures[[pop]][rownames(dfx), 'pv']
              df[rownames(dfx), 'p.adj'] <- object@signatures[[pop]][rownames(dfx), 'p.adj']
            }
            df$gene <- rownames(df)
            # Uses just positive values of fold changes.
            #
            # duvida: Não tinha que ser só os positivos se fosse up regulated genes?
            #
            df <- df[which(df$fc > 0),] #new line...
            return(df)
          }
)



#' Helper function of findSignatures.
#'
#' Find gene signatures based on a template-matching approach.
#' Find genes that are preferentially expressed in specified samples. Find it by
#' adjusting a line to the gene expression for each cell population, the getting
#' the most representative genes, the ones with better fit.
#'
#' @param expDat Expression matrix
#' @param sampTab Sample annotation table
#' @param qtile numeric; quantile (between 0 and 1) to select top genes
#' correlated with the idealized expression pattern.
#' @param remove logical; remove overlaping genes.
#' @param dLevel character; annotation level to group on, groups to compare.
#'
#' @return ; specificSets.
ranked_findSpecGenes<-function# find genes that are preferentially expressed in specified samples
(expDat, ### expression matrix
 sampTab, ### sample table
 qtile=0.95, ### quantile
 remove=FALSE,
 dLevel="population_name" #### annotation level to group on
){
  # aqui
  expDat <- as.data.frame(as.matrix(expDat))
  cat("Template matching...\n")
  myPatternG<-cn_sampR_to_pattern(as.vector(sampTab[,dLevel]));
  specificSets<-apply(myPatternG, 1, cn_testPattern, expDat=expDat);

  # Adaptively extract the best genes per lineage.
  cat("First pass identification of specific gene sets...\n")
  cvalT<-vector();
  ctGenes<-list();
  ctNames<-unique(as.vector(sampTab[,dLevel]));
  for(ctName in ctNames){
    x<-specificSets[[ctName]];
    # Selects the most representative genes according to the r squared value of
    # their fit regarding a gene population. The greater the r squared the
    # greater the fit.
    cval<-quantile(x$cval, qtile, na.rm=TRUE);
    tmp<-rownames(x[x$cval>cval,]);
    specificSets[[ctName]] <- specificSets[[ctName]][tmp,]
    # Most representative genes in a cell population.
    ctGenes[[ctName]]<-tmp;
    # Quantile value for that cell population.
    cvalT<-append(cvalT, cval);
  }
  if(remove){
    cat("Remove common genes...\n");
    # Limit to genes exclusive to each list.
    specGenes<-list();
    for(ctName in ctNames){
      # Other cell populations.
      others<-setdiff(ctNames, ctName);
      # Remove genes that are also representative to other cell populations.
      x<-setdiff( ctGenes[[ctName]], unlist(ctGenes[others]));
      specificSets[[ctName]] <- specificSets[[ctName]][x,]
      specificSets[[ctName]]$gene <- rownames(specificSets[[ctName]])
      specGenes[[ctName]]<-x;
    }
    result <- specGenes
  }else {
    result <- ctGenes;
  }
  specificSets <- lapply(specificSets, function(x){x[order(x$cval, decreasing=TRUE),]})
  specificSets <- lapply(specificSets, function(x){colnames(x) <- c('tm.pval', 'cval', 'tm.padj', 'gene'); x})

  specificSets
}



#' Helper function of ranked_findSpecGenes.
#'
#' Return a pattern for use in cn_testPattern (template matching). Convert the
#' cell population to matrix of 0 and 1, where the rows are the unique cell
#' populations and the columns are the population of each cell. Thus, the cell
#' values indicates the population of a cell.
#'
#' @param sampR character vector; vector of cell populations.
#'
#' @return dgCMatrix matrix; ans_sparse_matrix.
cn_sampR_to_pattern<-function# return a pattern for use in cn_testPattern (template matching)
(sampR){
  d_ids<-unique(as.vector(sampR));
  nnnc<-length(sampR);
  ans<-matrix(nrow=length(d_ids), ncol=nnnc);
  for(i in seq(length(d_ids))){
    x<-rep(0,nnnc);
    x[which(sampR==d_ids[i])]<-1;
    ans[i,]<-x;
  }
  colnames(ans)<-as.vector(sampR);
  rownames(ans)<-d_ids;
  ans;
}



#' Helper function of ranked_findSpecGenes.
#'
#' Get information about the line that represent each gene expression in each
#' cell population.
#'
#' @param pattern pattern
#' @param expDat expression matrix
#'
#' @return data frame; row.names=geneids, pval=pval, cval=cval, holm=holm.
cn_testPattern<-function(pattern, expDat){
  pval<-vector();
  cval<-vector();
  geneids<-rownames(expDat);
  # Calculate the lines that fits gene expressions for the cell group
  # represented by the pattern. Given a pattern of 0 and 1 representing the cell
  # population, it calculate the line of the gene expression for that cell type.
  # The y is the gene expression, the x is 1 if that cell belongs to the
  # population and 0 otherwise. The result is a list of matrix with the Estimate
  # of the intercept and the slope value, with Std.Err and t-value for each gene.
  llfit<-ls.print(lsfit(pattern, t(expDat)), digits=25, print=FALSE);
  xxx<-matrix( unlist(llfit$coef), ncol=8,byrow=TRUE);
  # t-value for X.
  ccorr<-xxx[,6];
  # R Squared.
  cval<- sqrt(as.numeric(llfit$summary[,2])) * sign(ccorr);
  # Pr(>|t|) for X.
  pval<-as.numeric(xxx[,8]);
  # Adjust p-values.
  holm<-p.adjust(pval, method='holm');
  data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm);
}
