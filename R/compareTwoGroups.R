#' Compute fold changes and find gene signatures.
#'
#' Compute mean expression of each gene in two populations and compares them
#' using fold changes, p-value of wilcox test and adjusted p-value using the
#' bonferroni method.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param column character; column in the metadata table to group cells for
#' differential expression. For example, if 'population' is specified,
#' population-specific gene signatures will be identified.
#' @param group1 character; name of group 1 used in comparison.
#' @param group2 character; name of group 2 used in comparison.
#' @param fc.threshold numeric; fold change threshold.
#' @param min.pct numeric; minimum percentage.
#'
#' @return data frame; a data frame with the mean expression of the genes in
#' group 1 (p column), the mean expression of the genes in group 2 (np column),
#' the fold changes (fc column), the p-value for the wilcox test (pv column),
#' and the adjusted p-value (p.adj column).
#'
#' @export
#' @docType methods
#' @rdname compareTwoGroups-methods
setGeneric("compareTwoGroups", function(object, assay.type='RNA',
                                        column = 'population',
                                        group1, group2, fc.threshold = 0.25, min.pct=0.1)
  standardGeneric("compareTwoGroups"))
#' @rdname compareTwoGroups-methods
#' @aliases compareTwoGroups
setMethod("compareTwoGroups",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                column = 'population', group1, group2,
                                fc.threshold, min.pct){
            print('discovering subpopulation-specific gene signatures')
            expDat <- slot(object, 'assays')[[assay.type]]@ndata
            membs <- as.vector(slot(object, 'assays')[[assay.type]]@sampTab[[column]])
            membs_df <- as.data.frame(slot(object, 'assays')[[assay.type]]@sampTab[ , c('sample_id', column), drop=FALSE])
            n_indexes <- membs_df[which(membs_df[[column]] == group1), 'sample_id']
            m_indexes <- membs_df[which(membs_df[[column]] == group2), 'sample_id']
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
            # m <- if(sum(membs == group2) > 1) apply(expDat[, membs == group2], 1, mean) else expDat[, membs == group2]
            # n <- if(sum(membs == group1) > 1) apply(expDat[, membs == group1], 1, mean) else expDat[, membs == group1]
            #m <- if(sum(membs == group2) > 1) apply1_sp(expDat[, membs == group2], sum)/ncol(expDat[, membs == group2]) else expDat[, membs == group2]
            #n <- if(sum(membs == group1) > 1) apply1_sp(expDat[, membs == group1], sum)/ncol(expDat[, membs == group1]) else expDat[, membs == group1]
            #names(m) <- rownames(expDat)
            #names(n) <- rownames(expDat)
            # Log scale
            #d <- data.frame(mean.np=m, mean.p=n, fc=n-m)
            genes.use <- rownames(d[which(abs(d$fc) >= fc.threshold),])
            if(length(genes.use) != 0){
              m <- m[genes.use]
              n <- n[genes.use]
              print(length(genes.use))
              # For wilcox test.
              coldata <- slot(object, 'assays')[[assay.type]]@sampTab
              coldata[n_indexes, "group"] <- "Group1"
              coldata[m_indexes, "group"] <- "Group2"
              coldata$group <- factor(x = coldata$group)
              countdata.test <- expDat[genes.use, rownames(x = coldata),
                                       drop=FALSE]
              p_val <- sapply(X = 1:nrow(x = countdata.test), FUN = function(x) {
                return(wilcox.test(countdata.test[x, ] ~ coldata$group)$p.value)
              })
              d2 <- data.frame(mean.np=m, mean.p=n, fc=n-m, pv=p_val,
                               p.adj=p.adjust(p_val, method='bonferroni'))
            } else {
              print('no genes were selected for this fold change threshold')
              d2 <- NULL
            }
            return(d2)
          }
)
