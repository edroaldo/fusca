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
                                        group1, group2, fc.threshold = 0.25)
  standardGeneric("compareTwoGroups"))
#' @rdname compareTwoGroups-methods
#' @aliases compareTwoGroups
setMethod("compareTwoGroups",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                column = 'population', group1, group2,
                                fc.threshold){
            print('discovering subpopulation-specific gene signatures')
            expDat <- slot(object, 'assays')[[assay.type]]@ndata
            membs <- as.vector(slot(object, 'assays')[[assay.type]]@sampTab[[column]])
            diffs <- list()
            # m <- if(sum(membs == group2) > 1) apply(expDat[, membs == group2], 1, mean) else expDat[, membs == group2]
            # n <- if(sum(membs == group1) > 1) apply(expDat[, membs == group1], 1, mean) else expDat[, membs == group1]
            m <- if(sum(membs == group2) > 1) apply1_sp(expDat[, membs == group2], sum)/ncol(expDat[, membs == group2]) else expDat[, membs == group2]
            n <- if(sum(membs == group1) > 1) apply1_sp(expDat[, membs == group1], sum)/ncol(expDat[, membs == group1]) else expDat[, membs == group1]
            names(m) <- rownames(expDat)
            names(n) <- rownames(expDat)
            # Log scale
            d <- data.frame(mean.np=m, mean.p=n, fc=n-m)
            genes.use <- rownames(d[which(abs(d$fc) > fc.threshold),])
            if(length(genes.use) != 0){
              m <- m[genes.use]
              n <- n[genes.use]
              print(length(genes.use))
              # For wilcox test.
              coldata <- slot(object, 'assays')[[assay.type]]@sampTab
              coldata[membs == group1, "group"] <- "Group1"
              coldata[membs == group2, "group"] <- "Group2"
              coldata$group <- factor(x = coldata$group)
              countdata.test <- expDat[genes.use, rownames(x = coldata)]
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
