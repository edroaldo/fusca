
# Find variable genes -----------------------------------------------------

#' Find Variable Genes
#'
#' Preprocessing step to find variable genes to use in the downstream analysis.
#' It reduces the number of genes to be considered and speeds up the workflow.
#' It also plots the genes by their p-value.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.]
#' @param method character; the method to perform the identification of variable
#' genes: coefficient_variation or vst (more efficient).
#' @param pvalue numeric; the p-value threshold for the genes to be considered.
#' Argument for the coefficient_variation method.
#' @param loess.span numeric; argument for the vst method.
#'
#'
#' @return dataframe; the information about the genes.
#'
#' @import ggplot2
#' @import Rcpp
#' @import RcppEigen
#' @export
FindVariableGenes <- function(object, assay.type='RNA', method='vst',
                              pvalue=1e-6, loess.span=0.3){
  if (method == 'coefficient_variation') {
    ndata <- slot(object, 'assays')[[assay.type]]@rawdata
    ndata <- Matrix::t(Matrix::t(ndata)/apply1_sp(Matrix::t(ndata), sum))
    means <- Matrix::rowMeans(ndata)
    vars <- apply(ndata, 1, var)
    cv2 <- vars/means^2
    minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
    useForFit <- means >= minMeanForFit # & spikeins
    fit <- statmod::glmgam.fit( cbind(a0 = 1, a1tilde = 1/means[useForFit] ),
                                cv2[useForFit] )
    a0 <- unname( fit$coefficients["a0"] )
    a1 <- unname( fit$coefficients["a1tilde"])
    #fit$coefficients
    #double-check...
    ###xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
    ##vfit <- a1/xg + a0
    df <- ncol(ndata) - 1
    afit <- a1/means+a0
    varFitRatio <- vars/(afit*means^2)
    varorder <- order(varFitRatio, decreasing=T)
    #oed <- slot(object, 'assays')[[assay.type]]@rawdata[varorder,]
    pval <- pchisq(varFitRatio*df, df=df, lower.tail=F)
    adj.pval <- p.adjust(pval, "fdr")
    var.genes <- names(adj.pval[which(adj.pval < pvalue)])
    #object@var.genes <- var.genes
    vargenes <- data.frame(mean=log(means), cv2=log(cv2), fit=afit,
                           adj.pvalue=adj.pval, size=-log(adj.pval))
    vargenes$gene <- rownames(vargenes)
    vargenes$var.genes <- 'no'
    vargenes[as.vector(vargenes$gene) %in% var.genes, 'var.genes'] <- 'yes'
    g <- ggplot(vargenes, aes(mean, cv2)) + geom_point(aes(colour=var.genes)) +
      theme_bw() + xlab("log(mean)") + ylab("log(CV2)") +
      #geom_line(aes(fit)) +
      theme(axis.text.x=element_text(size=12, angle=45, hjust=1),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border=element_rect(fill = NA, colour=alpha('black', 1), size=1)) +
      scale_color_brewer("", palette = 'Paired')
    print(g)
    return(vargenes)
  } else {
    ## Cpp code
    print('warnings supressed')
    suppressMessages(sourceCpp(code='
    #include <Rcpp.h>
    #include <RcppCommon.h>
    #include <RcppEigen.h>
    #include <iostream>
    #include <cmath>
    #include <unordered_map>
    #include <fstream>
    #include <string>
    #include <Rinternals.h>
    // [[Rcpp::depends(RcppEigen)]]
    using namespace std;
    using namespace Eigen;
    using namespace Rcpp;
    typedef Eigen::ArrayXd MapAr1;
    // [[Rcpp::export]]
     NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat,
                                  NumericVector mu,
                                  NumericVector sd,
                                  double vmax){
      mat = mat.transpose();
      NumericVector allVars(mat.cols());
      for (int k=0; k<mat.outerSize(); ++k){
        if (sd[k] == 0) continue;
        double colSum = 0;
        int nZero = mat.rows();
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          nZero -= 1;
          colSum += pow(std::min(vmax, (it.value() - mu[k]) / sd[k]), 2);
        }
        colSum += pow((0 - mu[k]) / sd[k], 2) * nZero;
        allVars[k] = colSum / (mat.rows() - 1);
      }
      return(allVars);
    }', echo=FALSE, showOutput=FALSE, verbose=FALSE))
    print('finished warnings?')
    ### R code
    object_group <- slot(object, 'assays')[[assay.type]]@rawdata
    clip.max <- sqrt(ncol(object_group))
    hvf.info <- data.frame(mean = Matrix::rowMeans(object_group))
    hvf.info$variance <- sparseMatrixStats::rowVars(object_group)
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
    ## loess curve fit
    fit <- loess(formula = log10(variance) ~ log10(mean),
                 data = hvf.info[not.const, ], span = loess.span)
    ## extract fitted variance
    hvf.info$variance.expected[not.const] <- 10^fit$fitted
    hvf.info$variance.standardized <- SparseRowVarStd(
      mat = object_group,
      mu = hvf.info$mean,
      sd = sqrt(hvf.info$variance.expected),
      vmax = clip.max
    )
    hvf.info <- hvf.info[order(hvf.info$variance.standardized, decreasing = TRUE),]
    return(hvf.info)
  }
}




# Scale Data --------------------------------------------------------------

#' Scale data
#'
#' Scale and center the data. Individually regress variables provided in
#' vars.regress using using a linear model. Other models are under development.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genes.use vector; genes to scale/center. Default is all genes in
#' normalized data.
#' @param scale.max numeric; max value in scaled data.
#' @param vars.regress character vector; variables to regress out.
#' @param blocksize numeric; size of the blocks in which genes will be scaled.
#'
#' @return CellRouter object with the scale.data slot updated.
#'
#' @export
#' @docType methods
#' @rdname scaleData-methods
setGeneric("scaleData", function(object, assay.type='RNA',
                                 genes.use = NULL, scale.max=10,
                                 vars.regress = NULL, blocksize=500)
  standardGeneric("scaleData"))
#' @rdname scaleData-methods
#' @aliases scaleData
setMethod("scaleData",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                genes.use, scale.max, vars.regress,
                                blocksize){
            # Select genes to scale.
            if(is.null(genes.use)){
              genes.use <- rownames(slot(object, 'assays')[[assay.type]]@ndata)
            }
            data.use <- slot(object, 'assays')[[assay.type]]@ndata[genes.use, ]
            if(!is.null(vars.regress)){
              print('Regression...')
              cat(vars.regress, '\n')
              gene.expr <-
                slot(object, 'assays')[[assay.type]]@ndata[genes.use, , drop = FALSE]
              latent.data <- as.data.frame(
                slot(object, 'assays')[[assay.type]]@sampTab[, vars.regress])
              rownames(latent.data) <- rownames(
                slot(object, 'assays')[[assay.type]]@sampTab
              )
              colnames(latent.data) <- vars.regress
              new.data <- sapply(X = genes.use, FUN=function(x){
                regression.mat <- cbind(latent.data, gene.expr[x,])
                colnames(regression.mat) <- c(colnames(latent.data), "GENE")
                fmla <- as.formula(
                  object = paste0(
                    "GENE ",
                    " ~ ",
                    paste(vars.regress, collapse = "+")
                  )
                )
                lm(formula = fmla, data = regression.mat)$residuals
              })
              data.use <- Matrix::t(new.data)
            }
            # scale.data <- Matrix::t(scale(Matrix::t(data.use),
            #                               center = rep(0, dim(data.use)[1])))
            #scale.data <- Matrix::t(scale.sparseMatrix(Matrix::t(data.use)))
            #scale.data[is.na(scale.data)] <- 0
            #scale.data[Matrix::which(scale.data > scale.max)] <- scale.max
            #slot(object, 'assays')[[assay.type]]@scale.data <- scale.data
            # blocksize = 25 #blocks of 25 genes
            #chunks <- split(rownames(data.use),
            #                ceiling(seq_along(rownames(data.use))/blocksize))
            #tmp.data <- lapply(chunks, function(x){
            #  Matrix::t(scale(Matrix::t(data.use[x,])))
            #})
            #scale.data <- do.call(rbind, tmp.data)
            # scale.data <- Matrix::t(scale.sparse(Matrix::t(data.use)))
            scale.data <- Matrix::t(scale.sparse(Matrix::t(data.use), blocksize=blocksize))
            scale.data[is.na(scale.data)] <- 0
            scale.data[Matrix::which(abs(scale.data) > scale.max)] <- scale.max
            slot(object, 'assays')[[assay.type]]@scale.data <- scale.data
            gc(verbose = FALSE)
            object
          }
)

