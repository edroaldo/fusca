#' Show gene expression in the space of reduced dimensionality.
#'
#' Plot selected gene expression in reduced dimension space and save to file.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genelist character vector; genes to show.
#' @param reduction.type character; the dimension reduction space to be used:
#' pca, tsne, DC of diffusion components, umap, or custom.
#' @param threshold numeric; threshold to rescale gene expression.
#' @param dims.use numeric vector; dimensions to use.
#' @param columns numeric; number of columns in the output figure.
#' @param dotsize numeric; dot size.
#' @param alpha numeric; transparency (between 0 and 1).
#' @param title character; the title of the plot.
#'
#' @return ggplot2; plot.
#'
#' @export
#' @docType methods
#' @rdname plotDRExpression-methods
setGeneric("plotDRExpression", function(object, assay.type='RNA', genelist,
                                        reduction.type = c('tsne', 'pca', 'DC',
                                                           'umap', 'custom'),
                                        threshold = 2, dims.use = c(1, 2),
                                        columns = 5, dotsize = 1, alpha = 0.5,
                                        title)
  standardGeneric("plotDRExpression"))
#' @rdname plotDRExpression-methods
#' @aliases plotDRExpression
setMethod("plotDRExpression",
          signature="CellRouter",
          definition=function(object, assay.type, genelist,
                              reduction.type = c('tsne', 'pca', 'DC', 'umap',
                                                 'custom'),
                              threshold = 2, dims.use = c(1, 2), columns = 5,
                              dotsize, alpha, title){
            reduction.type <- match.arg(reduction.type)
            matrix <- as.data.frame(slot(object, reduction.type)$
                                      cell.embeddings[ , dims.use])
            plots <- list()
            scores <- matrix
            colnames(scores) <- c('Dim_1', 'Dim_2')
            if(reduction.type == 'tsne'){
              xlab <- paste('tSNE ', dims.use[1], sep=' ')
              ylab <- paste('tSNE ', dims.use[2], sep=' ')
            } else if (reduction.type == 'pca'){
              xlab <- paste('PC', dims.use[1], sep='')
              ylab <- paste('PC', dims.use[2], sep='')
            } else if (reduction.type == 'DC'){
              xlab <- paste('DC', dims.use[1], sep='')
              ylab <- paste('DC', dims.use[2], sep='')
            }else if(reduction.type == 'umap'){
              xlab <- paste('UMAP', dims.use[1], sep='')
              ylab <- paste('UMAP', dims.use[2], sep='')
            } else {
              xlab <- 'Dim 1'
              ylab <- 'Dim 2'
            }
            # x <- slot(object, 'assays')[[assay.type]]@ndata[genelist, rownames(matrix)]
            x <- center_with_threshold(
              slot(object, 'assays')[[assay.type]]@ndata[genelist, rownames(matrix),
                                                         drop=FALSE],
              threshold)
            dfs <- data.frame()
            for(gene in genelist){
              expr <- x[gene, ]
              scores$GENE <- as.numeric(expr)
              scores$gene <- gene
              dfs <- rbind(dfs, scores)
            }
            dfs <- dfs[order(dfs$GENE), ]
            dfs$gene <- factor(dfs$gene, levels = genelist)
            # Create plot.
            p1 <- ggplot2::ggplot(dfs, ggplot2::aes(x = Dim_1, y = Dim_2,
                                                    colour = GENE)) +
              ggplot2::geom_point(ggplot2::aes(alpha = GENE), size = dotsize) +
              ggplot2::theme_bw() +
              ggplot2::scale_colour_gradientn("Relative expression",
                                              colours = c("midnightblue",
                                                          "white",
                                                          "orange")) +
              ggplot2::ylab(ylab) + ggplot2::xlab(xlab) +
              ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA,
                                                                  colour = "white"),
                             strip.background = ggplot2::element_blank()) +
              ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             panel.border = ggplot2::element_blank(),
                             panel.background = ggplot2::element_blank()) +
              ggplot2::theme(legend.position="none",
                             strip.background = ggplot2::element_rect(colour = "white",
                                                                      fill = "white")) +
              ggplot2::guides(colour = ggplot2::guide_colourbar(title.position = "top",
                                                                title.hjust = 0.5),
                              size = ggplot2::guide_legend(title.position = "top",
                                                           title.hjust = 0.5)) +
              ggplot2::facet_wrap(~gene, ncol = columns) +
              ggplot2::ggtitle(title)
            return(p1)
          }
)



#' Show gene expression in the space of reduced dimensionality.
#'
#' Plot selected gene expression in reduced dimension space and save to file.
#' Version with no title.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genelist character vector; genes to show.
#' @param reduction.type character; the dimension reduction space to be used:
#' pca, tsne, DC of diffusion components, umap, or custom.
#' @param threshold numeric; threshold to rescale gene expression.
#' @param dims.use numeric vector; dimensions to use.
#' @param columns numeric; number of columns in the output figure.
#' @param dotsize numeric; dot size.
#' @param alpha numeric; transparency (between 0 and 1).
#'
#' @return ggplot2; plot.
#'
#' @export
#' @docType methods
#' @rdname plotDRExpression2-methods
setGeneric("plotDRExpression2", function(object, assay.type='RNA', genelist,
                                         reduction.type = c('tsne', 'pca', 'DC',
                                                            'umap', 'custom'),
                                         threshold = 2, dims.use = c(1,2),
                                         columns = 5, dotsize = 1, alpha = 0.5)
  standardGeneric("plotDRExpression2"))
#' @rdname plotDRExpression2-methods
#' @aliases plotDRExpression2
setMethod("plotDRExpression2",
          signature="CellRouter",
          definition=function(object, assay.type, genelist,
                              reduction.type = c('tsne', 'pca', 'DC', 'umap',
                                                 'custom'),
                              threshold = 2, dims.use = c(1,2), columns = 5,
                              dotsize = 1, alpha = 0.5){
            reduction.type <- match.arg(reduction.type)
            matrix <- as.data.frame(slot(object, reduction.type)$
                                      cell.embeddings[,dims.use])
            plots <- list()
            scores <- matrix
            colnames(scores) <- c('Dim_1', 'Dim_2')
            if(reduction.type == 'tsne'){
              xlab <- paste('tSNE ', dims.use[1], sep=' ')
              ylab <- paste('tSNE ', dims.use[2], sep=' ')
            } else if (reduction.type == 'pca'){
              xlab <- paste('PC', dims.use[1], sep='')
              ylab <- paste('PC', dims.use[2], sep='')
            } else if (reduction.type == 'DC'){
              xlab <- paste('DC', dims.use[1], sep='')
              ylab <- paste('DC', dims.use[2], sep='')
            }else if(reduction.type == 'umap'){
              xlab <- paste('UMAP', dims.use[1], sep='')
              ylab <- paste('UMAP', dims.use[2], sep='')
            } else {
              xlab <- 'Dim 1'
              ylab <- 'Dim 2'
            }
            x <- slot(object, 'assays')[[assay.type]]@ndata[genelist, rownames(matrix)]
            gc()
            dfs <- data.frame()
            for(gene in genelist){
              expr <- x[gene,]
              scores$GENE <- as.numeric(expr)
              scores$gene <- gene
              dfs <- rbind(dfs, scores)
            }
            dfs <- dfs[order(dfs$GENE),]
            dfs$gene <- factor(dfs$gene, levels=genelist)
            # Create plot.
            p1 <- ggplot2::ggplot(dfs,ggplot2::aes(x = Dim_1, y = Dim_2,
                                                   colour = GENE)) +
              ggplot2::geom_point(ggplot2::aes(alpha = GENE), size = dotsize) +
              ggplot2::theme_bw() +
              ggplot2::scale_colour_gradientn("Relative expression",
                                              colours = c("midnightblue",
                                                          "white",
                                                          "orange")) +
              ggplot2::ylab(ylab) + ggplot2::xlab(xlab) +
              ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA,
                                                                  colour = "white"),
                             strip.background = ggplot2::element_blank()) +
              ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             panel.border = ggplot2::element_blank(),
                             panel.background = ggplot2::element_blank()) +
              ggplot2::theme(legend.position="bottom",
                             strip.background = ggplot2::element_rect(colour = "white",
                                                                      fill = "white")) +
              ggplot2::guides(colour = ggplot2::guide_colourbar(title.position = "top",
                                                                title.hjust = 0.5),
                              size = ggplot2::guide_legend(title.position = "top",
                                                           title.hjust = 0.5)) +
              ggplot2::facet_wrap(~gene, ncol = columns)
            return(p1)
          }
)
