#' Dot plot of genes in populations.
#'
#' Present a dot plot of selected genes in the specified populations and save to
#' file. The dot size represents the percentage of cells expressing the
#' indicated genes in each cluster and the dot color intensity represents the
#' average expression level of the indicated genes.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genes.use character vector; selected genes.
#' @param column character; column name in cellrouter@@sampTab indicating the
#' population.
#' @param thr numeric; threshold.
#' @param order character vector; ordered populations to be considered.
#'
#' @return ggplot2; plot.
#'
#' @export
#' @docType methods
#' @rdname dotplot-methods
setGeneric("dotplot", function(object, assay.type='RNA',
                               genes.use, column='population', thr,
                               order)
  standardGeneric("dotplot"))
#' @rdname dotplot-methods
#' @aliases dotplot
setMethod("dotplot",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                genes.use, column, thr, order){
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            perc <- data.frame(matrix(0, nrow=length(genes.use), ncol=0))
            exp <- perc
            rownames(perc) <- genes.use
            for(i in unique(sampTab[[column]])){
              cells.population <- rownames(sampTab[which(sampTab[[column]] == i),])
              p <- apply(slot(object, 'assays')[[assay.type]]@ndata[genes.use, cells.population],
                         1, function(x){sum(x>thr)/length(x)})
              perc <- cbind(perc, p)
              v <- apply(slot(object, 'assays')[[assay.type]]@ndata[genes.use, cells.population],
                         1, mean)
              exp <- cbind(exp, v)
            }
            colnames(exp) <- unique(sampTab[[column]])
            hc <- hclust(dist(exp, method='euclidean'), method = 'complete')
            rc <- hclust(dist(t(exp), method='euclidean'), method = 'complete')
            exp2 <- exp[hc$order, rc$order]
            order.c <- colnames(exp2)
            colnames(perc) <- unique(sampTab[[column]])
            colnames(exp) <- unique(sampTab[[column]])
            rownames(exp) <- sapply(strsplit(rownames(exp), split='__', fixed=TRUE), function(x){x[1]})
            perc$gene <- rownames(perc)
            perc = reshape2::melt(perc, id.vars  = 'gene')
            exp$gene <- rownames(exp)
            exp$gene <- factor(exp$gene, levels=genes.use)
            exp <- reshape2::melt(exp, id.vars='gene')
            exp$Percentage <- perc$value*100
            exp$variable <- factor(exp$variable, levels=order)
            g <- ggplot2::ggplot(exp, ggplot2::aes(gene, variable)) +
              ggplot2::geom_point(ggplot2::aes(fill=value, size=Percentage),
                                  colour="black", pch=21) +
              ggplot2::theme_bw() + ggplot2::xlab("") + ggplot2::ylab("") +
              ggplot2::theme(axis.text.x = ggplot2::element_text(size=10,
                                                                 angle=45, hjust=1),
                             legend.position = 'right',
                    panel.grid.minor = ggplot2::element_blank(),
                    panel.border = ggplot2::element_rect(fill = NA,
                                                         colour = ggplot2::alpha('black', 0.5),
                                                         size=1)) +
              ggplot2::scale_fill_gradientn("Average expression",
                                            colours = c("white","goldenrod1","brown")) +
              ggplot2::guides(col = ggplot2::guide_legend(direction="vertical",
                                                          keywidth = 0.75,
                                                          keyheight = 0.75,
                                                          override.aes = list(size=3))) +
              ggplot2::coord_flip()
            return(g)
          }
)


#' Dot plot of genes in populations.
#'
#' Present a dot plot of selected genes in the specified populations in the
#' column argument by the ones in the column2 argument and save to file. The
#' dot size represents the percentage of cells expressing the indicated genes in
#' each cluster and the dot color intensity represents the average expression
#' level of the indicated genes.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genes.use character vector; genes to show.
#' @param column character; column name in cellrouter@@sampTab indicating the
#' population.
#' @param column2 character; column name in cellrouter@@sampTab indicating the
#' groups that the population form column belongs.
#' @param thr numeric; threshold.
#' @param order character vector; ordered populations from column to be
#' considered.
#' @param order2 character vector; ordered populations from column2 to be
#' considered.
#' @param cols numeric; number of columns in the output figure.
#'
#' @return ggplot2; plot.
#'
#' @export
#' @docType methods
#' @rdname dotplotByMetadata-methods
setGeneric("dotplotByMetadata", function(object, assay.type='RNA',
                                         genes.use, column='population',
                                         column2, thr, order, order2, cols)
  standardGeneric("dotplotByMetadata"))
#' @rdname dotplotByMetadata-methods
#' @aliases dotplotByMetadata
setMethod("dotplotByMetadata",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                genes.use, column, column2, thr, order, cols){
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            combined <- data.frame()
            for(c in unique(sampTab[[column2]])){
              perc <- data.frame(matrix(0, nrow=length(genes.use), ncol=0))
              exp <- perc
              rownames(perc) <- genes.use
              xsampTab <- sampTab[which(sampTab[[column2]] == c),]
              for(i in unique(xsampTab[[column]])){
                cells.population <- rownames(xsampTab[which(xsampTab[[column]] == i),])
                p <- apply(slot(object, 'assays')[[assay.type]]@ndata[genes.use, cells.population],
                           1, function(x){sum(x>thr)/length(x)})
                perc <- cbind(perc, p)
                v <- apply(slot(object, 'assays')[[assay.type]]@ndata[genes.use, cells.population],
                           1, mean)
                exp <- cbind(exp, v)
              }
              colnames(exp) <- unique(xsampTab[[column]])
              colnames(perc) <- unique(xsampTab[[column]])
              colnames(exp) <- unique(xsampTab[[column]])
              rownames(exp) <- sapply(strsplit(rownames(exp), split='__', fixed=TRUE), function(x){x[1]})
              perc$gene <- rownames(perc)
              perc = reshape2::melt(perc, id.vars  = 'gene')
              exp$gene <- rownames(exp)
              exp$gene <- factor(exp$gene, levels=genes.use)
              exp <- reshape2::melt(exp, id.vars='gene')
              exp$Percentage <- perc$value*100
              exp$variable <- factor(exp$variable, levels=order)
              exp$facet <- c
              combined <- rbind(exp, combined)
            }
            combined$facet <- factor(combined$facet, levels=order2)
            g <- ggplot2::ggplot(combined, ggplot2::aes(variable, gene)) +
              ggplot2::geom_point(ggplot2::aes(fill=value, size=Percentage),
                                  colour="black",pch=21) +
              ggplot2::theme_bw() + ggplot2::xlab("") + ggplot2::ylab("") +
              ggplot2::theme(axis.text.x=ggplot2::element_text(size=10, angle=45,
                                                               hjust=1),
                             legend.position = 'right',
                    panel.grid.minor = ggplot2::element_blank(),
                    panel.border = ggplot2::element_rect(fill = NA, colour =
                                                           ggplot2::alpha('black', 0.5),
                                                         size=1),
                    strip.background = ggplot2::element_rect(colour="white", fill="white")) +
              ggplot2::scale_fill_gradientn("Average expression",
                                            colours=c("white","goldenrod1","brown")) +
              ggplot2::facet_wrap(~facet, ncol = cols, strip.position = "top", scales = "free_y") #+ coord_flip()
            ggplot2::guides(col = ggplot2::guide_legend(direction="vertical",
                                                        keywidth = 0.75,
                                                        keyheight = 0.75,
                                                        override.aes = list(size=3))) #+ coord_flip()
            return(g)
          }
)
