#' Compute and plot heatmap of GRN scores for transcriptional regulators.
#'
#' Compute and plot heatmap of the GRN scores of most representative genes in
#' trajectories and save to file. Integrate gene regulatory networks with gene
#' expression dynamics along the trajectories to calculate a GRN score and
#' identify putative regulators of these cell-fate transitions.
#'
#' @param object CellRouter object.
#' @param ggrn igraph; the gene regulatory network.
#' @param tfs character vector; gene names of transcriptional regulators.
#' @param transitions character vector; selected transitions of interest.
#' @param direction character; plot genes up-regulated, down-regulated or both
#' along trajectories.
#' @param dir.targets character; whether the predicted targets are up or
#' down-regulated.
#' @param q.up numeric; cutoff to select top q.up transcriptional regulators.
#' @param q.down numeric; cutoff to select top q.down transcriptional regulators.
#' @param flip boolean; apply coordinate flip (horizontal and vertical) to the
#' plot.
#'
#' @return ggplot2; plot.
#'
#' @export
#' @docType methods
#' @rdname grnScoresHeatmap-methods
setGeneric('grnScoresHeatmap', function(object, ggrn, tfs, transitions,
                                        direction = c("up", "down", "both"),
                                        dir.targets = c("up", "down"),
                                        q.up = 0.95, q.down = 0.05, flip)
  standardGeneric('grnScoresHeatmap'))
#' @rdname grnScoresHeatmap-methods
#' @aliases grnScoresHeatmap
setMethod('grnScoresHeatmap',
          signature="CellRouter",
          definition=function(object, ggrn, tfs, transitions,
                              direction = c("up", "down", "both"),
                              dir.targets = c("up", "down"), q.up = 0.95,
                              q.down = 0.05, flip){
            direction <- match.arg(direction)
            dir.targets <- match.arg(dir.targets)
            plots <- list()
            allscores <- list()
            heatmap <- list()
            for(p in transitions){
              if(direction=='up'){
                tfs.transition <-
                  intersect(tfs, names(object@top.correlations[['up']][[p]]))
              }else if(direction == 'down'){
                tfs.transition <-
                  intersect(tfs, names(object@top.correlations[['down']][[p]]))
              }else{
                tfs.transition.up <-
                  intersect(tfs, names(object@top.correlations[['up']][[p]]))
                tfs.transition.down <-
                  intersect(tfs, names(object@top.correlations[['down']][[p]]))
                tfs.transition <- c(tfs.transition.up, tfs.transition.down)
              }
              tfs.transition <- intersect(igraph::V(ggrn)$name, tfs.transition)
              tfs.transition <- object@correlation[[p]][tfs.transition]
              averages <- vector()
              num.genes <- vector()
              names <- vector()
              tf.targets <- list()
              for(r in names(tfs.transition)){
                rgrn <- igraph::induced.subgraph(
                  ggrn,
                  vids = unlist(igraph::neighborhood(graph = ggrn, order = 1,
                                                     nodes = r)))
                x <- object@top.correlations[[dir.targets]][[p]]
                # Subnetwork active during transition p.
                genes <- intersect(igraph::V(rgrn)$name, names(x))
                corrs <- object@correlation[[p]][genes]
                if(length(corrs) == 0){
                  #cat(r, 'has no targets\n')
                  next
                }else if(length(corrs) == 1 & names(corrs) == r){
                  #cat(r, 'regulates only itself\n')
                  next
                }else if(length(corrs) > 0){
                  # At least one target required.
                  tf.targets[[r]] <- names(corrs)
                  averages <- append(averages, mean(corrs))
                  num.genes <- append(num.genes, length(corrs))
                  names <- append(names, r)
                }
              }
              names(averages) <- names
              names(num.genes) <- names
              aux <- averages[is.na(averages)]
              averages <- averages[!is.na(averages)]
              # Order.
              averages <- averages[order(averages, decreasing = TRUE)]
              num.genes <- num.genes[names(averages)]
              # Rescale num.genes.
              # When it rescaled, it changes the scores when only up or
              # down-regulated genes are included.
              num.genes <- scales::rescale(num.genes, newrange = c(0.01,1))
              averages <- abs(averages)

              scores <- tfs.transition[names(averages)] * averages * num.genes
              # Different part from gnrscores.
              scores <- scores[order(scores, decreasing = TRUE)]
              targets <- tf.targets[names(scores)]
              # If up or down, q.up or q.down are the top genes.
              # If both, q.up or q.down are quantiles.
              if(direction=='up'){
                scores <- scores[order(scores, decreasing=TRUE)]
                scores <- scores[1:q.up]
                # In the gnrscore there is the color selection too.
              }else if(direction=='down'){
                scores <- scores[order(scores, decreasing=FALSE)]
                scores <- scores[1:q.down]
              }else{
                scores <- scores[which(scores > quantile(scores, q.up) |
                                         scores < quantile(scores, q.down))]
              }
              allscores[[p]] <- list(scores=scores, targets=targets)
              heatmap[[p]] <- scores
            }
            allgenes <- unique(as.vector(unlist(lapply(heatmap, names))))
            heatmap.m <- data.frame(matrix(0, nrow = length(allgenes),
                                           ncol = length(heatmap)))
            rownames(heatmap.m) <- allgenes
            colnames(heatmap.m) <- names(heatmap)
            for(g in rownames(heatmap.m)){
              for(p in colnames(heatmap.m)){
                x <- heatmap[[p]][g]
                heatmap.m[g, p] <- x
              }
            }
            heatmap.m[is.na(heatmap.m)] <- 0
            matrix <- as.data.frame(heatmap.m)
            hc <- hclust(dist(matrix), method='ward.D')
            matrix <- matrix[hc$order,]
            matrix$gene <- rownames(matrix)
            matrix.m <- reshape2::melt(matrix, id.var="gene")
            matrix.m$gene <- factor(rownames(matrix),
                                    levels = rev(rownames(matrix)))
            # Plot
            g1 <- ggplot2::ggplot(matrix.m, ggplot2::aes(variable, gene)) +
              ggplot2::geom_tile(ggplot2::aes(fill = value)) +
              ggplot2::scale_fill_gradientn("GRN score",
                                            colours = c("midnightblue",
                                                        "dodgerblue3",
                                                        'white',
                                                        "goldenrod1",
                                                        "darkorange2")) +
              ggplot2::theme_bw() +
              ggplot2::xlab("Transitions") + ggplot2::ylab("") +
              ggplot2::theme(legend.position = "right",
                    axis.text.y = ggplot2::element_text(size = ggplot2::rel(1),
                                                        angle = 00),
                    axis.text.x = ggplot2::element_text(size = ggplot2::rel(1),
                                                        angle = 30, hjust = 1),
                    panel.grid.major = ggplot2::element_blank(),
                    panel.grid.minor = ggplot2::element_blank(),
                    panel.border = ggplot2::element_rect(fill = NA,
                                                         colour = ggplot2::alpha('black', 1),
                                                         size = 1)) +
              ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_discrete(expand = c(0, 0))
            # Flip
            if(flip){
              g1 <- g1 + ggplot2::coord_flip()
            }
            return(g1)
          }
)
