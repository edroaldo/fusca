#' Compute and plot GRN scores for transcriptional regulators.
#'
#' Compute and plot the GRN scores of most representative genes in trajectories
#' and save to file. Integrate gene regulatory networks with gene expression
#' dynamics along the trajectories to calculate a GRN score and identify
#' putative regulators of these cell-fate transitions.
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
#' @return list; all scores list with scores and targets, and plots.
#'
#' @export
#' @docType methods
#' @rdname grnscores-methods
setGeneric('grnscores', function(object, ggrn, tfs, transitions,
                                 direction = c("up", "down", "both"),
                                 dir.targets = c("up", "down"),
                                 q.up = 0.95, q.down = 0.05, flip)
  standardGeneric('grnscores'))
#' @rdname grnscores-methods
#' @aliases grnscores
setMethod('grnscores',
          signature="CellRouter",
          definition=function(object, ggrn, tfs, transitions,
                              direction = c("up", "down", "both"),
                              dir.targets = c("up", "down"),
                              q.up = 0.95, q.down = 0.05, flip){
            direction <- match.arg(direction)
            dir.targets <- match.arg(dir.targets)
            plots <- list()
            allscores <- list()
            for (p in transitions){
              if (direction=='up') {
                tfs.transition <-
                  intersect(tfs, names(object@top.correlations[['up']][[p]]))
              } else if (direction == 'down'){
                tfs.transition <-
                  intersect(tfs, names(object@top.correlations[['down']][[p]]))
              } else {
                tfs.transition.up <-
                  intersect(tfs, names(object@top.correlations[['up']][[p]]))
                tfs.transition.down <-
                  intersect(tfs,names(object@top.correlations[['down']][[p]]))
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
                }else if((length(corrs) == 1) & (identical(names(corrs), r))){
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
              # When it is rescaled, it changes the scores when only up or
              # down-regulated genes are included.
              num.genes <- scales::rescale(num.genes, newrange = c(0.01,1))
              averages <- abs(averages)
              scores <- tfs.transition[names(averages)] * averages * num.genes
              # If up or down, q.up or q.down are the top genes.
              # If both, q.up or q.down are quantiles.
              if(direction=='up'){
                # Order.
                scores <- scores[order(scores, decreasing = TRUE)]
                scores <- scores[1:q.up]
                colors <- c('white', 'orange', 'red')
              }else if(direction=='down'){
                # Order.
                scores <- scores[order(scores, decreasing = FALSE)]
                scores <- scores[1:q.down]
                colors <- c('blue', 'yellow', 'green')
              }else{
                scores <- scores[which(scores > quantile(scores, q.up) |
                                         scores < quantile(scores, q.down))]
                colors <- c('blue3', 'white', 'chartreuse3')
              }
              scores <- scores[order(scores, decreasing = TRUE)]
              targets <- tf.targets[names(scores)]
              allscores[[p]] <- list(scores = scores, targets = targets)
              df <- data.frame(gene = names(scores), score = as.numeric(scores))
              df <- df[order(df$score, decreasing = TRUE), ]
              angle = 30
              if (flip){
                df$gene <- factor(df$gene, levels = rev(df$gene))
              } else {
                df$gene <- factor(df$gene, levels = df$gene)
                angle = 45
              }
              g <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = score,
                                                    fill = score)) +
                ggplot2::geom_bar(stat = 'identity', color='black') +
                ggplot2::scale_fill_gradientn("",colours=colors) +
                ggplot2::theme_bw() + ggplot2::xlab("") +
                ggplot2::ylab("GRN score") +
                ggplot2::theme(legend.position="none") +
                ggplot2::theme(axis.text.x =
                                 ggplot2::element_text(size = ggplot2::rel(1),
                                                       angle = angle,
                                                       hjust = 1),
                               axis.text.y =
                                 ggplot2::element_text(size = ggplot2::rel(0.9))) +
                ggplot2::ggtitle(paste(p)) +
                ggplot2::theme(panel.background = ggplot2::element_blank(),
                               panel.grid.major = ggplot2::element_blank(),
                               panel.grid.minor = ggplot2::element_blank(),
                               plot.background=ggplot2::element_blank(),
                               panel.border = ggplot2::element_blank(),
                               axis.line.x = ggplot2::element_line(size = 0.5,
                                                                   linetype = "solid",
                                                                   colour = "black"),
                               axis.line.y = ggplot2::element_line(size = 0.5,
                                                                   linetype = "solid",
                                                                   colour = "black"))
              if(flip){
                g <- g + ggplot2::coord_flip()
              }
              plots[[p]] <- g
            }
            return(list(grn.scores = allscores, plots = plots))
          }
)
