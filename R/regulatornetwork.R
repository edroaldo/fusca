#' Plot subnetworks for selected transcriptional regulators in the trajectories.
#'
#' In development, currently not working.
#' Plot a network module enriched around transcriptional regulators of
#' cell-fate transitions.
#'
#' @param object CellRouter object.
#' @param x list; GRN scores as calculated by the grnscores function.
#' @param ggrn igraph; the gene regulatory network.
#' @param genelist character vector; gene names containing the transcriptional
#' regulators from which the network modules will be identified.
#'
#' @return list; rgrn graphs for all names in x, and plot.
#'
#' @export
regulatornetwork <- function(object, x, ggrn, genelist){
  #browser()
  nets <- list()
  nets2 <- list()
  for(t in names(x)){
    for(n in genelist){
      genes <- x[[t]]$targets[[n]]
      if(length(genes) > 0){
        genes <- c(n, genes)
        rgrn <- igraph::induced.subgraph(ggrn,
                                         vids = unlist(
                                           igraph::neighborhood(graph = ggrn,
                                                                order = 0,
                                                                nodes = genes)))
        remove <- igraph::V(rgrn)$name[igraph::degree(rgrn) == 0]
        rgrn <- igraph::delete.vertices(rgrn, remove)
        igraph::V(rgrn)$color <- object@correlation[[t]][igraph::V(rgrn)$name]
        igraph::V(rgrn)$size <- log(igraph::degree(rgrn) + 1)
        g <- ggplot2::fortify(rgrn)
        g$network <- n
        g$transition <- t
        name <- paste(t, n, sep='_')
        nets2[[name]] <- rgrn
        nets[[name]] <- g
      }
    }
  }
  # duvida: Converter para dataframe passou o erro, mas agora da um novo erro
  # Faceting variables must have at least one value
  l <- as.data.frame(do.call(rbind, nets))
  l$network <- factor(l$network, levels = genelist)
  l$label[which(l$type == 'Target')] <- ""
  # Choose a threshold for each network independently.
  q <- quantile(l$color, 0.25)
  #browser()
  # duvida: não existe l$type, tem algo antes dando errado ou foi apagado?
  xxx <- l[which(l$type == 'Regulator'), ]
  xxx <- as.vector(xxx[which(xxx$color < q), 'label'])
  xxx <- xxx[!(xxx %in% genelist)]
  l$label[l$label %in% xxx] <- ""
  set.seed(1)
  g <- ggplot2::ggplot(data = l, ggplot2::aes(from_id = from, to_id = to)) +
    geomnet::geom_net(ggplot2::aes(colour = color, size = size, label = label),
                      layout.alg = "kamadakawai", labelon = TRUE, vjust = -0.6,
                      ecolour = "grey60", directed = FALSE, fontsize = 2,
                      ealpha = 0.1, labelcolour = 'black', fiteach=T,
                      arrowsize = 1) +
    ggplot2::xlim(c(-0.1, 1.1)) + ggplot2::ylim(-0.1, 1.1) +
    geomnet::theme_net() + ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_colour_gradientn(colours = c('cyan', 'white', 'red')) +
    ggplot2::facet_grid(network~transition, scales='free') +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA,
                                                        colour = "black"),
                   strip.background = ggplot2::element_rect(colour="black"))
  return(list(nets = nets2, plot = g))
}
