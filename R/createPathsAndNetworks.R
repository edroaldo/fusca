#' Create paths.
#'
#' Create paths and return network.
#'
#' @param summary summary
#' @param paths paths
#' @param data data
#' @param sampTab metadata table
#' @param top.paths top paths
#'
#' @return network
#'
#' @export
createPaths <- function(summary, paths, data, sampTab, top.paths){
  networks <- list()
  for(n in names(paths)){
    allpaths <- paths[[n]]
    p <- 1
    gpath <- list()
    for(path in as.vector(allpaths$path)){
      df <- data.frame(matrix(0, ncol=3))
      x <- sampTab[which(sampTab$population == n),]
      xx <- data[, rownames(x)]
      colnames(df) <- c('from','to', 'weight')
      i <- 1
      split_paths <- strsplit(as.vector(path), "->")[[1]]
      gene <- split_paths[c(-1, -length(split_paths))]
      if(length(gene) != 1){
        for(g in 1:(length(gene)-1)){
          gi <- gene[g]
          gj <- gene[g+1]
          df[i,'from'] <- gi
          df[i,'to'] <- gj
          df[i, 'weight'] <- 1
          i <- i + 1
        }
      }
      g <- igraph::graph.data.frame(df)
      g <- igraph::simplify(g)
      igraph::V(g)$label <- igraph::V(g)$name
      gg <- summary$gml[[n]]
      gs <- intersect(igraph::V(gg)$name, igraph::V(g)$name)
      igraph::V(g)[gs]$props <- igraph::V(gg)[gs]$props
      colrs <- c("cyan", "tomato", "gold")
      names(colrs) <- c('TARGET', 'SOURCE', 'INTER')
      igraph::V(g)$color <- colrs[igraph::V(g)$props]
      igraph::V(g)$size <- 10
      igraph::E(g)$arrow.size=0.1
      name <- paste('path', p, sep='_')
      gpath[[name]] <- g
      p <- p + 1
    }
    networks[[n]] <- gpath
  }
  return(networks)
}

#' Create networks from paths.
#'
#' Create and return networks.
#'
#' @param summary summary
#' @param paths paths
#' @param n.targets number of targets
#' @param data data
#' @param sampTab metadata
#' @param top.paths top paths
#'
#' @return networks
#'
#' @export
createNetworksFromPaths <- function(summary, paths, n.targets, data, sampTab, top.paths){
  networks <- list()
  for(n in names(paths)){
    cat("---------- ", n, ' ----------\n')
    df <- data.frame(matrix(0, ncol=3))
    x <- sampTab[which(sampTab$population == n),]
    xx <- data[ , rownames(x)]
    colnames(df) <- c('from','to', 'corr')
    i <- 1
    allpaths <- paths[[n]]
    cat('Number of paths: ', nrow(allpaths), '\n')
    cat('---------: ', top.paths)
    targets <- n.targets[[n]]
    if(top.paths > nrow(allpaths)){
      top.paths <- nrow(allpaths)
      cat(nrow(allpaths), '----\n')
      cat(top.paths, '-----\n')
    }
    for(path in as.vector(allpaths$path)){
      print(path)
      split_paths <- strsplit(as.vector(path), "->")[[1]]
      gene <- split_paths[c(-1,-length(split_paths))]
      if(length(gene) > 1){
        for(g in 1:(length(gene)-1)){
          print('here')
          gi <- gene[g]
          gj <- gene[g+1]
          df[i,'from'] <- gi
          df[i,'to'] <- gj
          i <- i + 1
        }
      }
    }
    print(df)
    g <- igraph::graph.data.frame(df)
    g <- igraph::simplify(g)
    igraph::V(g)$label <- igraph::V(g)$name
    gg <- summary$gml[[n]]
    gs <- intersect(igraph::V(gg)$name, igraph::V(g)$name)
    igraph::V(g)[gs]$props <- igraph::V(gg)[gs]$props
    igraph::V(g)$label[which(igraph::V(g)$props == 'INTER')] <- ""
    colrs <- c("cyan", "tomato", "gold")
    names(colrs) <- c('TARGET', 'SOURCE', 'INTER')
    igraph::V(g)$color <- colrs[igraph::V(g)$props]
    igraph::V(g)$size <- 10
    xxx <- intersect(names(targets), igraph::V(g)$name)
    targets <- targets[xxx]
    igraph::V(g)[names(targets)]$size <- targets
    igraph::E(g)$arrow.size=0.1
    networks[[n]] <- g
  }
  return(networks)
}
