#' Summarize flow network.
#'
#' Resume the results of the findpaths.simple function.
#'
#' @param files list; the filenames passed to the findpaths.simple function.
#' @param prefix character; the path to the files.
#'
#' @return list; the nodes, edges, paths and graphs (gml) for each celltype.
#'
#' @export
summarize.flow <- function(files, prefix){
  gs <-list()
  nodes <- list()
  edges <- list()
  paths <- list()
  for(sub in names(files)){
    sub <- gsub(('[ ]'), "_", sub)
    cat(sub, '\n')

    folder <- file.path(getwd(), prefix, sub, 'FlowNetwork_all_paths_subnet.gml')
    x <- file.path(getwd(), prefix, sub, 'FlowNetwork_all_paths.txt')

    path <- read.csv(x, sep='\t')
    path <- path[!duplicated(as.vector(path$path)),]
    path <- path[order(path$rank),]
    if(nrow(path) > 0){
      gml <- igraph::read.graph(folder, format="gml")
      gs[[sub]] <- gml
      category <- igraph::V(gml)$props
      flow <- igraph::V(gml)$totalFlows
      names(category) <- igraph::V(gml)$name
      names(flow) <- igraph::V(gml)$name

      node <- data.frame(category = category, flow = flow)
      node <- node[order(node$flow, decreasing = TRUE), ]

      # Try to integrate mean total flow to the score...
      df <- data.frame()
      i <- 1
      for(p in as.vector(path$path)){
        genes <- strsplit(p, split='->', fixed=TRUE)[[1]]
        genes <- genes[-c(1, length(genes))]
        df[i,'path'] <- p
        df[i, 'meanflow'] <- mean(node[genes, 'flow'])
        i <- i + 1
      }
      path$meanflow <- as.numeric(df$meanflow)

      edge <- as.data.frame(igraph::get.edgelist(gml))
      edge <- cbind(edge, weight = igraph::E(gml)$weight,
                    corr= igraph::E(gml)$correlation)
      edge <- edge[order(edge$weight, decreasing = TRUE),]

      nodes[[sub]] <- node
      edges[[sub]] <- edge
      paths[[sub]] <- path
    }
  }
  aux <- list(nodes=nodes, edges=edges, paths=paths, gml=gs)
  return(aux)
}
