#' Build kNN graph.
#'
#' Create kNN graph using the PCA components and cell-cell similarities using
#' the jaccard or invlog methods. Write graph to file.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param k numeric; number of nearest neighbors to build a k-nearest neighbors
#' graph for trajectory reconstruction.
#' @param column.ann character; column in the metadata table specifying the
#' transitions to be identified. For example, if 'population' is provided,
#' transitions will be identified between clusters previously identified.
#' However, sorted cell populations or customized states can also be used.
#' Check our tutorials for detailed examples.
#' @param num.pcs numeric; number of principal components that will define the
#' space from where the kNN graph is identified. For example, if num.pcs = 10,
#' the kNN graph will be created from a 10-dimensional PCA space.
#' @param sim.type character; updates the kNN graph to encode cell-cell
#' similarities using the jaccard or invlog methods.
#' @param filename character; name of gml file containing the kNN graph.
#'
#' @return CellRouter object with the graph, ndata, and sampTab slots updated.
#'
#' @export
#' @docType methods
#' @rdname buildKNN-methods
setGeneric("buildKNN", function(object, assay.type='RNA',
                                k = 5, column.ann, num.pcs = 20,
                                sim.type = "jaccard",
                                filename = "graph_clusters.gml")
  standardGeneric("buildKNN"))
#' @rdname buildKNN-methods
#' @aliases buildKNN
setMethod("buildKNN",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                k, column.ann,  num.pcs, sim.type,
                                filename){
            # Build a k-nn graph and create a matrix of similarities for the
            # k-nn cells.
            matrix <- object@pca$cell.embeddings[, 1:num.pcs]
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            smapTab <- sampTab[order(sampTab[[column.ann]]),]
            print('building k-nearest neighbors graph')
            dm <- as.matrix(proxy::dist(matrix))
            h <- cccd::nng(dx = dm, k = k)
            if(sim.type == 'jaccard'){
              sim <- as(object = igraph::similarity.jaccard(h,
                                                            vids = igraph::V(h),
                                                            loops = FALSE),
                        Class = 'dgCMatrix')
              # sim <- igraph::similarity.jaccard(h, vids=igraph::V(h), loops=FALSE)
            }else if(sim.type == 'invlog'){
              sim <- as(object = igraph::similarity.invlogweighted(h,
                                                                   vids = igraph::V(h)),
                        Class = 'dgCMatrix')
              # sim <- igraph::similarity.invlogweighted(h, vids=igraph::V(h))
            }
            el <- igraph::get.edgelist(h)
            weights <- sim[el]
            igraph::E(h)$weight  <- weights
            igraph::V(h)$name <- rownames(matrix)
            edges <- as.data.frame(igraph::get.edgelist(h))
            rownames(edges) <- paste(edges$V1, edges$V2, sep='_')
            edges$weight <- as.numeric(igraph::E(h)$weight)
            rownames(sim) <- igraph::V(h)$name
            colnames(sim) <- igraph::V(h)$name
            # Name the vertexes after the populations.
            igraph::V(h)[rownames(sampTab)]$comms <- as.vector(sampTab[[column.ann]])
            cell.comms <- split(sampTab, sampTab[[column.ann]])
            cell.comms <- lapply(cell.comms, rownames)
            # Useful information about the graph.
            graph <- list()
            graph[['network']] <- h
            graph[['edges']] <- edges
            graph[['similarity_matrix']] <- sim
            graph[['subpopulation']] <- cell.comms
            print('updating CellRouter object')
            object@graph <- graph
            slot(object, 'assays')[[assay.type]]@ndata <-
              slot(object, 'assays')[[assay.type]]@ndata[,rownames(sampTab)]
            slot(object, 'assays')[[assay.type]]@sampTab <- sampTab
            # Supressed the warning
            # At foreign.c:2616 :A boolean graph attribute was converted to numeric
            suppressWarnings(igraph::write.graph(graph = h, file = filename,
                                                 format = c('gml')))
            rm(h)
            rm(edges)
            rm(sim)
            return(object)
          }
)
