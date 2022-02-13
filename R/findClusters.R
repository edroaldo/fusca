#' Identify clusters.
#'
#' Identify clusters based on graph-clustering or model based clustering. For
#' model based clustering it is necessary to load the mclust package manually.
#'
#' @param object CellRouter object
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param method character; graph.clustering or model.clustering.
#' @param k numeric; number of nearest neighbors to build a k-nearest neighbors
#' graph.
#' @param num.pcs numeric; number of principal components that will define the
#' space from where the kNN graph is identified. For example, if num.pcs = 10,
#' the kNN graph will be created from a 10-dimensional PCA space.
#' @param nn.type character; method to find the k-nearest neighbor graph. If
#' 'nng', the code will use the nng function of the cccd package. If 'knn' or
#' 'snn', it will use the functions from the bluster package, which are
#' recommended for big data due to their efficiency.
#' @param sim.type character; updates the kNN graph to encode cell-cell
#' similarities using the jaccard or invlog methods.
#'
#' @return CellRouter object with the slots updated according to the chosen
#' method: graph.clustering updates graph, rawdata, ndata, and sampTab;
#' model.clustering updates sampTab.
#'
#' @export
#' @docType methods
#' @rdname findClusters-methods
setGeneric("findClusters", function(object, assay.type='RNA',
                                    sample.name='Sample1',
                                    method = 'graph.clustering',
                                    k = 20, num.pcs = 20,
                                    nn.type = "nng", sim.type = 'jaccard')
  standardGeneric("findClusters"))
#' @rdname findClusters-methods
#' @aliases findClusters
setMethod("findClusters",
          signature = "CellRouter",
          definition = function(object, assay.type, sample.name,
                                method, k, num.pcs, nn.type, sim.type){
            if(method == 'graph.clustering'){
              cat('Graph-based clustering\n')
              cat('k: ', k, '\n')
              cat('similarity type: ', sim.type, '\n')
              cat('number of principal components: ', num.pcs, '\n')
              object <- graphClustering(object, assay.type,
                                        k = k, num.pcs = num.pcs,
                                        nn.type, sim.type)
            }else if(method=='model.clustering'){
              cat('Model-based clustering\n')
              cat('number of principal components: ', num.pcs)
              object <- modelClustering(object, assay.type, num.pcs = num.pcs)
            }
            if (assay.type=='ST'){
              bcs <- slot(object, 'assays')[['ST']]@image$bcs[[sample.name]]
              clusters <- slot(object, 'assays')[['ST']]@
                sampTab[, c('sample_id','population','colors')]
              bcs <- merge(bcs, clusters,
                           by.x = "barcode", by.y = "sample_id", all = TRUE)
              bcs <- bcs[!is.na(bcs$tissue), ]
              slot(object, 'assays')[['ST']]@image$bcs[[sample.name]] <- bcs
            }
            object
          }
)



#' Graph-based clustering.
#'
#' Identify clusters using graph-based model. Use jaccard or invlog similarity.
#' Update graph, rawdata, ndata and sampTab.
#'
#' @param object CellRouter object
#' @param assay.type character; the type of data to use.
#' @param k numeric; number of nearest neighbors to build a k-nearest neighbors
#' graph.
#' @param num.pcs numeric; number of principal components that will define the
#' space from where the kNN graph is identified. For example, if num.pcs = 10,
#' the kNN graph will be created from a 10-dimensional PCA space.
#' @param nn.type character; method to find the k-nearest neighbor graph. If
#' 'nng', the code will use the nng function of the cccd package. If 'knn' or
#' 'snn', it will use the functions from the bluster package, which are
#' recommended for big data due to their efficiency.
#' @param sim.type character; updates the kNN graph to encode cell-cell
#' similarities using the jaccard or invlog methods.
#' @param filename character; save .gml file containing the kNN graph.
#'
#' @return CellRouter object with the graph, rawdata, ndata, and sampTab slots
#' updated.
#'
#' @export
#' @docType methods
#' @rdname graphClustering-methods
setGeneric("graphClustering", function(object, assay.type='RNA',
                                       k = 5, num.pcs,
                                       nn.type = "nng",
                                       sim.type = "jaccard",
                                       filename = "graph_subpopulations.gml")
  standardGeneric("graphClustering"))
#' @rdname graphClustering-methods
#' @aliases graphClustering
setMethod("graphClustering",
          signature = "CellRouter",
          definition = function(object, assay.type, k, num.pcs, nn.type, sim.type, filename){
            # Build a k-nn graph and create a matrix of similarities for the
            # k-nn cells.
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            matrix <- object@pca$cell.embeddings[, 1:num.pcs]
            print('building k-nearest neighbors graph')
            if (nn.type == 'knn'){
              h <- bluster::makeKNNGraph(x = matrix, k = k, directed = TRUE,
                                         BPPARAM = BiocParallel::MulticoreParam())
            } else if(nn.type == 'snn') {
              h <- bluster::makeSNNGraph(x = matrix, k = k, type="jaccard",
                                         BPPARAM = BiocParallel::MulticoreParam())
            } else {
              dm <- as.matrix(proxy::dist(matrix))
              h <- cccd::nng(dx = dm, k = k)
            }
            if (sim.type == 'jaccard') {
              # Convert to sparse matrix.
              sim <- as(object = igraph::similarity.jaccard(h,
                                                            vids = igraph::V(h),
                                                            loops = FALSE),
                        Class = 'dgCMatrix')
              # sim <- igraph::similarity.jaccard(h, vids=igraph::V(h), loops=FALSE)
            } else if (sim.type == 'invlog') {
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
            # Community detection to discover subpopulation structure.
            print('discoverying subpopulation structure')
            comms <- igraph::multilevel.community(igraph::as.undirected(h),
                                                  weights = igraph::E(h)$weight)
            igraph::V(h)$comms <- igraph::membership(comms)
            cell.comms <- commToNames(comms, '') #SP means SubPopulations
            allcells <- as.vector(unlist(cell.comms))
            # Making sure that color mappings are correct.
            sampTab <- sampTab[allcells,] #changesorder of cells in the table
            sampTab$population <- ''
            sampTab$colors <- ''
            comm.colors <- cRampClust(unique(igraph::membership(comms)), 8)
            names(comm.colors) <- names(cell.comms)
            for(comm in names(cell.comms)){
              sampTab[cell.comms[[comm]], 'population'] <- comm
              sampTab[cell.comms[[comm]], 'colors'] <- comm.colors[comm]
            }
            sampTab$community <- as.vector(sampTab$population)
            # Mapping information to the igraph object.
            igraph::V(h)[rownames(sampTab)]$subpopulation <- sampTab$colors
            igraph::V(h)[rownames(sampTab)]$colors <- sampTab$colors
            igraph::V(h)[names(nodeLabels(sampTab,'community'))]$label <-
              nodeLabels(sampTab, 'community')
            igraph::V(h)$size <- 5
            igraph::E(h)$arrow.size <- 0.01
            colors <- rainbow(max(igraph::membership(comms)))
            # Useful information about the graph.
            graph <- list()
            graph[['network']] <- h
            graph[['edges']] <- edges
            graph[['similarity_matrix']] <- sim
            graph[['subpopulation']] <- cell.comms
            graph[['communities']] <- comms
            print('updating CellRouter object')
            object@graph <- graph
            # Update rawdata, if it is stored.
            if(nrow(slot(object, 'assays')[[assay.type]]@rawdata) > 0){
              slot(object, 'assays')[[assay.type]]@rawdata <-
                slot(object, 'assays')[[assay.type]]@rawdata[, rownames(sampTab)]
            }
            # Update ndata and sampTab.
            slot(object, 'assays')[[assay.type]]@ndata <-
              slot(object, 'assays')[[assay.type]]@ndata[, rownames(sampTab)]
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




#' Assign names to subpopulations.
#'
#' Assign names to subpopulations. Found in graphClustering. Helper function of
#' findClusters.
#'
#' @param commObj igraph::multilevel.community; subpopulations.
#' @param prefix character; prefix.
commToNames <- function(commObj, prefix){
  ans <- list();
  comms <- igraph::communities(commObj);
  for(i in seq(length(comms))){
    nname<-paste(prefix, "", i, sep='');
    ans[[nname]] <- comms[[i]];
  }
  ans;
}

#' Assign labels for subpopulations.
#'
#' Assign labels for (for plotting). Found in graphClustering. Helper function
#' of findClusters.
#'
#' @param sampTab data frame; metadata information.
#' @param column character; column corresponding to subpopulation.
nodeLabels <- function(sampTab, column){
  x <- sampTab
  x$label <- ''
  xx <- vector()
  names <- vector()
  for(c in unique(as.vector(x[[column]]))){
    tmp <- x[which(x[[column]] == c), ]
    label <- c(unique(tmp[[column]]), rep('', nrow(tmp) - 1))
    tmp$label <- sample(label)
    xx <- append(xx, as.vector(tmp$label))
    names <- append(names, rownames(tmp))
  }
  names(xx) <- names
  xx
}



#' Model-based clustering using the Mclust package.
#'
#' Perform model-based clustering and updates sampTab.
#'
#' @param object CellRouter object
#' @param assay.type character; the type of data to use.
#' @param num.pcs numeric; number of principal components that will define the
#' space used as input to perform model-based clustering.
#'
#' @import mclust
setGeneric("modelClustering", function(object, assay.type = 'RNA', num.pcs)
  standardGeneric("modelClustering"))
setMethod("modelClustering", signature = "CellRouter",
          definition = function(object, assay.type, num.pcs){
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            colname <- 'population'
            matrix <- object@pca$cell.embeddings[,1:num.pcs]
            mclust <- mclust::Mclust(matrix)
            sampTab[names(mclust$classification),colname] <-
              as.character(mclust$classification)
            sampTab <- sampTab[order(as.numeric(sampTab[[colname]])),]
            colors <- cRampClust(1:length(unique(sampTab[[colname]])), 8)
            names(colors) <- unique(sampTab[[colname]])
            replicate_row <-
              as.vector(unlist(lapply(split(sampTab, sampTab[[colname]]), nrow)))
            colors_row <- rep(colors, times=replicate_row)
            sampTab[, 'colors'] <- colors_row
            slot(object, 'assays')[[assay.type]]@sampTab <- sampTab
            object
          }
)
