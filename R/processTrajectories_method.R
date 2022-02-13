#' Process trajectories.
#'
#' Process trajectories found in findPaths.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genes character vector; gene names that are used for trajectory
#' analysis.
#' @param path.rank character; method to rank paths: path_cost, path_flow, rank,
#' length. Available ranks are 'path_cost', 'path_flow', 'rank', 'length'.
#' @param num.cells numeric; minimum of cells contained in the trajectories.
#' @param neighs numeric; the size of the neighborhood in kNN graph used to
#' smoothen kinetic profiles.
#' @param column.ann character; transitions between the cell states provided
#' here are identified, such as clusters or sorted populations.
#' @param column.color character; the colors corresponding to the annotation
#' in column.ann.
#'
#' @return CellRouter object with the paths, the networks and the pathsinfo
#' slots updated.
#'
#' @export
#' @docType methods
#' @rdname processTrajectories-methods
setGeneric("processTrajectories", function(object, assay.type='RNA',
                                           genes, path.rank, num.cells,
                                           neighs, column.ann='population',
                                           column.color='colors')
  standardGeneric("processTrajectories"))
#' @rdname processTrajectories-methods
#' @aliases processTrajectories
setMethod("processTrajectories",
          signature="CellRouter",
          definition=function(object, assay.type,
                              genes, path.rank, num.cells, neighs,
                              column.ann, column.color){
            object@genes.trajectory <- genes
            print('parsing trajectory information')
            # Read the paths from each direcotry created in findPaths.
            paths <- lapply(object@directory,
                            function(x){
                              read.csv(
                                file.path(x, 'Cells_FlowNetwork_all_paths.txt'),
                                sep="\t", stringsAsFactors=FALSE)
                              })
            paths <- paths[lapply(paths, nrow) > 1]
            # Order paths.
            if(path.rank == "path_cost"){
              paths <- lapply(paths, function(x){x[order(x[, path.rank],
                                                         decreasing = FALSE), ]})
            }else{
              paths <- lapply(paths, function(x){x[order(x[, path.rank],
                                                         decreasing = TRUE), ]})
            }
            paths <- lapply(paths, function(x){x[!duplicated(x$path), ]})
            paths <- do.call(rbind, lapply(paths, function(x){x[1,]}))
            # Remove empty paths.
            paths <- paths[complete.cases(paths),]
            # Opening flow networks (graphs).
            networks <- lapply(object@directory,
                               function(x){
                                 file <- list.files(x, '*.gml*');
                                 if(length(file) > 0){
                                   #cat(file, '\n');
                                   igraph::read.graph(
                                     file.path(x, 'Cells_FlowNetwork_all_paths_subnet.gml'),
                                     format = 'gml')
                                 }
                               })
            networks <- networks[lapply(networks, length) > 0]
            # Create paths ID and populations.
            paths$ID <- paste('path', 1:nrow(paths), sep='_')
            paths$population <- rownames(paths)
            object@paths <- paths
            object@networks <- networks
            # Call pathinfo.
            object@pathsinfo <- pathsinfo(object, assay.type,
                                          num.cells = num.cells,
                                          neighs = neighs,
                                          column.ann = column.ann,
                                          column.color = column.color)
            return (object)
          }

)



#' Helper function of processTrajectories.
#'
#' Get information about the paths between populations.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param num.cells numeric; number of cells.
#' @param neighs numeric; neighs.
#' @param column.ann character; column corresponding to the population.
#' @param column.color character; column corresponding to the color.
#'
#' @return list; list with distr, path, pathInfo, pseudotime for each path, and
#' path_data.
setGeneric("pathsinfo", function(object, assay.type='RNA', num.cells, neighs,
                                 column.ann='community', column.color='colors')
  standardGeneric("pathsinfo"))
setMethod("pathsinfo",
          signature="CellRouter",
          definition=function(object, assay.type,
                              num.cells, neighs, column.ann,
                              column.color){
            paths <- object@paths
            keep <- intersect(
              rownames(slot(object, 'assays')[[assay.type]]@ndata),
              object@genes.trajectory)
            expDat <- slot(object, 'assays')[[assay.type]]@ndata[keep,]
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            print(neighs)
            networks <- object@networks
            means <- list()
            o <- neighs
            # For each transition, calculate the mean of gene expression of each
            # cell in the transition.
            for(transition in rownames(paths)){
              g <- networks[[transition]]
              path <- paths[transition, 'path']
              split_paths <- strsplit(as.vector(path), "->")[[1]]
              # Get cells in the path.
              cells <- split_paths[c(-1, -length(split_paths))]
              mean <- list()
              # For each cell in the path, calculates the mean of gene
              # expression in its neigbor cells.
              for(cell in cells){
                neighs <- igraph::induced.subgraph(
                  graph = g, vids = unlist(
                    igraph::neighborhood(graph = g, order = o, nodes = cell)))
                neigh.names <- igraph::V(neighs)$name
                # Get the mean of the expressed genes in the neighbor cells.
                # Does not work for large matrices.
                # mean[[cell]] <- apply(expDat[, neigh.names], 1, mean)
                # It had to be converted to numeric and named.
                mean[[cell]] <- as.numeric(
                  apply1_sp(expDat[, neigh.names], sum)/
                    ncol(expDat[, neigh.names]))
                names(mean[[cell]]) <- rownames(expDat)
              }
              # Do call returns a matrix.
              #
              # duvida: colocar sparse matrix.
              #
              mean.df <- do.call(cbind, mean)
              means[[transition]] <- mean.df
            }
            i <- 1
            path_info <- list()
            pathsDF <- data.frame()
            # For each transition, calculates the information about the path.
            for(transition in rownames(paths)){
              path <- paths[transition, 'path']
              split_paths <- strsplit(as.vector(path), "->")[[1]]
              cells <- split_paths[c(-1,-length(split_paths))]
              # Verifies if the number of cells in the transition is greater
              # than num.cells.
              if(length(cells) > num.cells){
                expression <- means[[transition]]
                path_name <- paths$population[i]
                print(path_name)
                path_info[['distr']][[path_name]] <- expression
                path_info[['path']][[path_name]] <- cells
                path_info[['pathInfo']][[path_name]] <-
                  data.frame(path = paste(as.vector(cells), collapse=','),
                             source = cells[1], sink = cells[length(cells)],
                             cost = paths$path_cost[i],
                             source_population =
                               sampTab[grep(paste('^', cells[1], '$', sep = ''),
                                            rownames(sampTab)), column.ann],
                             target_population =
                               sampTab[grep(paste('^', cells[length(cells)],
                                                  '$', sep=''),
                                            rownames(sampTab)), column.ann])
                pathsDF[path_name, 'source'] <- cells[1]
                pathsDF[path_name, 'sink'] <- cells[length(cells)]
                pathsDF[path_name, 'cost'] <- paths$path_cost[i]
                pathsDF[path_name, 'source_population'] <-
                  as.character(sampTab[grep(paste('^', cells[1], '$', sep = ''),
                                            rownames(sampTab)), column.ann])
                pathsDF[path_name, 'source_color'] <-
                  as.character(sampTab[grep(paste('^', cells[1], '$', sep = ''),
                                            rownames(sampTab)), column.color])
                pathsDF[path_name, 'target_population'] <-
                  as.character(sampTab[grep(paste('^', cells[length(cells)], '$',
                                                  sep = ''),
                                            rownames(sampTab)), column.ann])
                pathsDF[path_name, 'target_color'] <-
                  as.character(sampTab[grep(paste('^', cells[length(cells)],
                                                  '$', sep = ''),
                                            rownames(sampTab)), column.color])
                # The pseudotime is the diference the scalled length of the
                # cell vector.
                pseudotime <- (0:(length(cells)-1)) / length(cells)
                names(pseudotime) <- cells
                path_info[['pseudotime']][[path_name]] <- pseudotime
                path_info$path_data <- pathsDF
                i <- i + 1
              }
            }
            return (path_info)
          }
)
