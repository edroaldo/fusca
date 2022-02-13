#' Process curves
#'
#' Process curves found in findCurves.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param genes character vector; gene names that are used for trajectory
#' analysis.
#' @param num.cells numeric; minimum of cells contained in the trajectories.
#' @param neighs numeric; the size of the neighborhood in kNN graph used to
#' smoothen kinetic profiles.
#' @param column.ann character; transitions between the cell states provided
#' here are identified, such as clusters or sorted populations.
#' @param column.color character; the colors corresponding to the annotation
#' in column.ann.
#'
#' @return CellRouter object with the curves and the pathsinfo slots updated.
#'
#' @export
#' @docType methods
#' @rdname processCurves-methods
setGeneric("processCurves", function(object, assay.type='RNA',
                                     genes, num.cells,
                                     neighs, column.ann='population',
                                     column.color='colors')
  standardGeneric("processCurves"))
#' @rdname processCurves-methods
#' @aliases processCurves
setMethod("processCurves",
          signature = "CellRouter",
          definition = function(object, assay.type, genes, num.cells, neighs,
                                column.ann, column.color){
            #browser()
            object@genes.trajectory <- genes
            paths <- object@curves
            paths <- as.data.frame(do.call(rbind, paths))
            # Remove empty paths.
            # paths <- paths[complete.cases(paths), ]
            # Create paths ID and populations.
            paths$ID <- paste('path', 1:nrow(paths), sep='_')
            paths$population <- rownames(paths)
            object@curves <- paths
            # Call pathinfo.
            object@pathsinfo <- curves_info(object, assay.type,
                                            num.cells = num.cells,
                                            neighs = neighs,
                                            column.ann = column.ann,
                                            column.color = column.color)
            return (object)
          }

)



#' Helper function of processCurves.
#'
#' Get information about the paths between populations. This is the version of
#' pathsinfo for the curves found in findCurves.
#'
#' @param object CellRouter object
#' @param assay.type character; the type of data to use.
#' @param num.cells numeric; number of cells.
#' @param neighs numeric; neighs.
#' @param column.ann character; column corresponding to the population.
#' @param column.color character; column corresponding to the color.
#'
#' @return list; list with distr, path, pathInfo, pseudotime for each path, and
#' path_data.
setGeneric("curves_info", function(object, assay.type='RNA', num.cells, neighs,
                                   column.ann='community', column.color='colors')
  standardGeneric("curves_info"))
setMethod("curves_info",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                num.cells, neighs, column.ann,
                                column.color){
            #browser()
            paths <- object@curves
            keep <- intersect(rownames(slot(object, 'assays')[[assay.type]]@ndata),
                              object@genes.trajectory)
            expDat <- slot(object, 'assays')[[assay.type]]@ndata[keep, ]
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            means <- list()
            # For each transition, calculate the mean of gene expression of each
            # cell in the transition.
            for(transition in rownames(paths)){
              # Get cells in the path.
              cells <- paths[transition, 'cell'][[1]]
              mean <- list()
              # Get expression.
              mean.df <- as.matrix(expDat[ , cells])
              means[[transition]] <- mean.df
            }
            i <- 1
            path_info <- list()
            pathsDF <- data.frame()
            # For each transition, calculates the information about the path.
            for(transition in rownames(paths)){
              cells <- paths[transition, 'cell'][[1]]
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
                             cost = -99,
                             source_population =
                               sampTab[grep(paste('^', cells[1], '$', sep = ''),
                                            rownames(sampTab)), column.ann],
                             target_population =
                               sampTab[grep(paste('^', cells[length(cells)],
                                                  '$', sep=''),
                                            rownames(sampTab)), column.ann])
                pathsDF[path_name, 'source'] <- cells[1]
                pathsDF[path_name, 'sink'] <- cells[length(cells)]
                pathsDF[path_name, 'cost'] <- -99
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
                # pseudotime <- (0:(length(cells)-1)) / length(cells)
                pseudotime <- paths[transition, 'dist'][[1]]
                names(pseudotime) <- cells
                path_info[['pseudotime']][[path_name]] <- pseudotime
                path_info$path_data <- pathsDF
                i <- i + 1
              }
            }
            return (path_info)
          }
)
