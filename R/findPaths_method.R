#' Identify trajectories connecting source and target populations.
#'
#' Identify trajectories based on euclidean, maximum, manhattan, canberra,
#' binary or graph network. Writes to the computer terminal calling a java code
#' and that saves the data in files which address is stored in assay@@sampTab.
#' These files will be processed by the processTrajectories function. The java
#' code is written in the CellRouterJava folder in the package directory. The
#' java folder in the package directory has the code with some modifications to
#' be used in the findPathsRJava function.
#'
#' @param object CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sources character vector; names of the source populations.
#' @param targets character vector; names of the target populations.
#' @param column character; column in the metadata table specifying whether
#' transitions are between clusters or other annotations, such as sorted
#' populations.
#' @param maindir character; directory where the files will be saved.
#' @param method character; method used to determine the source and target cells
#' based on the source and target populations.
#'
#' @return CellRouter object with the directory slot updated.
#'
#' @export
#' @docType methods
#' @rdname findPaths-methods
setGeneric("findPaths", function(object, assay.type='RNA', sources, targets,
                                 column = 'population', maindir, method)
  standardGeneric("findPaths"))
#' @rdname findPaths-methods
#' @aliases findPaths
setMethod("findPaths",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                sources, targets, column = 'population',
                                maindir, method){
            curdir <- getwd()
            dirs <- list()
            # Find path to Java libraries in the package directory.
            libdir <- system.file("CellRouterJava", package = "fusca",
                                  mustWork = TRUE)
            # Calculate distance between cells.
            if(method %in% c("euclidean", "maximum", "manhattan", "canberra",
                             "binary")){
              tmp <- as.data.frame(as.matrix(dist(object@rdimension,
                                                  method = method)))
            } else {
              g <- object@graph$network
              tmp <- as.data.frame(igraph::distances(g,
                                                     v = igraph::V(g),
                                                     to = igraph::V(g),
                                                     weights = NULL,
                                                     algorithm = "bellman-ford"))
            }
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            b <- list()
            # Return a list b with the cells with the greates difference between
            # population1 and each other one.
            for (p1 in sources){
              cellsp1 <- as.vector(sampTab[which(sampTab[[column]] == p1),
                                           'sample_id'])
              for (p2 in targets){
                if (p1 != p2){
                  cellsp2 <- as.vector(sampTab[which(sampTab[[column]] == p2),
                                               'sample_id'])
                  # Create a table with cells from population1 as rows, cells
                  # from population2 as columns and their difference as values.
                  x <- tmp[cellsp1, cellsp2]
                  x$population1 <- p1
                  x$population2 <- p2
                  # The names of cells from population1 go to the merge column.
                  x$merge <- rownames(x)
                  # Organize the data with population1, population2 and merge as
                  # indexes and uses the cells from population2 as variables and
                  # the difference as values.
                  tmp2 <- reshape2::melt(x, id.vars = c('population1',
                                                        'population2',
                                                        'merge'))
                  tmp2 <- tmp2[order(tmp2$value, decreasing = TRUE), ]
                  # Get the cells with greatest difference.
                  b[[p1]][[p2]] <- tmp2[1, ]
                }
              }
            }
            # Uses java code.
            for(i in names(b)){
              for(j in names(b[[i]])){
                # s is the cell from population1 (source) and t from population2
                # (target).
                s <- as.vector(b[[i]][[j]][, 'merge'])
                t <- as.vector(b[[i]][[j]][, 'variable'])
                dir <- paste(i, j, sep='.')
                cat('-------------------Transition:', dir, ' -----------------------\n')
                # Create a directory.
                dir.create(file.path(maindir, dir), showWarnings = FALSE)
                # chmod 777 gives permission to read, write and execute files.
                system(paste('chmod 777 ', file.path(maindir, dir), sep=''))
                # Set wd to transition folder.
                setwd(file.path(maindir, dir))
                dirs[[dir]] <- file.path(maindir, dir)
                # Write table with source and target.
                write.table(s, file='cell_sources.txt', sep=',',
                            row.names = FALSE, quote = FALSE)
                write.table(t, file='cell_sinks.txt', sep=',',
                            row.names = FALSE, quote = FALSE)
                # Write bash file to be executed in unix system.
                sink("CellRouter.sh")
                cat('#!/bin/sh\n')
                cat('# run CellRouter, run!\n')
                cat(paste('LIB_DIR=', libdir, sep=''), '\n')
                #cat('LIB=$LIB_DIR/CellRouter.jar:$LIB_DIR/commons-cli-1.2.jar\n')
                cat('LIB=$LIB_DIR/dist/CellRouter.jar:$LIB_DIR/lib/commons-cli-1.2.jar\n')
                cat(paste('OUT_DIR=', maindir, sep=''), '\n')

                cat(paste('NET=', maindir, '/cell_edge_weighted_network.txt', sep=''), '\n')
                cat('network=Cells_FlowNetwork\n')
                cat('source=cell_sources.txt\n')
                cat('sink=cell_sinks.txt\n')
                cat('topPaths=25\n')

                cat('echo "1) Computing flow network"\n')
                cat('java -cp $LIB "cellrouter.CellRouter" -NET $NET -source $source -sink $sink -topPaths $topPaths -out $OUT_DIR -C -f $network\n')
                sink()

                system('chmod 777 CellRouter.sh')
                system('./CellRouter.sh')

                setwd(file.path(curdir))
              }
            }
            object@directory <- dirs

            return(object)
          }
)
