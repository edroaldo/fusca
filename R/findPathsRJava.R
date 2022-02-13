#' Identify trajectories connecting source and target populations.
#'
#' Identify trajectories based on euclidean, maximum, manhattan, canberra,
#' binary or graph network. The code uses the rJava library and does not
#' communicate with the terminal as the findPaths function does. The
#' computations are stored in the specified directories and can be read using
#' the summarize function. This function uses the java code stored in the java
#' folder in the package directory. If you are a linux user and the rJava
#' library is not installed, please check the topic at
#' \url{https://stackoverflow.com/questions/12872699/error-unable-to-load-installed-packages-just-now}.
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
#' based on the source and target populations. The available methods are
#' 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'graph'.
#'
#' @return CellRouter object with the directory slot updated.
#'
#' @export
#' @docType methods
#' @rdname findPathsRJava-methods
setGeneric("findPathsRJava", function(object, assay.type='RNA',
                                      sources, targets,
                                      column = 'population', maindir, method)
  standardGeneric("findPathsRJava"))
#' @rdname findPathsRJava-methods
#' @aliases findPathsRJava
setMethod("findPathsRJava",
          signature = "CellRouter",
          definition = function(object, assay.type,
                                sources, targets, column = 'population',
                                maindir, method){
            curdir <- getwd()
            dirs <- list()
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
            # Return a list b with the cells with the greater difference between
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
            # Creates the java virtual machine.
            # This step was done in .onLoad.R
            # rJava::.jinit()
            # class_paths <- c(paste(libdir, 'dist/CellRouter.jar', sep = '/'),
            #                  paste(libdir, 'lib/commons-cli-1.2.jar', sep = '/'))
            # class_paths <- c(paste(libdir, 'CellRouter.jar', sep = '/'),
            #                  paste(libdir, 'commons-cli-1.2.jar', sep = '/'))
            # rJava::.jaddClassPath(class_paths)
            # Creates the CellRouter object.
            jCR <- rJava::.jnew("cellrouter.CellRouter")
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
                # Set wd to transition folder.
                setwd(file.path(maindir, dir))
                dirs[[dir]] <- file.path(maindir, dir)
                # Write table with source and target.
                write.table(s, file='cell_sources.txt', sep=',',
                            row.names = FALSE, quote = FALSE)
                write.table(t, file='cell_sinks.txt', sep=',',
                            row.names = FALSE, quote = FALSE)
                # Create java variables.
                OUT_DIR <- rJava::.jnew("java.lang.String", maindir)
                #
                # duvida: pode escrever dentro da pasta?
                #
                NET <- rJava::.jnew("java.lang.String", file.path(maindir, 'cell_edge_weighted_network.txt'))
                network <- rJava::.jnew("java.lang.String", file.path(maindir, dir, 'Cells_FlowNetwork'))
                source <- rJava::.jnew("java.lang.String", file.path(maindir, dir, 'cell_sources.txt'))
                sink <- rJava::.jnew("java.lang.String", file.path(maindir, dir, 'cell_sinks.txt'))
                #
                # duvida: deixa 25 mesmo ou coloca como parâmetro?
                #
                topPaths <- as.integer(25)
                out <- rJava::.jnew("java.lang.String", maindir)
                f <- network
                # Call the createFlowNetwork method, attributes should be passed
                # in the correct order.
                # public void createFlowNetwork(String filename, String dir,
                # int top, String cellSources, String cellSinks, String file)
                # createFlowNetwork(f, out, topPaths, source, sink, file);
                rJava::.jcall(jCR, returnSig = "V", method = "createFlowNetwork",
                              f, out, topPaths, source, sink, NET)

                # Create gml file.
                subnet <- read.table('Cells_FlowNetwork_all_paths_subnet.txt',
                                     sep='	', header = TRUE);
                totalFlow <- read.table('Cells_FlowNetwork_all_paths_totalFlow.txt',
                                        sep='	', header = TRUE);
                if (nrow(subnet) > 0){
                  colnames(subnet) <- c('from', 'to', 'weight');
                  iG <- igraph::graph.data.frame(subnet, directed = FALSE);
                  igraph::V(iG)$label <- igraph::V(iG)$name;
                  igraph::V(iG)[as.vector(totalFlow$cell)]$totalFlows <- totalFlow$flows;
                  igraph::V(iG)$degree <- igraph::degree(iG);
                  neighS <-
                    igraph::V(igraph::induced.subgraph(graph = iG,
                                                       vids = unlist(igraph::neighborhood(
                                                         graph = iG, order = 1,
                                                         nodes = 's'))))$name;
                  neighT <-
                    igraph::V(igraph::induced.subgraph(graph = iG,
                                                       vids = unlist(igraph::neighborhood(
                                                         graph = iG, order = 1,
                                                         nodes = 't'))))$name;
                  igraph::V(iG)$props = rep("INTER", length(igraph::V(iG)$name));
                  igraph::V(iG)[match(neighS, igraph::V(iG)$name)]$props = "SOURCE";
                  igraph::V(iG)[match(neighT, igraph::V(iG)$name)]$props = "TARGET";
                  iG2 <- igraph::delete.vertices(iG, c('s','t'));
                  igraph::write.graph(iG2,
                                      'Cells_FlowNetwork_all_paths_subnet.gml',
                                      format = 'gml');
                } else {
                  cat('No cells in the path', dir,'. The gaph was not created.\n')
                }
                # Change directory back to the original one.
                setwd(file.path(curdir))
              }
            }
            object@directory <- dirs
            return(object)
          }
)
