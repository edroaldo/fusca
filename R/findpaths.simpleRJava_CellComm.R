#' Identify signaling pathways connecting cell surface receptors to transcriptional regulators in the protein-protein interaction network.
#'
#' Identify trajectories surface receptors to transcriptional regulators in the
#' protein-protein interaction network. The code uses the rJava library and does
#' not write and run a .sh file as the findpaths.simple function does.
#' The java code saves the data in files which address is stored in
#' assay@@sampTab. These files will be processed by the processTrajectories
#' function. This function uses the java code stored in the java folder in the
#' package directory. If you are a linux user and the rJava
#' library is not installed, please check the topic at
#' \url{https://stackoverflow.com/questions/12872699/error-unable-to-load-installed-packages-just-now}.
#'
#'
#' @param object CellComm object.
#' @param sources.targets list; the sources and targets of the interaction
#' receptor as calculated by the sourcestargets function:
#' \code{sts[['receptor']][[interaction]]}.
#' @param file character; the .txt file corresponding to the cell population. ?
#' @param dir character; directory for the selected interaction.
#' @param maindir character; directory where the directories for each
#' interaction will be saved.
#'
#' @return Message indicating conclusion.
#'
#' @export
#' @docType methods
#' @rdname findpaths.simpleRJava-methods
setGeneric("findpaths.simpleRJava", function(object, sources.targets, file, dir,
                                        maindir)
  standardGeneric("findpaths.simpleRJava"))
#' @rdname findpaths.simpleRJava-methods
#' @aliases findpaths.simpleRJava
setMethod("findpaths.simpleRJava",
          signature="CellComm",
          definition = function(object, sources.targets, file, dir, maindir){
            curdir <- getwd()
            print(file)
            print(maindir)
            dir.create(file.path(maindir, dir), showWarnings = FALSE)
            setwd(file.path(maindir, dir))
            #dirs[[dir]] <- file.path(maindir, dir)

            sources <- sources.targets[['sources']]
            targets <- sources.targets[['targets']]
            write.table(sources, file='sources.txt', sep=',', row.names=FALSE,
                        quote=FALSE)
            write.table(targets, file='sinks.txt', sep=',', row.names=FALSE,
                        quote=FALSE)
            print(paste('NET=', maindir, file, sep='/'))

            # Use Java Code.
            jCC <- rJava::.jnew("flow.Flow")
            # Create java variables.
            OUT_DIR <- rJava::.jnew("java.lang.String", maindir)
            #
            # Duvida: Esse NET ta certo? Nao precisaria ter o nome de algum arquivo?
            #
            NET <- rJava::.jnew("java.lang.String", file.path(maindir, file))
            network <- rJava::.jnew("java.lang.String", file.path(maindir, dir,
                                                                  'FlowNetwork'))
            source <- rJava::.jnew("java.lang.String", file.path(maindir, dir,
                                                                 'sources.txt'))
            sink <- rJava::.jnew("java.lang.String", file.path(maindir, dir,
                                                               'sinks.txt'))
            #
            # duvida: deixa 25 mesmo ou coloca como parâmetro?
            #
            topPaths <- as.integer(25)

            # How the java code was called.
            # cat('java -cp $LIB "flow.Flow" -NET $NET -source $source
            # -sink $sink -topPaths $topPaths -out $OUT_DIR -C -f $network\n')

            # How the arguments were read in the Java method.
            # String filename = cmd.getOptionValue("f");
            # String dir = cmd.getOptionValue("out");
            # int top = Integer.valueOf(cmd.getOptionValue("topPaths"));
            #
            # String cellSources = cmd.getOptionValue("source");
            # String cellSinks = cmd.getOptionValue("sink");

            # Arguments from runAnalysis method.
            # file = cmd.getOptionValue("NET");
            file <- NET

            # How the method is called now.
            # createFlowNetwork(String filename, String dir, int top,
            # String cellSources, String cellSinks, String file)
            filename <- network
            dir <- OUT_DIR
            top <- topPaths
            cellSources <- source
            cellSinks <- sink
            rJava::.jcall(jCC, returnSig = "V", method = "createFlowNetwork",
                          filename, dir, top, cellSources, cellSinks,
                          file)

            # Create gml file.
            subnet <- read.table('FlowNetwork_all_paths_subnet.txt', sep='	',
                                 header = TRUE);
            totalFlow <- read.table('FlowNetwork_all_paths_totalFlow.txt',
                                    sep='	', header = TRUE);
            if (nrow(subnet) > 0){
              colnames(subnet) <- c('from', 'to', 'weight', 'correlation',
                                    'sign');
              iG <- igraph::graph.data.frame(subnet, directed=FALSE);
              igraph::V(iG)$label <- igraph::V(iG)$name;
              igraph::V(iG)[as.vector(totalFlow$cell)]$totalFlows <- totalFlow$flows;
              igraph::V(iG)$degree <- igraph::degree(iG);
              neighS <- igraph::V(
                igraph::induced.subgraph(graph=iG,
                                         vids=unlist(
                                           igraph::neighborhood(graph=iG,
                                                                order=1,
                                                                nodes='s'))))$name;
              neighT <- igraph::V(
                igraph::induced.subgraph(graph=iG,
                                         vids=unlist(
                                           igraph::neighborhood(graph=iG,
                                                                order=1,
                                                                nodes='t'))))$name;
              igraph::V(iG)$props = rep("INTER", length(igraph::V(iG)$name));
              igraph::V(iG)[match(neighS, igraph::V(iG)$name)]$props = "SOURCE";
              igraph::V(iG)[match(neighT, igraph::V(iG)$name)]$props = "TARGET";
              iG2 <- igraph::delete.vertices(iG, c('s','t'));
              igraph::write.graph(iG2,
                                  'FlowNetwork_all_paths_subnet.gml',
                                  format = 'gml');
            } else {
              cat('No cells in the path', dir,'. The gaph was not created.\n')
            }

            # Change directory back to the original one.
            setwd(file.path(curdir))
            #object@directory <- dirs
            return(object)
          }
)
