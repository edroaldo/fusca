#' Identify trajectories connecting populations in the interaction.
#'
#' Identify signaling pathway source and target genes (cell surface receptors
#' and transcriptional regulators, specifically) in the weighted protein
#' interaction network. Create a .sh file and run it using the terminal to call
#' the java code in Flow.jar.
#' The computations are stored in the specified directories and can
#' be read using the summarize function.
#' The java code is written in the CellCommJava folder in the package directory.
#' The java folder in the package directory has the code with some
#' modifications to be used in the findpaths.simpleRJava function.
#'
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
#' @rdname findpaths.simple-methods
setGeneric("findpaths.simple", function(sources.targets, file, dir,
                                        maindir)
  standardGeneric("findpaths.simple"))
#' @rdname findpaths.simple-methods
#' @aliases findpaths.simple
setMethod("findpaths.simple",
          signature = "CellComm",
          definition = function(sources.targets, file, dir, maindir){
            curdir <- getwd()
            # Find path to Java libraries in the package directory.
            libdir <- system.file("CellCommJava", package = "fusca",
                                  mustWork = TRUE)
            print(file)
            print(maindir)

            dir.create(file.path(maindir, dir), showWarnings = FALSE)
            system(paste('chmod 777 ', file.path(maindir, dir), sep=''))
            setwd(file.path(maindir, dir))
            #dirs[[dir]] <- file.path(maindir, dir)

            sources <- sources.targets[['sources']]
            targets <- sources.targets[['targets']]
            write.table(sources, file='sources.txt', sep=',', row.names=FALSE,
                        quote=FALSE)
            write.table(targets, file='sinks.txt', sep=',', row.names=FALSE,
                        quote=FALSE)

            print(paste('NET=', maindir, file, sep='/'))

            sink("CellComm.sh")
            cat('#!/bin/sh\n')
            cat('# run CellComm, run!\n')
            cat(paste('LIB_DIR=', libdir, sep=''), '\n')
            #cat('LIB=$LIB_DIR/dist/Flow.jar:$LIB_DIR/lib/commons-cli-1.2.jar\n')
            cat('LIB=$LIB_DIR/dist/Flow.jar:$LIB_DIR/lib/commons-cli-1.2.jar\n')
            cat(paste('OUT_DIR=', maindir, sep=''), '\n')
            #
            # Duvida: Esse NET ta certo? Nao precisaria ter o nome de algum arquivo?
            #
            #cat(paste('NET=', maindir, '/cell_edge_weighted_network.txt', sep=''), '\n') #-->need to build this network first
            cat(paste('NET=', maindir, file, sep='/'), '\n') #-->need to build this network first
            #cat(paste('NET=', maindir, '/edge_weighted_network.txt', sep=''), '\n') #-->need to build this network first
            cat('network=FlowNetwork\n')
            cat('source=sources.txt\n')
            cat('sink=sinks.txt\n')
            cat('topPaths=25\n')

            cat('echo "1) Computing flow network"\n')
            cat('java -cp $LIB "flow.Flow" -NET $NET -source $source -sink $sink -topPaths $topPaths -out $OUT_DIR -C -f $network\n')
            #cat('java -cp $LIB "cellrouter.CellRouter" -NET $NET -source $source -sink $sink -topPaths $topPaths -out $OUT_DIR -C -f $network\n')
            sink()

            system('chmod 777 CellComm.sh')
            system('./CellComm.sh')

            setwd(file.path(curdir))
            #object@directory <- dirs

            return("Done.")
          }
)
