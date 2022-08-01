# Helper functions --------------------------------------------------------

# The geom_spatial function is defined to make plotting your tissue image in
# ggplot a simple task.
# https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/rkit
geom_spatial <-  function(mapping=NULL,
                          data=NULL,
                          stat="identity",
                          position="identity",
                          na.rm=FALSE,
                          show.legend=NA,
                          inherit.aes=FALSE,
                          ...) {

  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data=function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },

    draw_group=function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },

    required_aes=c("grob","x","y")

  )

  layer(
    geom=GeomCustom,
    mapping=mapping,
    data=data,
    stat=stat,
    position=position,
    show.legend=show.legend,
    inherit.aes=inherit.aes,
    params=list(na.rm=na.rm, ...)
  )
}



# Read data ---------------------------------------------------------------

#' Read 10X expression data.
#'
#' Read the spatial transcriptomics expression matrix given a directory with
#' the barcodes (barcodes.tsv.gz), features (features.tsv.gz), and matrix
#' (matrix.mtx.gz) files.
#'
#' @param matrix_dir character; the path to the directory with the barcodes,
#' features, and matrix files. The string should NOT end with your system's path
#' separator.
#'
#' @export
read10X <- function(matrix_dir){
  barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
  features.path <- file.path(matrix_dir, "features.tsv.gz")
  matrix.path <- file.path(matrix_dir, "matrix.mtx.gz")
  matrix <- Matrix::t(Matrix::readMM(file=matrix.path))
  feature.names=read.delim(features.path,
                           header=FALSE,
                           stringsAsFactors=FALSE)
  barcode.names=read.delim(barcode.path,
                           header=FALSE,
                           stringsAsFactors=FALSE)
  rownames(matrix)=gsub('-','.',barcode.names$V1) #rownames(matrix)=barcode.names$V1
  colnames(matrix)=make.unique(feature.names$V2)
  # The matrix had to be transposed to have spots as columns and genes as rows,
  # the default format of cellrouter for expression data.
  return(as(Matrix::t(matrix), "dgCMatrix"))
}



#' Read ST metadata.
#'
#' Read the tissue image, the scale factor and the barcodes of the spatial
#' transcriptomics metadata and add them to the respective assay metadata.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample_names character; names of the tissue samples.
#' @param image_paths character; file path of the tissue image (png format).
#' @param scalefactor_paths character; file path of the scale factor (json
#' format).
#' @param tissue_paths character; file path of the positions of tissue barcodes
#' (csv format).
#'
#' @return the CellRouter object.
#'
#' @export
#' @docType methods
#' @rdname read10XImage-methods
setGeneric("read10XImage", function(object, assay.type='ST', sample_names,
                                       image_paths, scalefactor_paths,
                                       tissue_paths)
  standardGeneric("read10XImage"))
#' @rdname read10XImage-methods
#' @aliases read10XImage
setMethod("read10XImage",
          signature="CellRouter",
          definition=function(object, assay.type, sample_names, image_paths,
                              scalefactor_paths, tissue_paths){
            # Read png image.
            img_bmp <- list()
            for (i in 1:length(sample_names)) {
              img_bmp[[i]] <- readbitmap::read.bitmap(image_paths[i])
            }
            height <- list()
            for (i in 1:length(sample_names)) {
              height[[i]] <-  data.frame(height=nrow(img_bmp[[i]]))
            }
            height <- dplyr::bind_rows(height)
            width <- list()
            for (i in 1:length(sample_names)) {
              width[[i]] <- data.frame(width=ncol(img_bmp[[i]]))
            }
            width <- dplyr::bind_rows(width)
            # Read image tibble.
            # Convert the Images to Grobs to be compatible with ggplot2.
            grobs <- list()
            for (i in 1:length(sample_names)) {
              grobs[[i]] <- grid::rasterGrob(img_bmp[[i]],
                                             width=unit(1,"npc"),
                                             height=unit(1,"npc"))
            }
            # Inserted the tibble to keep the format for multiple assays.
            images_tibble <- dplyr::tibble(sample=factor(sample_names),
                                           grob=grobs)
            images_tibble$height <- height$height
            images_tibble$width <- width$width
            # Read scales.
            scales <- list()
            for (i in 1:length(sample_names)) {
              scales[[i]] <- rjson::fromJSON(file=scalefactor_paths[i])
            }
            # Read barcodes.
            bcs <- list()
            for (i in 1:length(sample_names)) {
              bcs[[i]] <- read.csv(tissue_paths[i],
                                   col.names=c("barcode","tissue","row","col",
                                               "imagerow","imagecol"),
                                   header=FALSE)
              # Scale tissue coordinates for low resolution image.
              bcs[[i]]$imagerow <- bcs[[i]]$imagerow * scales[[i]]$tissue_lowres_scalef
              bcs[[i]]$imagecol <- bcs[[i]]$imagecol * scales[[i]]$tissue_lowres_scalef
              bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue)
              bcs[[i]]$height <- height$height[i]
              bcs[[i]]$width <- width$width[i]
              bcs[[i]]$barcode =  gsub('-','.',bcs[[i]]$barcode) #!
              rownames(bcs[[i]]) <- bcs[[i]]$barcode
            }
            names(bcs) <- sample_names
            # Create metadata and add to cellrouter object.
            metadata <- list(width=width, height=height,
                             images_tibble=images_tibble,
                             scales=scales, bcs=bcs)
            slot(object, 'assays')[[assay.type]]@image <- metadata
            return(object)
          }
)



# Subclusters -------------------------------------------------------------

#' Find subclusters
#'
#' Calculate spatial clusters inside the clusters already identified by
#' gene expression information.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param cluster.column character; the name of the column where the
#' clustering information is stored.
#' @param subcluster.column character; the name of the column where the
#' subclustering information will be stored.
#' @param method character; method to perform clustering (mclust or dbscan).
#' @param eps numeric; eps parameter for the DBSCAN method.
#' @param minPts numeric; minPts parameter for the DBSCAN method.
#'
#' @return the CellRouter object.
#'
#' @import mclust
#'
#' @export
#' @docType methods
#' @rdname findSubclusters-methods
setGeneric("findSubclusters", function(object, assay.type='ST',
                                       sample.name='Sample1',
                                       cluster.column='population',
                                       subcluster.column='Subpopulation',
                                       method=c('mclust', 'dbscan'),
                                       eps=2, minPts=3)
  standardGeneric("findSubclusters"))
#' @rdname findSubclusters-methods
#' @aliases findSubclusters
setMethod("findSubclusters",
          signature="CellRouter",
          definition=function(object, assay.type, sample.name, cluster.column,
                              subcluster.column,
                              method=c('mclust', 'dbscan'),
                              eps=2, minPts=3){
            method <- match.arg(method)
            bcs <- slot(object, 'assays')[[assay.type]]@image$bcs
            # Calculate subclusters.
            sample_names <- list(sample.name)
            for (i in 1:length(sample_names)) {
              # Get cluster names.
              unique_clusters <- unique(bcs[[i]][cluster.column])
              unique_clusters <- sort(unique_clusters[!is.na(unique_clusters)])
              # Create empty dataframe with Cluster, row, and col.
              bcs[[i]] <- bcs[[i]] %>%
                tibble::add_column(!!(subcluster.column) := NA)
              # Calculate the subclusters for each cluster and add it to the dataframe.
              for (cluster in unique_clusters) {
                if (method == 'mclust'){
                  mc <- bcs[[i]] %>%
                    dplyr::filter(!!as.symbol(cluster.column) == cluster) %>%
                    dplyr::select(row, col) %>%
                    mclust::Mclust()
                  subclusters <- mc$classification
                } else if (method == 'dbscan') {
                  db <- bcs[[i]] %>%
                    dplyr::filter(!!as.symbol(cluster.column) == cluster) %>%
                    dplyr::select(row, col) %>%
                    dbscan::dbscan(eps=eps, minPts=minPts)
                  subclusters <- db$cluster
                }
                bcs[[i]][which(bcs[[i]][cluster.column] == cluster), ][subcluster.column] <-
                  paste(cluster, subclusters, sep='-')
              }
            }
            slot(object, 'assays')[[assay.type]]@image$bcs <- bcs
            if (assay.type=='ST'){
              metadata <- slot(object, 'assays')[['ST']]@sampTab
              clusters <- bcs[[1]][, c('barcode', subcluster.column)]
              metadata <- merge(metadata, clusters,
                                by.x="sample_id", by.y="barcode",
                                all=TRUE)
              metadata <- metadata[which(!is.na(metadata[[cluster.column]])), ]
              rownames(metadata) <- metadata$sample_id
              slot(object, 'assays')[['ST']]@sampTab <- metadata
            }
            return(object)
          }
)



# Centroids ---------------------------------------------------------------

#' Calculate cluster centroids
#'
#' Calculate spatial centroids for the clusters or the subclusters.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param cluster.column character; the name of the column where the clustering
#' information is stored.
#' @param cluster.type character; the column where the clusters are
#' indicated. It could be either 'Cluster' or 'Subcluster'.
#'
#' @return list; the centroids dataframe for each tissue
#' sample. Each dataframe has the columns cluster, row, col, imagerow, imagecol.
#' The row and col are the centroid coordinates while the imagerow and imagecol
#' are the coordinates for plotting, as scaled in the bcs list.
#'
#' @export
#' @docType methods
#' @rdname calculateCentroids-methods
setGeneric("calculateCentroids", function(object, assay.type='ST',
                                          sample.name='Sample1', cluster.column,
                                          cluster.type=c('Cluster', 'Subcluster'))
  standardGeneric("calculateCentroids"))
#' @rdname calculateCentroids-methods
#' @aliases calculateCentroids
setMethod("calculateCentroids",
          signature="CellRouter",
          definition=function(object, assay.type, sample.name,
                              cluster.column, cluster.type){
            cluster.type <- match.arg(cluster.type)
            bcs <- slot(object, 'assays')[[assay.type]]@image$bcs
            sample_names <- list(sample.name)
            # Calculate centroids.
            centroids <- list()
            for (i in 1:length(sample_names)) {
              # Get cluster names.
              unique_clusters <- unique(bcs[[i]][cluster.column])
              unique_clusters <- sort(unique_clusters[!is.na(unique_clusters)])
              # Create empty dataframe with Cluster, row, and col.
              if (cluster.type == 'Cluster'){
                centroids[[i]] <- data.frame(matrix(ncol=5, nrow=0))
                names(centroids[[i]]) <- c('Cluster', 'row', 'col',
                                           'imagerow', 'imagecol')
              } else if (cluster.type == 'Subcluster') {
                centroids[[i]] <- data.frame(matrix(ncol=6, nrow=0))
                names(centroids[[i]]) <- c('Cluster', 'Subcluster',
                                           'row', 'col', 'imagerow', 'imagecol')
              }
              # Calculate the centroid for each cluster and add it to the dataframe.
              for (cluster in unique_clusters) {
                c_row <- mean(bcs[[i]][which(bcs[[i]][cluster.column] == cluster), ]$row)
                c_col <- mean(bcs[[i]][which(bcs[[i]][cluster.column] == cluster), ]$col)
                # imagerow and imagecol have been already scaled to image resolution.
                i_row <- bcs[[i]][which(bcs[[i]]['row'] == round(c_row)), 'imagerow'][1]
                i_col <- bcs[[i]][which(bcs[[i]]['col'] == round(c_col)), 'imagecol'][1]
                if (cluster.type == 'Cluster'){
                  centroids[[i]][nrow(centroids[[i]]) + 1, ] <-
                    c(cluster, c_row, c_col, i_row, i_col)
                } else if (cluster.type == 'Subcluster') {
                  centroids[[i]][nrow(centroids[[i]]) + 1, ] <-
                    c(strsplit(cluster, '-')[[1]][1], cluster,
                      c_row, c_col, i_row, i_col)
                }
              }
              centroids[[i]]$row <- as.numeric(centroids[[i]]$row)
              centroids[[i]]$col <- as.numeric(centroids[[i]]$col)
              centroids[[i]]$imagerow <- as.numeric(centroids[[i]]$imagerow)
              centroids[[i]]$imagecol <- as.numeric(centroids[[i]]$imagecol)
            }
            names(centroids) <- sample_names
            slot(object, 'assays')[[assay.type]]@image$
              centroids[[cluster.type]] <- centroids
            return(object)
        }
)



# Distance between clusters -----------------------------------------------

#' Calculate distance matrix
#'
#' Calculate distance matrix for the clusters in each tissue sample.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param cluster.type character; the column where the clusters are
#' indicated. It could be either 'Cluster' or 'Subcluster'.
#' @param spot_distance numeric; the distance between the center of the spots.
#' The default is 100 micrometers. If this parameters is 1, the distance
#' matrix will be relative to the rows and columns of the spots.
#' @param normalize boolean; normalize the distances dividing by the
#' total size of the tissue.
#'
#' @return the CellRouter object.
#'
#' @export
#' @docType methods
#' @rdname calculateDistanceMatrix-methods
setGeneric("calculateDistanceMatrix", function(object, assay.type='ST',
                                               sample.name='Sample1',
                                               cluster.type=c('Cluster',
                                                              'Subcluster'),
                                               spot_distance=100,
                                               normalize=TRUE)
  standardGeneric("calculateDistanceMatrix"))
#' @rdname calculateDistanceMatrix-methods
#' @aliases calculateDistanceMatrix
setMethod("calculateDistanceMatrix",
          signature="CellRouter",
          definition=function(object, assay.type, sample.name, cluster.type,
                              spot_distance, normalize){
            cluster.type <- match.arg(cluster.type)
            centroids <- slot(object, 'assays')[[assay.type]]@image$
              centroids[[cluster.type]]
            distance_matrices <- list()
            sample_names <- list(sample.name)
            for (i in 1:length(sample_names)) {
              if (normalize){
                # Select positions from sample.
                positions <- centroids[[i]] %>% dplyr::select(imagerow, imagecol)
                rownames(positions) <- centroids[[i]][[cluster.type]]
                # Calculate distance in pixels and normalize.
                # The maximum distance is the diagonal.
                bcs <- slot(object, 'assays')[[assay.type]]@image$bcs
                max_dist <- sqrt((bcs[[i]]$height[[1]]^2) + (bcs[[i]]$width[[1]])^2)
                distance_matrices[[i]] <- as.matrix(dist(positions))/(max_dist)
              } else {
                # Select positions from sample.
                positions <- centroids[[i]] %>% dplyr::select(row, col)
                rownames(positions) <- centroids[[i]][[cluster.type]]
                # Calculate distance matrix and convert to actual distance.
                # The spot distance is divided by 2 because there is one space
                # between spots in the same line since the rows/cols are not
                # aligned.
                distance_matrices[[i]] <- as.matrix(dist(positions)) * (spot_distance/2)
              }
            }
            names(distance_matrices) <- sample_names
            slot(object, 'assays')[[assay.type]]@image$
              distances[[cluster.type]] <- distance_matrices
            return(object)
          }
)



# Visualization -----------------------------------------------------------

#' Plot spatial clusters
#'
#' Plot clusters in tissue image.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param cluster.column character; the name of the column where the clustering
#' information is stored.
#' @param colors.column character; the name of the column where the color of each
#' cluster is stored.
#' @param annotate_centroids boolean; annotate the cluster centroids.
#' @param point.size numeric; the size of the spot in the image.
#' @param title character; figure title.
#'
#' @return ggplot2; plot.
#'
#' @export
#' @import ggplot2
plotSpatialClusters <- function(object, assay.type='ST',
                                sample.name='Sample1',
                                cluster.column='population',
                                colors.column='colors',
                                annotate_centroids=FALSE,
                                point.size=1.75, title=NULL){
  if (is.null(title)){
    title <- paste0('Clusters of ', sample.name)
  }
  bcs <- slot(object, 'assays')[[assay.type]]@image$bcs[[sample.name]]
  # To annotate the centroids.
  cluster.type <- 'Cluster'
  centroids <- slot(object, 'assays')[[assay.type]]@image$
    centroids[[cluster.type]]
  # cluster_colors <- unique(bcs[[colors.column]])
  # Excluded the NAs for better plotting.
  cluster_colors <- na.exclude(unique(slot(object, 'assays')[[assay.type]]@sampTab[[colors.column]]))
  unique_clusters <- na.exclude(unique(slot(object, 'assays')[[assay.type]]@sampTab[[cluster.column]]))
  names(cluster_colors) <- unique_clusters
  # For code reusing.
  original_bcs_columns <- c('barcode', 'tissue', 'row', 'col',
                            'imagerow', 'imagecol', 'height', 'width',
                            cluster.column)
  bcs <- bcs[, original_bcs_columns]
  names(bcs)[names(bcs) == cluster.column] <- 'Cluster'
  bcs <- list(bcs)
  names(bcs) <- sample.name
  # Plot data.
  bcs_merge <- dplyr::bind_rows(bcs, .id="sample")
  images_tibble <- slot(object, 'assays')[[assay.type]]@image$images_tibble
  sample_names <- list(sample.name)
  plots <- list()
  for (i in 1:length(sample_names)) {
    plots[[i]] <- bcs_merge %>%
      dplyr::filter(sample == sample_names[i]) %>%
      dplyr::filter(tissue == "1") %>%
      dplyr::filter(Cluster %in% unique_clusters) %>%
      ggplot(aes(x=imagecol, y=imagerow, fill=factor(Cluster))) +
      geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape=21, colour="black", size=point.size, stroke=0.1)+
      coord_cartesian(expand=FALSE)+
      scale_fill_manual("Cluster", values=cluster_colors)+
      xlim(0,max(bcs_merge %>%
                   dplyr::filter(sample == sample_names[i]) %>%
                   dplyr::select(width)))+
      ylim(max(bcs_merge %>%
                 dplyr::filter(sample == sample_names[i]) %>%
                 dplyr::select(height)),0)+
      xlab("") +
      ylab("") +
      ggtitle(title)+
      labs(fill="Cluster")+
      guides(fill=guide_legend(override.aes=list(size=3)))+
      theme_set(theme_bw(base_size=10))+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background=element_blank(),
            axis.line=element_line(colour="black"),
            axis.text=element_blank(),
            axis.ticks=element_blank())
    # Insert centroids.
    # x is col and y is row.
    if (annotate_centroids){
      plots[[i]] <- plots[[i]] +
        annotate("text", color='white', size=4, fontface=2,
                 y=centroids[[i]][, 'imagerow'],
                 x=centroids[[i]][, 'imagecol'],
                 label=centroids[[i]][, cluster.type])
    }
    # filename_sample=paste0(filename, '_', sample_names[[i]], '.png')
    # ggsave(filename=filename_sample, plot=plots[[i]])
  }
  # if (display){
  #   cowplot::plot_grid(plotlist=plots)
  # }
  return(plots[[i]])
}



#' Plot selected spatial clusters
#'
#' Plot selected clusters in tissue image.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param cluster.column character; the name of the column where the clustering
#' information is stored.
#' @param colors.column character; the name of the column where the color of each
#' cluster is stored.
#' @param selected_clusters list; a vector containing the name of the clusters
#' that will be plotted for each tissue sample.
#' @param facets boolean; plot selected clusters in different facets.
#' @param annotate_centroids boolean; annotate the cluster centroids.
#' @param num.cols numeric; the number of columns in the output figure.
#' @param point.size numeric; the size of the spot in the image.
#' @param title character; figure title.
#'
#' @return ggplot2; plot.
#'
#' @export
#' @import ggplot2
plotSelectedClusters <- function(object, assay.type='ST', sample.name='Sample1',
                                 cluster.column='population',
                                 colors.column='colors', selected_clusters,
                                 facets=FALSE, annotate_centroids=FALSE,
                                 num.cols=2, point.size=1.75,
                                 title=NULL){

  if (is.null(title)){
    title <- paste0('Selected clusters of ', sample.name)
  }
  bcs <- slot(object, 'assays')[[assay.type]]@image$bcs[[sample.name]]
  # To annotate the centroids.
  cluster.type <- 'Cluster'
  centroids <- slot(object, 'assays')[[assay.type]]@image$
    centroids[[cluster.type]]
  # cluster_colors <- unique(bcs[[colors.column]])
  cluster_colors <- unique(slot(object, 'assays')[[assay.type]]@sampTab[[colors.column]])
  names(cluster_colors) <- unique(slot(object, 'assays')[[assay.type]]@sampTab[[cluster.column]])
  # For code reusing.
  original_bcs_columns <- c('barcode', 'tissue', 'row', 'col',
                            'imagerow', 'imagecol', 'height', 'width',
                            cluster.column)
  bcs <- bcs[, original_bcs_columns]
  names(bcs)[names(bcs) == cluster.column] <- 'Cluster'
  bcs <- list(bcs)
  names(bcs) <- sample.name
  # Plot data.
  bcs_merge <- dplyr::bind_rows(bcs, .id="sample")
  images_tibble <- slot(object, 'assays')[[assay.type]]@image$images_tibble
  sample_names <- list(sample.name)
  plots <- list()
  for (i in 1:length(sample_names)) {
    cluster_colors <- cluster_colors[selected_clusters[[i]]]
    plots[[i]] <- bcs_merge %>%
      dplyr::filter(sample == sample_names[i]) %>%
      dplyr::filter(tissue == "1") %>%
      dplyr::filter(Cluster %in% selected_clusters[[i]]) %>%
      ggplot(aes(x=imagecol,y=imagerow,fill=factor(Cluster))) +
      geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape=21, colour="black", size=point.size, stroke=0.1)+
      coord_cartesian(expand=FALSE)+
      scale_fill_manual("Cluster", values=cluster_colors)+
      xlim(0,max(bcs_merge %>%
                   dplyr::filter(sample == sample_names[i]) %>%
                   dplyr::select(width)))+
      ylim(max(bcs_merge %>%
                 dplyr::filter(sample == sample_names[i]) %>%
                 dplyr::select(height)),0)+
      xlab("") +
      ylab("") +
      ggtitle(title)+
      labs(fill="Cluster")+
      guides(fill=guide_legend(override.aes=list(size=3)))+
      theme_set(theme_bw(base_size=10))+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background=element_blank(),
            axis.line=element_line(colour="black"),
            axis.text=element_blank(),
            axis.ticks=element_blank())
    # Insert centroids.
    # x is col and y is row.
    if (annotate_centroids){
      plots[[i]] <- plots[[i]] +
        annotate("text", color='white', size=4, fontface=2,
                 y=centroids[[i]][which(centroids[[i]][[cluster.type]] %in%
                                          selected_clusters[[i]]),'imagerow'],
                 x=centroids[[i]][which(centroids[[i]][[cluster.type]] %in%
                                          selected_clusters[[i]]),'imagecol'],
                 label=selected_clusters[[i]])
    }
    # Plot facets
    if (facets){
      plots[[i]] <- plots[[i]] + facet_wrap(~Cluster, ncol=num.cols)
    }
  }
  return(plots[[i]])
}



#' Plot subclusters
#'
#' Plot subclusters from selected clusters.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param cluster.column character; the name of the column where the
#' clustering information is stored.
#' @param subcluster.column character; the name of the column where the
#' subclustering information is stored.
#' @param selected_clusters list; a vector containing the name of the clusters
#' that will be plotted for each tissue sample.
#' @param annotate_centroids boolean; annotate the subcluster centroids.
#' @param point.size numeric; the size of the spot in the image.
#' @param title character; figure title.
#'
#' @return ggplot2; plot.
#'
#' @export
#' @import ggplot2
plotSpatialSubclusters <- function(object, assay.type='ST',
                                   sample.name='Sample1',
                                   cluster.column='population',
                                   subcluster.column='Subpopulation',
                                   selected_clusters,
                                   annotate_centroids=FALSE,
                                   point.size=1.75,
                                   title=NULL){
  if (is.null(title)){
    title <- paste0('Selected subclusters of ', sample.name)
  }
  bcs <- slot(object, 'assays')[[assay.type]]@image$bcs[[sample.name]]
  cluster.type <- 'Subcluster'
  centroids <- slot(object, 'assays')[[assay.type]]@image$
    centroids[[cluster.type]]
  # cluster_colors <- unique(bcs[[colors.column]])
  # For code reusing.
  original_bcs_columns <- c('barcode', 'tissue', 'row', 'col',
                            'imagerow', 'imagecol', 'height', 'width',
                            cluster.column, subcluster.column)
  bcs <- bcs[, original_bcs_columns]
  names(bcs)[names(bcs) == cluster.column] <- 'Cluster'
  names(bcs)[names(bcs) == subcluster.column] <- 'Subcluster'
  bcs <- list(bcs)
  names(bcs) <- sample.name
  # Plot data.
  bcs_merge <- dplyr::bind_rows(bcs, .id="sample")
  images_tibble <- slot(object, 'assays')[[assay.type]]@image$images_tibble
  sample_names <- list(sample.name)
  plots <- list()
  for (i in 1:length(sample_names)) {
    plots[[i]] <- bcs_merge %>%
      dplyr::filter(sample == sample_names[i]) %>%
      dplyr::filter(tissue == "1") %>%
      dplyr::filter(Cluster %in% selected_clusters[[i]]) %>%
      ggplot(aes(x=imagecol,y=imagerow,fill=factor(Subcluster))) +
      geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(shape=21, colour="black", size=point.size, stroke=0.1)+
      coord_cartesian(expand=FALSE)+
      scale_fill_manual(values=c(RColorBrewer::brewer.pal(name="Dark2", n=8),
                                   RColorBrewer::brewer.pal(name="Paired", n=8)))+
      xlim(0,max(bcs_merge %>%
                   dplyr::filter(sample == sample_names[i]) %>%
                   dplyr::select(width)))+
      ylim(max(bcs_merge %>%
                 dplyr::filter(sample == sample_names[i]) %>%
                 dplyr::select(height)),0)+
      xlab("") +
      ylab("") +
      ggtitle(title)+
      labs(fill="Subcluster")+
      guides(fill=guide_legend(override.aes=list(size=3)))+
      theme_set(theme_bw(base_size=10))+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background=element_blank(),
            axis.line=element_line(colour="black"),
            axis.text=element_blank(),
            axis.ticks=element_blank())
    # Insert centroids.
    # x is col and y is row.
    if (annotate_centroids){
      plots[[i]] <- plots[[i]] +
        annotate("text", color='white', size=4, fontface=2,
                 y=centroids[[i]][which(centroids[[i]][['Cluster']] %in%
                                            selected_clusters[[i]]),'imagerow'],
                 x=centroids[[i]][which(centroids[[i]][['Cluster']] %in%
                                            selected_clusters[[i]]),'imagecol'],
                 label=centroids[[i]][which(centroids[[i]][['Cluster']] %in%
                                                selected_clusters[[i]]),
                                        cluster.type])
    }
  }
  return(plots[[i]])
}


#' Plot spatial expression
#'
#' Plot spatial expression for selected genes.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param cluster.column character; the name of the column where the
#' clustering information is stored.
#' @param selected_clusters list; a vector containing the name of the clusters
#' that will be plotted for each tissue sample.
#' @param genes character; a vector containing the name of the genes.
#' that will be plotted for each tissue sample.
#' @param point.size numeric; the size of the spot in the image.
#'
#' @return list; ggplot2 plots.
#'
#' @export
plotSpatialExpression <- function(object, assay.type='ST', sample.name='Sample1',
                                  cluster.column='population',
                                  selected_clusters=NULL,
                                  genes, point.size=1.75){
  myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
  # Get metadata.
  bcs <- slot(object, 'assays')[[assay.type]]@image$bcs[[sample.name]]
  # For code reusing.
  bcs <- list(bcs)
  names(bcs) <- sample.name
  # Plot data.
  bcs_merge <- dplyr::bind_rows(bcs, .id="sample")
  images_tibble <- slot(object, 'assays')[[assay.type]]@image$images_tibble
  sample_names <- list(sample.name)
  for (i in 1:length(sample_names)) {
    bcs_merge_sample <- bcs_merge %>%
      dplyr::filter(sample == sample_names[i]) %>%
      dplyr::filter(tissue == "1")
    # Select clusters.
    if (is.null(selected_clusters)){
      selected_clusters <- unique(bcs_merge_sample[[cluster.column]])
      selected_clusters <- list(selected_clusters)
    }
    bcs_merge_sample <- bcs_merge_sample[bcs_merge_sample[[cluster.column]] %in% selected_clusters[[i]], , drop=FALSE]
    # The addInfo funciton remove some barcodes from the expression matrix,
    # so it is necessary to ensure that both the expression matrix and the
    # metadata information have the same barcodes.
    names_exp <- colnames(slot(object, 'assays')[[assay.type]]@ndata)
    names_bcs <- bcs_merge_sample[['barcode']]
    mutual_bcs <- intersect(names_exp, names_bcs)
    # Get metadata.
    bcs_merge_sample <- bcs_merge_sample[which(bcs_merge_sample[['barcode']] %in%
                                                 mutual_bcs), ]
    # Get expression data.
    exp_mtx <- slot(object, 'assays')[[assay.type]]@ndata[genes, mutual_bcs,
                                                          drop=FALSE]
    exp_mtx <- exp_mtx[, bcs_merge_sample[['barcode']], drop=FALSE]
    # Gene plots.
    plots <- list()
    dfs <- data.frame()
    for(gene in genes){
      expr <- exp_mtx[gene, ]
      bcs_merge_sample$GENE <- as.numeric(expr)
      # bcs_merge_sample$gene <- gene
      dfs <- rbind(dfs, bcs_merge_sample)
      # Plots
      plots[[gene]] <- ggplot(dfs, aes(x=imagecol, y=imagerow, fill=GENE)) +
        geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
        geom_point(shape=21, colour="black", size=point.size, stroke=0.1)+
        coord_cartesian(expand=FALSE)+
        scale_fill_gradientn(colours=myPalette(100))+
        xlim(0, max(bcs_merge %>%
                     dplyr::filter(sample==sample_names[i]) %>%
                     dplyr::select(width)))+
        ylim(max(bcs_merge %>%
                   dplyr::filter(sample==sample_names[i]) %>%
                   dplyr::select(height)),0)+
        xlab("") +
        ylab("") + ggtitle(gene)
      theme_set(theme_bw(base_size=10))+
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              axis.line=element_line(colour="black"),
              axis.text=element_blank(),
              axis.ticks=element_blank())
    }
  }
  return(plots)
}



#' Plot spatial expression
#'
#' Plot spatial expression for selected genes.
#'
#' @param object the CellRouter object.
#' @param assay.type character; the type of data to use.
#' @param sample.name character; the name of the tissue sample.
#' @param cluster.column character; the name of the column where the
#' clustering information is stored.
#' @param selected_clusters list; a vector containing the name of the clusters
#' that will be plotted for each tissue sample.
#' @param genelist character vector; genes to show.
#' @param point.size numeric; the size of the spot in the image.
#'
#' @return list; ggplot2 plots.
#'
#' @export
#' @docType methods
#' @rdname plotSpatialScore-methods
setGeneric("plotSpatialScore", function(object, assay.type='ST',
                                        sample.name='Sample1',
                                        cluster.column='population',
                                        selected_clusters=NULL, genelist,
                                        point.size=1.75)
  standardGeneric("plotSpatialScore"))
#' @rdname plotSpatialScore-methods
#' @aliases plotSpatialScore
setMethod("plotSpatialScore",
          signature="CellRouter",
          definition=function(object, assay.type, sample.name, cluster.column,
                              selected_clusters, genelist, point.size){
            # Get image data.
            image_data <-  slot(object, 'assays')[[assay.type]]@image$images_tibble
            # Get sampTab.
            sampTab <- slot(object, 'assays')[[assay.type]]@sampTab
            # Color palette.
            myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
            # Get scores.
            if (is.null(selected_clusters)){
              selected_clusters <- unique(sampTab[[cluster.column]])
              selected_clusters <- list(selected_clusters)
            }
            scores <- sampTab[which(sampTab[[cluster.column]] %in% selected_clusters[[1]]), , drop=FALSE]
            #x <- 2^(object@ndata[genelist, rownames(scores)])-1
            x <- as.data.frame(t(sampTab[, genelist, drop=FALSE]))
            plots <- list()
            # Get barcodes and merge for location of spots.
            bcs <- slot(object, 'assays')[[assay.type]]@image$bcs[[sample.name]]
            bcs <- bcs %>%
              dplyr::filter(tissue == "1")
            scores <- merge(scores, bcs, by.x="sample_id", by.y="barcode",
                            all.x=TRUE)
            # Order the expression data and the scores in the same way.
            x <- x[ , scores[['sample_id']]]
            # Plot.
            dfs <- data.frame()
            for(gene in genelist){
              expr <- x[gene, ]
              scores$GENE <- as.numeric(expr)
              scores$gene <- gene
              dfs <- rbind(dfs, scores)
              # Plot data.
              p1 <- ggplot(dfs, aes(x=imagecol, y=imagerow, fill=GENE)) +
                geom_spatial(data=image_data, aes(grob=grob), x=0.5, y=0.5)+
                geom_point(shape=21, colour="black", size=point.size, stroke=0.1)+
                coord_cartesian(expand=FALSE)+
                scale_fill_gradientn("Expression", colours=myPalette(100))+
                # scale_colour_gradientn("Relative expression",
                #                        colours=c("midnightblue","white", "orange")) +
                xlim(0, 600) + ylim(600,0) +
                xlab("") +
                ylab("") + ggtitle(gene)
              theme_set(theme_bw(base_size=10))+
                theme(panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      panel.background=element_blank(),
                      axis.line=element_line(colour="black"),
                      axis.text=element_blank(),
                      axis.ticks=element_blank())
              p1 <- p1 + theme(legend.position="bottom")
              plots[[gene]] <- p1
            }
            return(plots)
          }
)
