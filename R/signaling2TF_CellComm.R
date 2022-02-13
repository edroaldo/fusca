#' Signaling to transcription factors.
#'
#' Select the paths in the specified interaction and plot the activity score.
#' Save the plot to file and return the data frame used to generate it.
#'
#' @param npaths list; the ranked pathways and the number of number of targets
#' in each interaction as calculated by the rankpaths function.
#' @param interaction character; the name of the selected interaction.
#' @param num.pathways numeric; the number of pathways considered in the
#' signaling.
#'
#' @return list; a path and interaction data frame, and the plot.
#'
#' @export
signaling2TF <- function(npaths, interaction, num.pathways = 10){
  df <- npaths$npaths[[interaction]]
  df <- df[1:num.pathways, ]
  df$path <- gsub("(s->|->t)", "", as.vector(df$path))
  a <- sapply(strsplit(df$path, split='->', fixed=TRUE), function(x){x[1]})
  b <- sapply(strsplit(df$path, split='->', fixed=TRUE), function(x){x[length(x)]})
  df$path <- paste(a, '->...->', b, sep='')
  #df <- df[!duplicated(as.vector(df$path)),]
  df$path <- factor(as.vector(df$path), levels=rev(unique(as.vector(df$path))))
  df$interaction <- interaction
  g <- ggplot2::ggplot(df, ggplot2::aes(x = path, y = score, fill=score)) +
    ggplot2::geom_bar(stat = 'identity', colour='black') +
    ggplot2::theme_bw() +
    ggplot2::xlab("Signaling-to-transcription paths") +
    ggplot2::ylab("Activity score") +
    ggplot2::theme(legend.position="none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = ggplot2::rel(1),
                                                       angle = 00, hjust = 1)) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(size = 0.5, linetype = "solid",
                                              colour = "black"),
          axis.line.y = ggplot2::element_line(size = 0.5, linetype = "solid",
                                              colour = "black")) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradientn(colours = c('white', 'orange', 'brown')) +
    ggplot2::ggtitle(interaction)
  return(list(df = df, plot = g))
}
