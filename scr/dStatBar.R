dStatBar <- function(targets, metrics, annotation) {
  
  #' plots horizontal bar chart of D statistic
  #'
  #' @param targets sample information dataframe
  #' @param metrics metrics to plot, here the D statistic
  #' @param annotation annotation variable to colour code bars
  
  annotation <- enquo(annotation)
  
  dstat.plot.data <- data.frame(targets %>% dplyr::select(analysisID, all_of((!!annotation))),
                                'd_stat' = metrics$d.stats)
  
  d.stat.barplot <- ggplot(dstat.plot.data, aes(x=reorder(analysisID, d_stat), y = d_stat, fill=get(names(dstat.plot.data)[2]))) +
    geom_bar(stat='identity') + 
    coord_flip() +
    xlab("sampleID") + 
    ylab("Hoeffding's D-statistic") + 
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.title = element_text(size = 18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = metrics$dstat.outlier.threshold) +
    scale_fill_manual(values = cbp1)  +
    scale_y_continuous(limits = c(-0.01 , (max(dstat.plot.data$d_stat) *1.05))) +
    
    # annotate outliers with an asterisk and labels in red
    geom_point(data = dstat.plot.data %>% dplyr::filter(d_stat > metrics$dstat.outlier.threshold),
               aes(x = analysisID, y = -0.003), shape = 8, show.legend = F) + 
    annotate(geom = "text",
             x = dstat.plot.data %>% dplyr::filter(d_stat > metrics$dstat.outlier.threshold) %>% .$analysisID,
             y = -0.005,
             hjust = 1,
             size = 3,
             colour = "red",
             label = dstat.plot.data %>% dplyr::filter(d_stat > metrics$dstat.outlier.threshold) %>% .$analysisID) +
    
    # annotate the rest of the samples in black
    annotate(geom = "text",
             x = dstat.plot.data %>% dplyr::filter(d_stat <= metrics$dstat.outlier.threshold) %>% .$analysisID,
             y = -0.005,
             hjust = 1,
             size = 3,
             colour = "black",
             label = dstat.plot.data %>% dplyr::filter(d_stat <= metrics$dstat.outlier.threshold) %>% .$analysisID)
  
  return(d.stat.barplot)
}


# previous version without sample annotation


#' dStatBar <- function(targets, metrics, annotation) {
#'   
#'   #' plots horizontal bar chart of D statistic
#'   #'
#'   #' @param targets sample information dataframe
#'   #' @param metrics metrics to plot, here the D statistic
#'   #' @param annotation annotation variable to colour code bars
#'   
#'   annotation <- enquo(annotation)
#'   
#'   dstat.plot.data <- data.frame(targets %>% dplyr::select(analysisID, all_of((!!annotation))),
#'                                 'd_stat' = metrics$d.stats)
#'   
#'   d.stat.barplot <- ggplot(dstat.plot.data, aes(x=reorder(analysisID, d_stat), y = d_stat, fill=get(names(dstat.plot.data)[2]))) +
#'     geom_bar(stat='identity') + 
#'     coord_flip() +
#'     ggtitle("Hoeffding's d-statistic on joint dist. of M and A values") + 
#'     xlab("sampleID") + 
#'     ylab("Hoeffding's d-statistic") + 
#'     theme(legend.title = element_blank(),
#'           axis.text.y = element_blank(),
#'           axis.ticks.y = element_blank()) +
#'     geom_hline(yintercept = metrics$dstat.outlier.threshold) +
#'     scale_fill_manual(values = cbp1)
#'   
#'   return(d.stat.barplot)
#' }