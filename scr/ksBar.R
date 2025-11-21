ksBar <- function(targets, metrics, dist.outliers, annotation){
  
  #' plots horizontal bar chart of K-S statistic
  #'
  #' @param targets sample information dataframe
  #' @param metrics metrics to plot, here the  K-S statistic
  #' @param dist.outliers list of outliers identified based on K-S distance
  #' @param annotation annotation variable to colour code bars
  
  annotation <- enquo(annotation)
  
  # ks.plot.data <- data.frame(targets %>% dplyr::select(analysisID, all_of((!!annotation))),
  #                            'k_stat' = metrics$k_stat)
  
  ks.plot.data <- left_join(targets %>% dplyr::select(analysisID, all_of((!!annotation))),
                            metrics, 
                            by = join_by(analysisID))
  
  
  ks.barchart <- ggplot(ks.plot.data, aes(x=reorder(analysisID, k_stat), y = k_stat, fill=get(names(ks.plot.data)[2]))) +
    geom_bar(stat='identity') +
    coord_flip() +
    xlab("sampleID") + 
    ylab("k-s statistic") + 
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.title = element_text(size = 18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = dist.outliers$threshold) +
    scale_fill_manual(values = cbp1) +
    scale_y_continuous(limits = c(-0.04 , (max(ks.plot.data$k_stat) *1.05))) +
    
    # annotate outliers with an asterisk and labels in red
    geom_point(data = ks.plot.data %>% dplyr::filter(k_stat > dist.outliers$threshold),
               aes(x = analysisID, y = -0.01), shape = 8, show.legend = F) + 
    annotate(geom = "text",
             x = ks.plot.data %>% dplyr::filter(k_stat > dist.outliers$threshold) %>% .$analysisID,
             y = -0.02,
             hjust = 1,
             size = 3,
             colour = "red",
             label = ks.plot.data %>% dplyr::filter(k_stat > dist.outliers$threshold) %>% .$analysisID) + 
    
    # annotate the rest of the samples in black
    annotate(geom = "text",
             x = ks.plot.data %>% dplyr::filter(k_stat  <= dist.outliers$threshold) %>% .$analysisID,
             y = -0.02,
             hjust = 1,
             size = 3,
             colour = "black",
             label = ks.plot.data %>% dplyr::filter(k_stat  <= dist.outliers$threshold) %>% .$analysisID)
    
  
  return(ks.barchart)
}

# previous version without sample annotation

#' ksBar <- function(targets, metrics, dist.outliers, annotation){
#'   
#'   #' plots horizontal bar chart of K-S statistic
#'   #'
#'   #' @param targets sample information dataframe
#'   #' @param metrics metrics to plot, here the  K-S statistic
#'   #' @param dist.outliers list of outliers identified based on K-S distance
#'   #' @param annotation annotation variable to colour code bars
#'   
#'   annotation <- enquo(annotation)
#'   
#'   ks.plot.data <- data.frame(targets %>% dplyr::select(analysisID, all_of((!!annotation))),
#'                              'k_stat' = metrics$k_stat)
#'   
#'   
#'   ks.barchart <- ggplot(ks.plot.data, aes(x=reorder(analysisID, k_stat), y = k_stat, fill=get(names(ks.plot.data)[2]))) +
#'     geom_bar(stat='identity') +
#'     coord_flip() +
#'     ggtitle("Ranked K-S Statistic") + 
#'     xlab("sampleID") + 
#'     ylab("k-s statistic") + 
#'     theme(legend.title = element_blank(),
#'           axis.text.y = element_blank(),
#'           axis.ticks.y = element_blank()) +
#'     geom_hline(yintercept = dist.outliers$threshold) +
#'     scale_fill_manual(values = cbp1)
#'   
#'   
#'   return(ks.barchart)
#' }