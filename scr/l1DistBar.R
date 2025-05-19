l1DistBar <- function(targets, metrics, annotation){
  
  #' plots horizontal bar chart of L1 distance
  #'
  #' @param targets sample information dataframe
  #' @param metrics metrics to plot, here the L1 distance
  #' @param annotation annotation variable to colour code bars
  
  l1.plot.data <- data.frame(targets %>% dplyr::select(analysisID, all_of((!!annotation))),
                             'SumDist' = metrics$sum.dist)
  
  l1.barchart <-  ggplot(l1.plot.data, aes(x=reorder(analysisID, SumDist), y = SumDist, fill=get(names(l1.plot.data)[2]))) +
    geom_bar(stat='identity') +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.title = element_text(size = 18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("SampleID") + 
    ylab("Sum of L1 Distance") + 
    scale_y_continuous(limits = c(-7 , (max(l1.plot.data$SumDist) *1.05))) +
    geom_hline(yintercept = metrics$l1.outlier.threshold) +
    scale_fill_manual(values = cbp1) +
    geom_point(data = l1.plot.data %>% dplyr::filter(SumDist > metrics$l1.outlier.threshold),
               aes(x = analysisID, y = -1), shape = 8, show.legend = F) + 
    annotate(geom = "text",
             x = l1.plot.data %>% dplyr::filter(SumDist > metrics$l1.outlier.threshold) %>% .$analysisID,
             y = -2,
             hjust = 1,
             size = 3,
             label = l1.plot.data %>% dplyr::filter(SumDist > metrics$l1.outlier.threshold) %>% .$analysisID)
  
  
  return(l1.barchart)
}

# previous version without sample annotation

#' l1DistBar <- function(targets, metrics, annotation){
#'   
#'   #' plots horizontal bar chart of L1 distance
#'   #'
#'   #' @param targets sample information dataframe
#'   #' @param metrics metrics to plot, here the L1 distance
#'   #' @param annotation annotation variable to colour code bars
#'   
#'   l1.plot.data <- data.frame(targets %>% dplyr::select(analysisID, all_of((!!annotation))),
#'                              'SumDist' = metrics$sum.dist)
#'   
#'   l1.barchart <- ggplot(l1.plot.data, aes(x=reorder(analysisID, SumDist), y = SumDist, fill=get(names(l1.plot.data)[2]))) +
#'     geom_bar(stat='identity') +
#'     coord_flip() +
#'     theme(legend.title = element_blank(),
#'           axis.text.y = element_blank(),
#'           axis.ticks.y = element_blank()) +
#'     xlab("SampleID") + 
#'     ylab("Sum of L1 Distance") + 
#'     geom_hline(yintercept = metrics$l1.outlier.threshold) +
#'     scale_fill_manual(values = cbp1)
#'   
#'   return(l1.barchart)
#' }
