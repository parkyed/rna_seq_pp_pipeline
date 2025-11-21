maPlots <- function(vst.counts, targets, metrics, annotation, top_x, bottom_x, subsample) {
  
  #' plots top x and bottom x MA plots, based on D statistic
  #'
  #' @param vst.counts matrix of log2 vst transformed counts, with row = ensembl gene IDs, col =  sample
  #' @param targets sample information dataframe
  #' @param metrics metrics used in the plot, here the D statistic
  #' @param annotation annotation variable to colour code bars
  #' @param top_x the n of plots with the highest D stat values
  #' @param bottom_x the n plots with the lowest D stat values
  #' @param subsample the size of the subsample of randomly selected transcripts to use to calculate M and A
  
  
  # Compute the M and A matrices and assign to variables
  ma_list <- maCalc(vst.counts, subsample = 10000)
  M=ma_list$M
  A=ma_list$A
  
  # select the top x and bottom x rows of the M and A matrices
  
  annotation <- enquo(annotation)
  
  dstat.plot.data <- data.frame(targets %>% dplyr::select(analysisID, all_of((!!annotation))),
                                'd_stat' = metrics$d.stats)
  
  colnames(dstat.plot.data)[2] <- "annotation"
  
  ma_stats_id_sel <- rbind(dstat.plot.data %>% slice_max(order_by = d_stat, n = top_x),
                           dstat.plot.data %>% slice_min(order_by = d_stat, n = bottom_x)) %>% 
    dplyr::mutate(group = c(rep('high', times=top_x), rep('low', times=bottom_x))) %>% 
    dplyr::arrange(-d_stat)
  
  M_sel <- M[,ma_stats_id_sel$analysisID]
  
  A_sel <- A[,ma_stats_id_sel$analysisID]
  
  # create data frame for the plots
  ma_plot_df <- data.frame(
    sample_id=rep(ma_stats_id_sel$analysisID,each=nrow(M_sel)),
    annotation = rep(ma_stats_id_sel$annotation,each=nrow(M_sel)),
    d_stat = round(rep(ma_stats_id_sel$d_stat, each = nrow(M_sel)),3),
    group = factor(rep(ma_stats_id_sel$group,each=nrow(M_sel)), levels = c('high', 'low')),
    M=as.numeric(M_sel),
    A=as.numeric(A_sel)
  )
  
  # plots
  ma_plots <- ggplot(ma_plot_df,aes(x=A, y=M, color=annotation)) +
    geom_point(alpha=0.1) +
    facet_wrap(~group + d_stat + sample_id, ncol=5) +
    scale_colour_manual(values = cbp1) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(size = 14, margin = margin())) +
    guides(colour = guide_legend(override.aes = list(size=7, alpha = 1))) # increase the size of the legend points, override the alpha above
  ma_plots
  
  return(ma_plots)
}
