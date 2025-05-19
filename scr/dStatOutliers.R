dStatOutliers <- function(vst.counts, n_samp = 5000){
  
  #' computes the D statistic and returns outliers based on the statistic and Tukey's distance
  #'
  #' @param vst.counts matrix of log2 vst transformed counts, with row = ensembl gene IDs, col =  sample
  
  # calculate the M and A matrices
  
  ma_list <- maCalc(vst.counts, subsample = n_samp)
  M=ma_list$M
  A=ma_list$A
  
  # Compute the Hoeffding's statistic (Da) statistics, and identify outliers using Tukey's distance
  
  d.stats =  sapply(1:ncol(M), function(x){wdm::indep_test(x = A[,x],
                                                           y = M[,x],
                                                           method = 'hoeffding')$statistic})
  
  (dstat.outlier.threshold <- boxplot.stats(d.stats)$stats[5] %>% unname)
  
  dstat.outliers <- which(d.stats > dstat.outlier.threshold)
  
  return(list('d.stats' = d.stats, 'dstat.outlier.threshold' = dstat.outlier.threshold, 'dstat.outliers' = dstat.outliers))
  
}
