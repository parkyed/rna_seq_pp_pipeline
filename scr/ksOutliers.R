ksOutliers <- function(k.stats){
  
  #' identify outlier values that lie more than 1.5IQR beyond 25th or 75th percentiles
  #'
  #' @param k.stats vector of K-S statistics
  
  ks.outlier.threshold <- boxplot.stats(k.stats$k_stat)$stats[5]
  
  dist.outliers <- which(tibble::deframe(k.stats) > ks.outlier.threshold)
  
  return(list(threshold = ks.outlier.threshold, outliers = dist.outliers))
}