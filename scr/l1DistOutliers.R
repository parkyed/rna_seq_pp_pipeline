l1DistOutliers <- function(mat.dist.matrix){
  
  #' identify outliers based on Tukey's distance, and return the total Manhattan distance, the tukey threhold and list of outliers
  #'
  #' @param mat.dist.matrix scaled matrix of L1 (Manhattan) distances between columns
  
  sum.dist <- rowSums(mat.dist.matrix)
  
  l1.outlier.threshold <- boxplot.stats(sum.dist)$stats[5] %>% unname
  
  l1.outliers <- which(sum.dist > l1.outlier.threshold)
  
  return(list('sum.dist' = sum.dist, 'l1.outlier.threshold' = l1.outlier.threshold, 'l1.outliers' = l1.outliers))
  
}