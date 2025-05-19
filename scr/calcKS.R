calcKS <- function(vst.counts, n_samp = 5000){
  
  #' calculate k-s statistic for each sample
  #'
  #' @param vst.counts matrix vst transformed counts
  
  pooled.vst.counts <- sample(c(as.matrix(vst.counts)), n_samp)
  
  ks.results <- apply(vst.counts, 2, stats::ks.test, y=pooled.vst.counts, exact = NULL)
  
  k.stats <- lapply(ks.results, function(x){x$statistic})
  
  k.stats <- lapply(k.stats, unname) %>% unlist() %>% as.data.frame() %>% rownames_to_column()
  
  colnames(k.stats) <- c('analysisID', 'k_stat')
  
  return(k.stats)
}