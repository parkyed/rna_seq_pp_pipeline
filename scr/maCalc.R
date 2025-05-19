maCalc <- function(x, subsample) {
  
  #' computes M and A matrices using the counts for a sub-sample of randomly selected transcripts
  #' output: returns a list with M and A data matrices
  #'
  #' @param x matrix of log2 normalised counts (e.g. vst counts), with row = ensembl gene IDs, col =  sample
  #' @param subsample the size of the subsample of randomly selected transcripts to use to calculate M and A
  
  medArray = rowMedians(x, na.rm=TRUE)
  M =  x - medArray
  A = (x + medArray)/2
  if(nrow(M)>subsample) {
    set.seed(12345)
    sel = sample(nrow(M), subsample)
    sM = M[sel, ]
    sA = A[sel, ]
  } else {
    sM = M
    sA = A
  }
  output <- list(M=sM,A=sA) 
  
  return(output)
}