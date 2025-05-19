l1DistanceMatrix <- function(mat.dist) {
  
  #' Calculate L1-distance between each row (sample) in the input matrix. Use the as.matrix function to include zeros and upper triangle 
  #' scale the output to between zero and 1 by dividing by the largest value in the matrix
  #' return a scaled matrix of L1 distances between columns in the input matrix
  #'
  #' @param mat.dist The matrix of log2 transformed, DESeq2 median-ratio normalised counts
  
  mat.dist.matrix <- as.matrix(dist(t(mat.dist), method='manhattan'))
  
  mat.dist.matrix = mat.dist.matrix/max(mat.dist.matrix)
  
  return(mat.dist.matrix)
}