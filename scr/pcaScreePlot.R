pcaScreePlot <- function(pca_results){
  
  #' Generate a scree plot showing the % of variation in dataset explained by each principal component
  #'
  #' @param pca_results the matrix of eigen values output from the PCA cacluation
  
  ## compute the variance explained by each PC for the n PCs calculated
  pca_var_explained <- (pca_results$sdev^2/sum(pca_results$sdev^2))[1:ncol(pca_results$x)]
  
  ## Scree plot of the proportion of variance explained by each PC
  scree_plot <- data.frame('principal_component'=colnames(pca_results$x), 'variance_explained'= pca_var_explained) %>%
    ggplot(aes(x=principal_component, y=pca_var_explained)) +
    geom_bar(stat='identity') +
    ggtitle("Scree Plot - % Variance Explained")
  
  return(scree_plot)
  
}