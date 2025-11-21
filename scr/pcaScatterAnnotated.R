pcaScatterAnnotated <- function(pca_results, targets, annotation) {
  
  #' Generate a 2D scatter plot of the first two PCs, annotated with a given column in the phenotype data
  #'
  #' @param pca_results the matrix of eigen values output from the PCA cacluation
  #' @param targets the targets matrix
  #' @param annotation The the name of the column in the targets dataframe used to annotate plot
  
  annotation <- enquo(annotation)  
  
  ## compute the variance explained by each PC for the n PCs calculated
  pca_var_explained <- (pca_results$sdev^2/sum(pca_results$sdev^2))
  
  ## Scatter plots of pairs of principle components
  pca_results_merged <- pca_results$x %>%
    as_tibble(rownames = "analysisID") %>%
    merge(targets, by="analysisID",  sort=FALSE)    # may need to left join this...
  
  # calculate centroids
  (centroids.df <- pca_results_merged %>% dplyr::select(PC1, PC2, !!annotation))
  colnames(centroids.df) <- c("PC1", "PC2", "anno")
  centroids <- aggregate(cbind(PC1, PC2)~anno,centroids.df,mean)
  
  pca_plot <- ggplot(pca_results_merged, aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(!!annotation)),size=4) +
    theme_bw(base_size=32) +
    labs(x=paste0("PC1: ",round(pca_var_explained[1]*100,1),"%"),
         y=paste0("PC2: ",round(pca_var_explained[2]*100,1),"%"),
         color = element_blank()) +
    theme(legend.position="bottom",
          legend.margin=ggplot2::margin(-5,-5,-5,-5),
          legend.box.margin=ggplot2::margin(-5,-5,-5,-5),
          axis.text=element_text(size=18),
          axis.title=element_text(size=18),
          plot.title=element_text(size=18),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18)) +
    stat_ellipse(aes(color=(!!annotation)), level = 0.85) +
    geom_point(data = centroids, size = 4, shape = 4) +
    scale_colour_manual(values = cbp1)

  
  return(pca_plot)
}