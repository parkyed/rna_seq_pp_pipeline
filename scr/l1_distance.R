# Calculate the normalised L1 distance matrix

l1_distance_matrix <- function(mat.dist) {
  
  #' Calculate L1-distance between each row (sample) in the input matrix. Use the as.matrix function to include zeros and upper triangle 
  #' scale the output to between zero and 1 by dividing by the largest value in the matrix
  #' return a scaled matrix of L1 distances between columns in the input matrix
  #'
  #' @param mat.dist The matrix of log2 transformed, DESeq2 median-ratio normalised counts
  
  mat.dist.matrix <- as.matrix(dist(t(mat.dist), method='manhattan'))
  
  mat.dist.matrix = mat.dist.matrix/max(mat.dist.matrix)
  
  return(mat.dist.matrix)
}

# identify outliers based on Tukey's distance, and return the total manhattan distance, the tukey threhold and list of outliers

l1DistOutliers <- function(mat.dist.matrix){
  
  sum.dist <- rowSums(mat.dist.matrix)
  
  l1.outlier.threshold <- boxplot.stats(sum.dist)$stats[5] %>% unname
  
  l1.outliers <- which(sum.dist > l1.outlier.threshold)
  
  return(list('sum.dist' = sum.dist, 'l1.outlier.threshold' = l1.outlier.threshold, 'l1.outliers' = l1.outliers))
}


# plot bar charts for alternative 

l1DistBarChart <- function(targets, metrics, annotation){
  
  l1.plot.data <- data.frame(targets %>% dplyr::select(analysisID, all_of((!!annotation))),
                             'SumDist' = metrics$sum.dist)
  
  l1.barchart <- ggplot(l1.plot.data, aes(x=reorder(analysisID, SumDist), y = SumDist, fill=get(names(l1.plot.data)[2]))) +
    geom_bar(stat='identity') +
    coord_flip() +
    theme(legend.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    xlab("SampleID") + 
    ylab("Sum of L1 Distance") + 
    geom_hline(yintercept = metrics$l1.outlier.threshold)
  
  return(l1.barchart)
}


# Create a heatmap based on the L1 distances

l1_distance_heatmap <- function(mat.dist, targets, annotation) {
  
  #' Generate a labelled heatmap from a matrix of log2 transformed, normalised counts, annotated by a factor of interest
  #'
  #' @param mat.dist The matrix of log2 transformed, DESeq2 median-ratio normalised counts
  #' @param targets The targets dataframe
  #' @param annotation The factor in the targets dataframe used to annotate the heatmap
  
  # enquote annotation
  annotation <- enquo(annotation)  
  
  # calculate scaled L1 distance matrix
  mat.dist.matrix <- l1_distance_matrix(mat.dist)
  
  # create dataframe of annotations and move analysisID to rownames
  targets.filtered <- dplyr::select(targets, analysisID, (!!annotation))
  hm.anno.df <- targets.filtered
  rownames(hm.anno.df) <- hm.anno.df[,1]
  hm.anno.df[,1] <- NULL
  
  #anno.df <-  data.frame(targets.filtered, row.names = targets$analysisID)
  if (!identical(rownames(hm.anno.df), rownames(mat.dist.matrix))) { stop() }
  
  # plot heatmap
  l1_heatmap <- pheatmap(mat.dist.matrix,
                         cluster_rows = T, cluster_cols = T, show_rownames = T, show_colnames = T,
                         border_color = NA, scale = "none", ylab = "", main = "L1 distances between samples", 
                         col = colorRampPalette(brewer.pal(9, "YlOrRd"))(20),
                         annotation_row = hm.anno.df,
                         fontsize = 6, cellwidth = 7, cellheight = 7,
                         height = 6,
                         width = 6)
  
  return(l1_heatmap)
}

### PREVIOUS l1 bar chart function - doesn't scale to alternative annotations as repeats the expensive l1 distance matrix calculation ###

# Create a barchart of the sum of L1 distances for each sample from all other samples

#' l1_distance_barchart <- function(mat.dist, targets, annotation) {
#'   
#'   #' Generate a horizontal bar chart of the sum of the L1 distances for each sample from all other samples, annotated by a factor of interest
#'   #'
#'   #' @param mat.dist The matrix of log2 transformed, DESeq2 median-ratio normalised counts
#'   #' @param targets The targets dataframe
#'   #' @param annotation The factor in the targets dataframe used to annotate the chart
#'   
#'   # enquote annotation
#'   annotation <- enquo(annotation)  
#'   
#'   # calculate scaled L1 distance matrix
#'   mat.dist.matrix <- l1_distance_matrix(mat.dist)
#'   
#'   # create dataframe of annotations given unquoted annotation arguement
#'   targets.filtered <- dplyr::select(targets, analysisID, all_of((!!annotation)))
#'   
#'   # sum distances between each array and all other arrays 
#'   l1_dist_df <- data.frame(targets.filtered, 'SumDist' = rowSums(mat.dist.matrix))
#'   
#'   # identify outlier values that lie more than 1.5IQR beyond 25th or 75th percentiles
#'   l1_dist_outlier_threshold <- boxplot.stats(l1_dist_df$SumDist)$stats[5]
#'   l1_dist_outliers <- filter(l1_dist_df, SumDist > l1_dist_outlier_threshold)
#'   
#'   # plot horizontal bar chart
#'   l1_barchart <- ggplot(l1_dist_df, aes(x=reorder(analysisID, SumDist), y = SumDist, fill=get(names(l1_dist_df)[2]))) +
#'     geom_bar(stat='identity') +
#'     coord_flip() +
#'     theme(legend.title = element_blank(),
#'           axis.text.y = element_blank(),
#'           axis.ticks.y = element_blank()) +
#'     xlab("SampleID") + 
#'     ylab("Sum of L1 Distance") + 
#'     geom_hline(yintercept = l1_dist_outlier_threshold)
#'   
#'   output <- list('dist_outliers' = l1_dist_outliers, 'dist_barchart' = l1_barchart)
#'   
#'   return(output)
#' }
