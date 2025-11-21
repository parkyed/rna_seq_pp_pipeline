

# create melted dataframe with all analysis IDs, and annotation added - box plot input

meltAnnotate <- function(vst.counts, sample.info, annotation){
  
  # filter targets file based on the column name in annotation argument
  annotation <- enquo(annotation)
  
  si.filtered <- dplyr::select(sample.info, analysisID, (!!annotation))
  
  if (!identical(colnames(vst.counts), si.filtered$analysisID)) { stop() }
  
  # melt dataframe or normalised counts, add the annotation column
  
  vst.melt <- as.data.frame(vst.counts) %>%               
    tibble::rownames_to_column("ensemblID")  %>%      
    melt(id.vars = 'ensemblID', value.name='count', variable.name=c('analysisID')) %>%  
    left_join(si.filtered, by = 'analysisID')
  
  return(vst.melt)
}

# create box plot from melted annotated dataframe

distBoxPlot <- function(logNormCount.melt){
  
  ggplot(logNormCount.melt, aes(x=reorder(analysisID, count, FUN = median),
                                y = count,
                                fill=get(names(logNormCount.melt)[4]))) + 
    geom_boxplot(outlier.shape = NULL) + 
    ggtitle("Box Plot of Distribution of vst. counts")+ 
    xlab("analysisID") + 
    coord_flip() +
    theme(legend.title = element_blank()) +
    ylab("vst counts")
  
}

# calculate k-s statistic for each sample

calcKS <- function(vst.counts){
  
  pooled.vst.counts = sample(c(vst.counts), 100000)
  
  ks.results <- apply(vst.counts, 2, stats::ks.test, y=pooled.vst.counts, exact = NULL)
  
  k.stats <- lapply(ks.results, function(x){x$statistic})
  
  k.stats <- lapply(k.stats, unname) %>% unlist() %>% as.data.frame() %>% rownames_to_column()
  
  colnames(k.stats) <- c('analysisID', 'k_stat')
  
  return(k.stats)
}


# identify outlier values that lie more than 1.5IQR beyond 25th or 75th percentiles

distOutliers <- function(k.stats){
  ks.outlier.threshold <- boxplot.stats(k.stats$k_stat)$stats[5]
  dist.outliers <- dplyr::filter(k.stats, k_stat > ks.outlier.threshold) %>% dplyr::arrange(desc(k_stat))
  return(list(threshold = ks.outlier.threshold, outliers = dist.outliers))
}

# plot bar chart of D statistic for each sample

ksBar <- function(sample.info, k.stats, dist.outliers, annotation){
  annotation <- enquo(annotation)
  sample.info %>%
    dplyr::left_join(k.stats, by = 'analysisID') %>% 
    ggplot(aes(x=reorder(analysisID, k_stat), y = k_stat, fill=(!!annotation))) +
    geom_bar(stat='identity') +
    ggtitle("Ranked K-S Statistic") + 
    xlab("sampleID") + 
    ylab("k-s statistic") + 
    coord_flip() +
    theme(legend.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    geom_hline(yintercept = dist.outliers$threshold)
}

# =======================================================================

## all in one function - OLD and very heavy

#' box_density_k_s <- function(logNormCount, targets, annotation){
#'   
#'   #' Generate a box plot, density plots, calculate K-S statistic for each sample and show in a horizontal bar. Identify outliers based on k-s.
#'   #'
#'   #' @param logNormCount The matrix of log2 transformed, DESeq2 median-ratio normalised counts
#'   #' @param targets The targets dataframe
#'   #' @param annotation The factor in the targets dataframe used to annotate the chart
#'   
#'   # filter targets file based on the column name in annotation argument
#'   annotation <- enquo(annotation)
#'   targets.filtered <- dplyr::select(targets, analysisID, (!!annotation))
#'   if (!identical(colnames(logNormCount), targets.filtered$analysisID)) { stop() }
#'   
#'   # melt dataframe or normalised counts, add the annotation column
#'   logNormCount.melt <- as.data.frame(logNormCount) %>%               
#'     tibble::rownames_to_column("ensemblID")  %>%      
#'     melt(id.vars = 'ensemblID', value.name='count', variable.name=c('analysisID')) %>%  
#'     left_join(targets.filtered, by = 'analysisID')
#'   
#'   # boxplot of the distribution of each sample
#'   dist_box_plot <- ggplot(logNormCount.melt, aes(x=reorder(analysisID, count, FUN = median), y = count, fill=get(names(logNormCount.melt)[4]))) + 
#'     geom_boxplot() + 
#'     ggtitle("Box Plot of Distribution of log2 normalised counts")+ 
#'     xlab("sampleID") + 
#'     coord_flip() +
#'     theme(legend.title = element_blank()) +
#'     ylab("log[2](count + 1)")
#'   
#'   # density plot of the distribution of each sample, split out by annotation
#'   dist_density_plot <- ggplot(logNormCount.melt, aes(x = count, fill=analysisID)) + 
#'     geom_density(alpha = 0.2, size = 1.25) + 
#'     ggtitle("Density Plot of log2 normalised counts") + 
#'     xlab("counts") + 
#'     ylab("density") + 
#'     facet_wrap(~ get(names(logNormCount.melt)[4])) + 
#'     theme(legend.position="none")
#'   
#'   # pool all the log counts
#'   pooled_logNormCount = c(logNormCount)                                             
#'   
#'   # loop over every column in the log2Normalised counts and perform k-s test
#'   # matricise this using purrr - get rid of the for loop!!
#'   k_stats <- character()
#'   for (i in 1:nrow(targets.filtered)){
#'     k_stat <- ks.test(logNormCount[,i], pooled_logNormCount, exact = FALSE)
#'     k_stats[i] <- k_stat$statistic
#'   }
#'   targets.filtered$k_stat <- as.numeric(k_stats)
#'   
#'   # identify outlier values that lie more than 1.5IQR beyond 25th or 75th percentiles
#'   ks_outlier_threshold <- boxplot.stats(targets.filtered$k_stat)$stats[5]
#'   dist_outliers <- filter(targets.filtered, k_stat > ks_outlier_threshold)
#'   
#'   # bar chart of D statistic for each sample
#'   dist_outlier_plot <- ggplot(targets.filtered, aes(x=reorder(analysisID, k_stat), y = k_stat, fill=get(names(logNormCount.melt)[4]))) +
#'     geom_bar(stat='identity') +
#'     ggtitle("Bar plot of k-s statistic by sample") + 
#'     xlab("sampleID") + 
#'     ylab("k-s statistic") + 
#'     coord_flip() +
#'     theme(legend.title = element_blank()) +
#'     geom_hline(yintercept = ks_outlier_threshold)
#'   
#'   output <- list('box_plot'= dist_box_plot, 'density_plot' = dist_density_plot, 'outliers' = dist_outliers, 'k_stat_bar' = dist_outlier_plot)
#'   
#'   return(output)
#' }