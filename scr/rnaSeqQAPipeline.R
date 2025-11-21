rnaSeqQAPipeline <- function(norm.counts, vst.counts, sample.info, annotation_list, out_dir, ma_topx, n_samp = 5000, pqual = 300){
  
  # L1 distance: calculate scaled L1 distance matrix, identify outliers based on Tukey's distance
  
  mat.dist.matrix <- l1DistanceMatrix(norm.counts) # [find a faster function for this]
  
  l1.metrics <- l1DistOutliers(mat.dist.matrix)
  
  # PCA generate PCA results using the full set of features with vst counts, save pca scree plot
  
  pca_results <- stats::prcomp(t(vst.counts), center = TRUE, scale = FALSE, rank = 9)
  
  ggsave(filename = "pca.scree.pdf", plot = pcaScreePlot(pca_results), dpi = pqual, path = out_dir)
  
  # KS distance: calculate k-s statistic for each sample, outliers based on 1.5 IQR
 
  k.stats <- calcKS(norm.counts)
  
  dist.outliers <- ksOutliers(k.stats)
  
  # Hoeffding's D statistic
  
  d.stats <- dStatOutliers(vst.counts)
  
  # Sex check - pca plot of y linked genes only
  
  y.linked.ensembl <- read.table(here("resources", "y_linked_ensemblID.txt"), sep="", header=F) %>% deframe

  y.linked.ensembl <- as.list(rownames(vst.counts)[rownames(vst.counts) %in% y.linked.ensembl]) # ensure all y-linked genes are in the dataset
  
  vst.y.linked <- subset(vst.counts, row.names(vst.counts) %in% y.linked.ensembl)

  pca.y.linked <- stats::prcomp(t(vst.y.linked), center = TRUE, scale = FALSE, rank = 9)
  
  pca.sex.check <- pcaScatterAnnotated(pca_results = pca.y.linked, targets = sample.info, annotation = sex)
  
  ggsave(filename = "pca.sex.check.pdf", plot = pca.sex.check, dpi = pqual, path = out_dir)
  
  # loop through annotations to create the plots
  
  for(anno in annotation_list){
    
    l1.bar <- l1DistBar(targets = sample.info, metrics = l1.metrics, annotation = anno)
    ggsave(filename = paste("l1.bar", anno, "pdf", sep = "."), plot = l1.bar , dpi = pqual, path = out_dir)
    
    pca.scatter <- pcaScatterAnnotated(pca_results = pca_results, targets = sample.info, annotation = !!sym(anno))
    ggsave(filename = paste("pca", anno, "pdf", sep="."), plot = pca.scatter, dpi = pqual, path = out_dir)
    
    ks.bar.plot <- ksBar(targets = sample.info, metrics = k.stats, dist.outliers = dist.outliers, annotation = anno)
    ggsave(filename = paste("ks.bar", anno, "pdf", sep = "."), plot = ks.bar.plot, dpi = pqual, path = out_dir)
    
    d.stat.bar <- dStatBar(targets = sample.info, metrics = d.stats, annotation = anno)
    ggsave(filename = paste("d.stat.bar", anno, "pdf", sep = "."), plot = d.stat.bar, dpi = pqual, path = out_dir)
    
    ma.plot <- maPlots(vst.counts = vst.counts, targets = sample.info, metrics = d.stats, annotation = anno, top_x= ma_topx, bottom_x= ma_topx, subsample = n_samp)
    ggsave(filename = paste("ma.plot", anno, "pdf", sep = "."), plot = ma.plot, dpi = pqual, path = out_dir)
    
  }
  
  outliers <- list('l1' = l1.metrics$l1.outliers, 'ks' = dist.outliers$outliers, 'dstat' = d.stats$dstat.outliers)
  
  return(outliers)
}