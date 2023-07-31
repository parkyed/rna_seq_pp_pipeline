# multset Jaccard filter threshold calculation

msJaccFilter <- function(x, condition, s.min, s.max, s.len){
  
  # create a sequence of thresholds to test, sequence generated based on equally spaced log values
  
  (t.seq <- exp(seq(log(s.min), log(s.max), length = s.len)))
  
  # Function to calculate the Multiset Jaccard Index' over all samples
  Jaccard = function(mat){
    row.sums = rowSums(mat)
    inter = sum(row.sums == ncol(mat))
    union = sum(row.sums > 0)
    return(inter/union)
  }
  
  # Calculate vector of the minimum Multiset Jaccard Index over the classes, for each value of t
  
  ms.jac = sapply(t.seq, function(t){
    filt.counts = x >=t
    group.jac = sapply(unique(condition), function(cond){
      Jaccard(filt.counts[,condition==cond])
    })
    return(min(group.jac))
  })
  
  # calculate the value of t that maximises the multiset Jaccard
  (t.hold <- t.seq[which.max(ms.jac)])
  
  # plot the threshold value against the value of the Multiset Jaccard index to visualise
  
  ms.jacc.plot <- ggplot(data=data.frame(t = t.seq, jacc = ms.jac)) +
    geom_line(aes(x=t, y=jacc)) +
    geom_hline(yintercept = max(ms.jac), lty=2,col='gray') +
    geom_point(aes(x=t.seq[which.max(ms.jac)], y=max(ms.jac)), col="red", cex=6, pch=1) +
    xlab("Low Count Threshold") + 
    ylab("Multiset Jaccard Index")
  
  return(list('threshold' = t.hold, 'ms.jacc.plot' = ms.jacc.plot))
}
