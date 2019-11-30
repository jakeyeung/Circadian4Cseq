DeseqNormalizeLong <- function(counts.long, dist.min, dist.max, plot.diags.dir = FALSE){
  # Input:
  #   counts.long by bait (can include WT, KO, Kidney, Time)
  #   should include sample.i i and samp.name ZT08;WT;i
  # dist.min:
  #   minimum distance from bait in order to be included in analysis
  # dist.max
  #  max dist from bait to be included in analysis
  jbait <- unique(counts.long$bait)
  print(jbait)
  print(counts.long)
  if (length(jbait) > 1) warning("Warning: more than one bait found")
  
  
  # redo factor or samp.name to prevent a column of NA if no kidney was performed
  
  samp.names <- unique(counts.long$samp.name)
  # new factors
  new.factors <- levels(samp.names)[levels(samp.names) %in% samp.names[which(!is.na(samp.names))]]
  
  counts.mat <- dcast(counts.long, formula = coord ~ samp.name, value.var = "counts")
  rownames(counts.mat) <- counts.mat$coord; counts.mat$coord <- NULL
  
  # remove rows where no counts across all replicates
  counts.mat <- counts.mat[which(rowSums(counts.mat) > 0), ]
  
  print(colnames(counts.mat))
  des.mat <- data.frame(time = sapply(colnames(counts.mat), 
                                      function(cname) strsplit(cname, ";")[[1]][[1]]),
                        bkgr = sapply(colnames(counts.mat), 
                                      function(cname) strsplit(cname, ";")[[1]][[2]]))
  
  dds <- DESeqDataSetFromMatrix(countData = round(counts.mat),
                                colData = des.mat,
                                design = ~ time + bkgr)
  
  vsd <- varianceStabilizingTransformation(dds)
  counts.vsd <- assays(vsd)@listData[[1]]
  
  vsd.long <- melt(counts.vsd)
  colnames(vsd.long) <- c("rel.pos", "samp.name", "counts.vsd")
  
  #   print(head(vsd.long$samp.name))
  vsd.long$rel.pos <- as.numeric(vsd.long$rel.pos)
  
  #   print(as.character(vsd.long$samp.name[[1]]))
  #    print(strsplit(as.character(vsd.long$samp.name[[1]]), ";")[[1]])
  vsd.long$time <- sapply(as.character(vsd.long$samp.name), function(s) strsplit(s, ";")[[1]][[1]])
  vsd.long$bkgr <- sapply(as.character(vsd.long$samp.name), function(s) strsplit(s, ";")[[1]][[2]])
  vsd.long$samp.i <- sapply(as.character(vsd.long$samp.name), function(s) strsplit(s, ";")[[1]][[3]])
  
  if (plot.diags.dir != FALSE){
    jbait <- vsd.long
    plots.diags.path <- file.path(plot.diags.dir, paste0(jbait, ".pdf"))
    pdf(plot.diags.path)
    res <- results(dds)
    resMLE <- results(dds, addMLE=TRUE)
    vsd <- varianceStabilizingTransformation(dds)
    plotMA(data.frame(resMLE$baseMean, resMLE$lfcMLE, resMLE$padj < 0.1))
    plotMA(res, main = "Shrunken LFC")
    
    vsd.long.mean <- vsd.long %>%
      group_by(bkgr, time, rel.pos) %>%
      summarise(counts.vsd.mean = mean(counts.vsd))
    
    ggplot(subset(vsd.long.mean, abs(rel.pos) < 100000), 
           aes(x = rel.pos, y = counts.vsd.mean, colour = time)) + 
      geom_line() + facet_wrap(~bkgr) + ggtitle(paste(jbait, mindist, maxdist))
    
    dev.off()
  }
  return(vsd.long)
}