DoDESEQ <- function(counts.long.cis, min.median.counts, jbkgrs, jbait, mindist, maxdist, remove.samp, time.only = TRUE, verbose=FALSE){
  
  counts.sub <- subset(counts.long.cis, bait %in% jbait & abs(rel.pos) > mindist & abs(rel.pos) < maxdist)
  
  
  # Add sample.i ------------------------------------------------------------
  
  counts.sub <- counts.sub %>%
    group_by(time, bkgr, rel.pos) %>%
    mutate(sample.i = seq(length(counts)))
  
  # filter out noisy fragments
  counts.sub <- counts.sub %>%
    group_by(time, bkgr, rel.pos) %>%
    filter(median(counts) > min.median.counts)
  
  counts.sub$coord <- as.character(counts.sub$rel.pos)
  
  # takes ~7min
  start <- Sys.time()
  
  counts.sub <- counts.sub %>%
    mutate(samp.name = paste(time, bkgr, sample.i, sep = ";"))
  
  # counts.sub$samp.name <- as.factor(mapply(function(time, bk, samp.i) paste(time, bk, samp.i, sep = ";"), 
  #                                as.character(counts.sub$time), 
  #                                as.character(counts.sub$bkgr), 
  #                                as.character(counts.sub$sample.i)))
  if (verbose) print(Sys.time() - start)
  
  if (length(unique(counts.sub$bkgr)) == 3){
    counts.sub$bkgr <- factor(x = as.character(counts.sub$bkgr), levels = c("WT", "KO", "Kidney"))
    counts.sub$samp.name <- factor(x = as.character(counts.sub$samp.name), 
                                   levels = c(paste("ZT08;WT", seq(4), sep = ";"), 
                                              paste("ZT20;WT", seq(4), sep = ";"), 
                                              paste("ZT08;KO", seq(3), sep = ";"), 
                                              paste("ZT20;KO", seq(3), sep = ";"), 
                                              paste("ZT08;Kidney", seq(4), sep = ";"), 
                                              paste("ZT20;Kidney", seq(4), sep = ";")))
  } else {
    counts.sub$bkgr <- factor(x = as.character(counts.sub$bkgr), levels = c("WT", "KO"))
    counts.sub$samp.name <- factor(x = as.character(counts.sub$samp.name), 
                                   levels = c(paste("ZT08;WT", seq(4), sep = ";"), 
                                              paste("ZT20;WT", seq(4), sep = ";"), 
                                              paste("ZT08;KO", seq(3), sep = ";"), 
                                              paste("ZT20;KO", seq(3), sep = ";")))
  }
  
  if (!is.null(remove.samp)){
    if (verbose)  print("Removing samples")
    if (verbose)  print(remove.samp)
    counts.sub <- subset(counts.sub, !samp.name %in% remove.samp)
  }
  
  counts.mat <- dcast(counts.sub, formula = coord ~ samp.name, value.var = "counts")
  
  rownames(counts.mat) <- counts.mat$coord; counts.mat$coord <- NULL
  
  counts.mat <- counts.mat[which(rowSums(counts.mat) > 0), ]
  
  
  if (!is.null(jbkgrs)){
    if (verbose)  print(paste("Keep bkgr:", jbkgrs))
    counts.mat <- counts.mat[, grepl(jbkgrs, colnames(counts.mat))]
    if (verbose)  print(colnames(counts.mat))
  }
  
  des.mat <- data.frame(time = sapply(colnames(counts.mat), function(cname) strsplit(cname, ";")[[1]][[1]]),
                        bkgr = sapply(colnames(counts.mat), function(cname) strsplit(cname, ";")[[1]][[2]]))
  
  # Run DESeq2 --------------------------------------------------------------
  
  # ceiling rounds up 0.5s
  if (time.only){
    dds <- DESeqDataSetFromMatrix(countData = round(counts.mat),
                                  colData = des.mat,
                                  design = ~ time)  
  } else {
    dds <- DESeqDataSetFromMatrix(countData = round(counts.mat),
                                  colData = des.mat,
                                  design = ~ time + bkgr) 
  }
  
  start <- Sys.time()
  dds <- DESeq(dds)
  res <- results(dds)
  resMLE <- results(dds, addMLE=TRUE)
  # vsd <- varianceStabilizingTransformation(dds)
  if (verbose) print(Sys.time() - start)
  
  
  norm.long.res <- data.frame(rel.pos = as.numeric(rownames(res)), baseMean = res$baseMean, log2FC = res$log2FoldChange, padj = res$padj, pvalue = res$pvalue)
  norm.long <- melt(counts(dds, normalized=TRUE), value.name = "norm.counts", varnames = c("rel.pos", "samp"))
  
  # add times, genotypes, samp.i
  norm.long$time <- sapply(as.character(norm.long$samp), function(s) strsplit(s, ";")[[1]][[1]])
  norm.long$genotype <- sapply(as.character(norm.long$samp), function(s) strsplit(s, ";")[[1]][[2]])
  norm.long$samp.i <- sapply(as.character(norm.long$samp), function(s) strsplit(s, ";")[[1]][[3]])
  
  
  norm.long.avg <- norm.long %>%
    group_by(rel.pos, time, genotype) %>%
    summarise(norm.counts.avg = mean(norm.counts)) %>%
    mutate(norm.counts.avg.log10 = log10(norm.counts.avg + 0.5))
  
  return(list(norm.long = norm.long, norm.long.res = norm.long.res, norm.long.avg = norm.long.avg))
}