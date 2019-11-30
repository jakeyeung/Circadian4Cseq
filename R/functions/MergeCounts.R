MergeCounts <- function(counts.all, baitpair, baits.long){
  counts.sub <- subset(counts.all, bait %in% baitpair)
  baits.long.sub <- subset(baits.long, bait %in% baitpair)
  bait.ref <- baitpair[[1]]
  bait.pos <- hash(as.character(baits.long.sub$bait), baits.long.sub$pos)
  
  baits.long.sub$pos.rel <- mapply(function(bait, pos){
    return(pos - bait.pos[[bait.ref]])
  }, as.character(baits.long.sub$bait), baits.long.sub$pos)
  
  bait.pos.rel <- hash(as.character(baits.long.sub$bait), baits.long.sub$pos.rel)
  
  counts.sub$pos.rel <- mapply(function(bait, pos){
    return(pos + bait.pos.rel[[bait]])
  }, as.character(counts.sub$bait), counts.sub$pos)
  
  return(counts.sub)
}

MergeCountsWTKO <- function(jbait, counts.delt, baits.long, jtissue = "Liver", max.dist = Inf, show.plot=TRUE, jtitle = "",
                            jxlab="Position relative to bait", jshow.legend=TRUE, flip.y.axis=FALSE){
  bpair <- c(paste(jbait, "WT", sep="-"), paste(jbait, "KO", sep="-"))
  
  baitsub <- subset(counts.delt, tissue == jtissue & bait == jbait)
  
  baitsub$bait <- factor(sapply(baitsub$genotype, function(geno){
    return(paste(jbait, geno, sep = "-"))
  }), levels = bpair)
  
  # add fake info for baitlocs
  bait.locs.wtko <- data.frame()
  for (pair in bpair){
    chromo <- subset(baits.long, bait == jbait)$chromo
    pos <- subset(baits.long, bait == jbait)$pos
    bait.locs.wtko <- rbind(bait.locs.wtko, data.frame(bait=pair, chromo=chromo, pos=pos))
  }
  
  merged <- MergeCounts(baitsub, bpair, bait.locs.wtko)
  jsub <- subset(merged, abs(pos.rel) < max.dist)
  if (show.plot==FALSE) return(jsub)
  if (show.plot == "Zscore"){
    jplot <- PlotZscoresOverlay(jsub, jtitle = jtitle, jbait, xlab = jxlab, show.legend = jshow.legend, flip.y.axis = flip.y.axis)
  } else if (show.plot == "Pvalue"){
    jplot <- PlotPvaluesOverlay(jsub, jtitle = jtitle, jbait, xlab = jxlab, show.legend = jshow.legend, flip.y.axis = flip.y.axis)
  } else {
    warning("Defaulting to Zscore")
    # default to Zscore
    jplot <- PlotZscoresOverlay(jsub, jtitle = jtitle, jbait, xlab = jxlab, show.legend = jshow.legend, flip.y.axis = flip.y.axis)
  }
  return(jplot)
}

MergeCountsLivKidTemporal <- function(jbait, counts.delt, baits.long, tissues = c("Kidney", "Liver"), max.dist = Inf, show.plot=TRUE, jtitle = "",
                            jxlab="Position relative to bait", jshow.legend=TRUE, flip.y.axis=FALSE){
  bpair <- c(paste(jbait, tissues[[1]], sep="-"), paste(jbait, tissues[[2]], sep="-"))
  
  baitsub <- subset(counts.delt, tissue %in% tissues & bait == jbait & genotype == "WT")
  
  baitsub$bait <- factor(sapply(baitsub$tissue, function(tiss){
    return(paste(jbait, tiss, sep = "-"))
  }), levels = bpair)
  
  # add fake info for baitlocs
  bait.locs.wtko <- data.frame()
  for (pair in bpair){
    chromo <- subset(baits.long, bait == jbait)$chromo
    pos <- subset(baits.long, bait == jbait)$pos
    bait.locs.wtko <- rbind(bait.locs.wtko, data.frame(bait=pair, chromo=chromo, pos=pos))
  }
  
  merged <- MergeCounts(baitsub, bpair, bait.locs.wtko)
  jsub <- subset(merged, abs(pos.rel) < max.dist)
  if (show.plot==FALSE) return(jsub)
  if (show.plot == "Zscore"){
    jplot <- PlotZscoresOverlay(jsub, jtitle = jtitle, jbait, xlab = jxlab, show.legend = jshow.legend, flip.y.axis=flip.y.axis)
  } else if (show.plot == "Pvalue"){
    jplot <- PlotPvaluesOverlay(jsub, jtitle = jtitle, jbait, xlab = jxlab, show.legend = jshow.legend, flip.y.axis=flip.y.axis)
  } else {
    warning("Defaulting to Zscore")
    # default to Zscore
    jplot <- PlotZscoresOverlay(jsub, jtitle = jtitle, jbait, xlab = jxlab, show.legend = jshow.legend, flip.y.axis=flip.y.axis)
  }
  return(jplot)
}

MergeCountsLivKid <- function(jbait, counts.delt.lk, baits.long, zt="ZT08", geno="WT", max.dist = Inf, show.plot="Zscore", jtitle="", jxlab="Position relative to bait", jshow.legend=TRUE, flip.y.axis=FALSE, reorder.zt20vzt08=TRUE, in.kb=FALSE){
  
  if (!reorder.zt20vzt08){
    timepair <- c("ZT08", "ZT20")
  } else {
    timepair <- c("ZT20", "ZT08")
  }
  bpair <- c(paste(jbait, timepair[[1]], sep="-"), paste(jbait, timepair[[2]], sep="-"))
  
  baitsub <- subset(counts.delt.lk, genotype == geno & bait == jbait)
  
  baitsub$bait <- factor(sapply(as.character(baitsub$time), function(time){
    return(paste(jbait, time, sep = "-"))
  }), levels = bpair)
  
  # add fake info for baitlocs
  bait.locs.livkid <- data.frame()
  for (pair in bpair){
    chromo <- subset(baits.long, bait == jbait)$chromo
    pos <- subset(baits.long, bait == jbait)$pos
    bait.locs.livkid <- rbind(bait.locs.livkid, data.frame(bait=pair, chromo=chromo, pos=pos))
  }
  
  merged <- MergeCounts(baitsub, bpair, bait.locs.livkid)
  jsub <- subset(merged, abs(pos.rel) < max.dist)
  if (show.plot==FALSE) return(jsub)
  if (show.plot == "Zscore"){
    jplot <- PlotZscoresOverlay(jsub, jtitle, xlab=jxlab, show.legend=jshow.legend, flip.y.axis=flip.y.axis, in.kb=in.kb)
  } else if (show.plot == "Pvalue"){
    jplot <- PlotPvaluesOverlay(jsub, jtitle, xlab=jxlab, show.legend=jshow.legend, flip.y.axis=flip.y.axis, in.kb=in.kb)
  } else {
    # default to Zscore
    warning("Defaulting to Zscore")
    jplot <- PlotZscoresOverlay(jsub, jtitle, xlab=jxlab, show.legend=jshow.legend, flip.y.axis=flip.y.axis, in.kb=in.kb)
  }
  return(jplot)
}

PlotZscoresOverlay <- function(jsub, jtitle, xlab="Position relative to bait", ylab="Zscore", show.legend=TRUE, flip.y.axis=FALSE, in.kb=FALSE, 
                               jcols = c("black", "brown"),
                               ltypes = c("solid", "solid")){
  if (show.legend){
    leg.pos <- "bottom"
  } else {
    leg.pos <- "none"
  }
  if (flip.y.axis){
    jsub$zscore.row <- jsub$zscore.row * -1
  }
  if (in.kb){
    jsub$pos.rel <- jsub$pos.rel/1000
  }
  
  Zscores <- ggplot(jsub, aes(x = pos.rel, y = -zscore.row, linetype = bait, color = bait)) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos.rel, y = -zscore.row)) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos.rel, y = -zscore.row)) + 
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0), linetype="dotted") +
    theme_bw() + 
    ggtitle(jtitle) + 
    xlab(xlab) +
    ylab(ylab) + 
    theme(aspect.ratio = 1/8, legend.position = leg.pos) + 
    scale_color_manual(values = jcols) + 
    scale_linetype_manual(values = ltypes)
  return(Zscores)
}

PlotPvaluesOverlay <- function(jsub, jtitle, xlab="Position relative to bait", ylab="-log10(Pvalue)", show.legend=TRUE, flip.y.axis=FALSE, in.kb=FALSE,
                               jcols = c("black", "brown"), ltypes = c("solid", "solid")){
  if (show.legend){
    leg.pos <- "bottom"
  } else {
    leg.pos <- "none"
  }
  if (flip.y.axis){
    jsub$A.delta <- jsub$A.delta * -1 
  }
  if (in.kb){
    jsub$pos.rel <- jsub$pos.rel/1000
  }
  
  jsub$pval.sign <- mapply(function(delta, pval) sign(delta) * -log10(pval), jsub$A.delta, jsub$pval.row)  # -log10
  
  Pvals <- ggplot(jsub, aes(x = pos.rel, y = pval.sign, linetype = bait, color = bait)) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos.rel, y = pval.sign)) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos.rel, y = pval.sign)) + 
    geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0), linetype="dotted") +
    theme_bw() + 
    xlab(xlab) + 
    ylab(ylab) + 
    ggtitle(jtitle) + 
    theme(aspect.ratio = 1/8, legend.position = leg.pos) + 
    scale_color_manual(values = jcols) + 
    scale_linetype_manual(values = ltypes)
  return(Pvals)
}

