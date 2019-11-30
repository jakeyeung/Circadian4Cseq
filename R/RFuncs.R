# Jake Yeung
# Date of Creation: 2019-06-13
# File: ~/projects/4cseq_part2/functions/RFuncs.R
# Load R functions


ConvertToSingleDay <- function(dat, time.cname = "time", exprs.cname = "exprs",
                               geno.cname = "geno", exper.cname = "experiment",
                               tissue.cname = "tissue", gene.cname = "gene", transcript.cname = NA, use.se=FALSE){
  library(lazyeval)
  if (any(dat$time < 0)){
    warning("Function cannot handle negative times")
  }
  mutate_call = lazyeval::interp(~jtime - 24 * floor( jtime / 24 ), jtime = as.name(time.cname))
  dat.sub <- dat %>%
    filter_(paste(time.cname, ">=", 24)) %>%
    mutate_(.dots = setNames(list(mutate_call), time.cname))
  
  if (!use.se){
    sem.cname <- "sem"
    sem_call = lazyeval::interp(~sd( exprs ), exprs = as.name(exprs.cname))
  } else {
    se.cname <- "se"
    sem.cname <- "se"
    sem_call = lazyeval::interp(~mean( se ), exprs = as.name(se.cname))
  }
  signal_call = lazyeval::interp(~mean ( exprs ), exprs = as.name(exprs.cname))
  
  if (is.na(transcript.cname)){
    dat.new <- bind_rows(dat.sub, subset(dat, time < 24)) %>%
      group_by_(gene.cname, geno.cname, exper.cname, tissue.cname, time.cname) %>%
      summarise_(.dots = setNames(list(sem_call, signal_call), c(sem.cname, exprs.cname)))
  } else {
    dat.new <- bind_rows(dat.sub, subset(dat, time < 24)) %>%
      group_by_(gene.cname, geno.cname, exper.cname, tissue.cname, time.cname, transcript.cname) %>%
      summarise_(.dots = setNames(list(sem_call, signal_call), c(sem.cname, exprs.cname)))
  }
  return(dat.new)
}

PlotGeneTissuesWTKO2 <- function(dat, timelabel="ZT", jtitle="", split.by="geno", ncols = 2, center = FALSE, convert.linear = FALSE, jsize = 24, single.day=FALSE, pretty.geno.names=FALSE){
  if (pretty.geno.names){
    dat$geno <- gsub("BmalKO", "Bmal1 KO", dat$geno)
    dat$geno <- gsub("SV129", "WT", dat$geno)
    dat$geno <- factor(dat$geno, levels = c("WT", "Bmal1 KO"))
  }
  if (convert.linear){
    dat$exprs <- (2 ^ dat$exprs) - 1
    jylab <- "TPM"
  } else {
    jylab <- "Log2 mRNA Abundance"
  }
  # split by geno or tissue
  if (center){
    dat <- dat %>%
      group_by(tissue, geno) %>%
      mutate(exprs = scale(exprs, center = TRUE, scale = FALSE))
  }
  if (single.day){
    dat <- ConvertToSingleDay(dat)
  }
  # m <- ggplot(dat, aes(x = time, colour = tissue, linetype = geno, y = exprs)) +
  m <- ggplot(dat, aes(x = time, linetype = geno, y = exprs)) +
    geom_point() + geom_line() + xlab(timelabel) + ylab(jylab) +
    theme_bw(jsize) + ggtitle(jtitle) +
    theme(aspect.ratio = 1, legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  if (split.by == "geno"){
    m <- m + facet_wrap(~geno, ncol = ncols)
  } else if (split.by == "tissue"){
    m <- m + facet_wrap(~tissue, ncol = ncols)
  } else {
    warning("Split by must be geno or tissue")
  }
  if (single.day){
    # m <- m + geom_errorbar(data = dat, aes(ymin=exprs-sem, ymax=exprs+sem, colour = tissue), linetype = "solid", size=0.5, width=0.5)
    m <- m + geom_errorbar(data = dat, aes(ymin=exprs-sem, ymax=exprs+sem), linetype = "solid", size=0.5, width=0.5)
  }
  return(m)
}


PlotSignalLR2 <- function(jbait, counts.long, pseudo.low, mindist, jtitle=NULL, do.facet=TRUE, ltype="solid", convert.kb = TRUE){
  if (is.null(jtitle)){
    jtitle <- paste(jbait, "Pseudo:", pseudo.low)
  }
  jsub <- subset(counts.long, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  if (convert.kb){
    jsub <- jsub %>%
      mutate(pos = pos / 1000)
    jxlab <- "Position relative to bait [kb]"
  } else {
    jxlab <- "Position relative to bait"
  }
  AvsPosT <- ggplot(jsub, aes(x = pos, y = A, colour = time)) +
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = A, colour = time), linetype = ltype) +
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = A, colour = time), linetype = ltype) +
    theme(aspect.ratio = 0.25, panel.grid.minor = element_blank()) + 
    geom_hline(aes(yintercept=0)) +
    theme_bw() + ggtitle(jtitle) + geom_vline(aes(xintercept=0), linetype="dotted") +
    theme(legend.position = "bottom") +
    xlab(jxlab) +
    ylab("Log10 Signal")
    # scale_y_continuous(breaks = NULL)
  if (do.facet){
    AvsPosT <- AvsPosT + facet_grid(genotype ~ tissue)
  }
  return(AvsPosT)
}

PlotZscorePvalLR2 <- function(jbait, counts.delt, pseudo.low, mindist, jylim.zscore = c(-5.5, 5.5), jylim.pval = c(-10, 10), convert.kb = TRUE, jtitle = ""){
  jsub <- subset(counts.delt, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  if (convert.kb){
    jsub <- jsub %>%
      mutate(pos = pos / 1000)
    jxlab <- "Position relative to bait [kb]"
  } else {
    jxlab <- "Position relative to bait"
  }
  if (jylim.zscore == "auto"){
    # print(range(jsub$zscore.row))
    jylim.zscore <- c(floor(min(-jsub$zscore.row)), ceiling(max(-jsub$zscore.row)))
    # print(jylim.zscore)
  } 
  if (jylim.pval == "auto"){
    # print(range(-log10(jsub$pval)))
    jylim.pval <- c(floor(min(-log10(jsub$pval) * sign(-jsub$zscore.row))), ceiling(max(-log10(jsub$pval) * sign(-jsub$zscore.row))))
    # print(jylim.pval)
  }
  Zscores <- ggplot(jsub, aes(x = pos, y = -zscore.row)) +
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = -zscore.row)) +
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = -zscore.row)) +
    facet_grid(genotype ~ tissue) + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0), linetype="dotted") + 
    # theme(aspect.ratio = 0.25, panel.grid.minor = element_blank()) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab(jxlab) + ylab("Zscore") + 
    theme_bw() +
    ggtitle(jtitle) + 
    ylim(jylim.zscore)

  Pvals <- ggplot(jsub, aes(x = pos, y = -log10(pval.row) * sign(-zscore.row))) +
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = -log10(pval.row) * sign(-zscore.row))) +
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = -log10(pval.row) * sign(-zscore.row))) +
    facet_grid(genotype ~ tissue) + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0), linetype="dotted") + 
    # theme(aspect.ratio = 0.25, panel.grid.minor = element_blank()) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme_bw() +
    xlab(jxlab) + ylab(expression("Signed -log[10](Pval)")) + 
    ylim(jylim.pval)
	multiplot(Zscores, Pvals)
}

