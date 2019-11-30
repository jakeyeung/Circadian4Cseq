# PlotGeneAcrossTissues.R

PlotGeneLivKidKO <- function(dat, jtitle, add.vline=TRUE){
  if (missing(jtitle)) jtitle <- as.character(dat$gene[[1]])
  m1 <- ggplot(dat, 
               aes(x = time, y = exprs.adj, group = condition, linetype = genotype, shape = tissue)) + 
    geom_line() + geom_point(size = 3) + 
    theme_bw(24) + xlab("Time (ZT)") + 
    ylab("log2 expression") + 
    ggtitle(jtitle) + 
    theme(aspect.ratio = 0.5, legend.position = "bottom")
  if (add.vline){
    m1 <- m1 + geom_vline(xintercept = c(8, 20, 8+24, 20+24), linetype="dotted")
  }
  return(m1)
}

PlotGeneAcrossTissues <- function(dat, jtitle, convert.linear = FALSE){
  library(ggplot2)
  if (missing(jtitle)){
    jtitle = unique(dat$gene)
  }
  if (convert.linear){
    dat$exprs <- (2 ^ dat$exprs) - 1 
    jylab <- "mRNA expression (linear scale)"
  } else {
    jylab <- "mRNA expression (log2 scale)"
  }
  m <- ggplot(dat, aes(x = time, y = exprs,
                       group = experiment, 
                       colour = experiment)) + 
    geom_point() + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = jylab)
  return(m)
}

PlotGeneAcrossTissuesRnaseq <- function(dat, jtitle, convert.linear = FALSE){
  library(ggplot2)
  if (missing(jtitle)){
    jtitle = unique(dat$gene)
  }
  if (convert.linear){
    dat$exprs <- (2 ^ dat$exprs) - 1 
    jylab <- "mRNA expression (linear scale)"
  } else {
    jylab <- "mRNA expression (log2 scale)"
  }
  
  m <- ggplot(dat, aes(x = time, y = exprs)) + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = jylab) +
    xlab("CT") +
    scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12)) +
    theme(axis.text.x=element_text(angle=90,vjust = 0)) + theme_bw(24)
}

PlotTpmAcrossTissues <- function(dat, jtitle, log2.transform=FALSE){
  library(ggplot2)
  
  if (missing(jtitle)){
    jgene <- unique(dat$gene_name)
    jtranscript <- unique(dat$transcript_id)
    jtitle = paste(jgene, jtranscript)
  }
  if (log2.transform == FALSE){
    m <- ggplot(dat, aes(x = time, y = tpm, group = transcript_id, colour = transcript_id, shape = transcript_id)) 
    jylab <- "TPM expression"
  } else {
    m <- ggplot(dat, aes(x = time, y = log2(tpm + 0.01), group = transcript_id, colour = transcript_id, shape = transcript_id))
    jylab <- "log2 TPM expression"
  }
  m <- m + theme(legend.position="bottom") +
    geom_point() + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = jylab)
  return(m)
}

PlotRnaseqAcrossTissues <- function(dat, jtitle){
  p <- ggplot(dat, aes(x = time, y = exprs)) + 
    geom_line() + 
    facet_wrap(~tissue) +
    ggtitle(jtitle) + 
    ylab(label = "log2 mRNA expression") +
    xlab("CT") +
    scale_x_continuous(limits = c(18, 64), breaks = seq(24, 64, 12))  +
    theme_bw(24) + 
    theme(axis.text.x=element_text(angle=90,vjust = 0))
  return(p)
}

Center <- function(x){
  # Running scale on dplyr doesnt work, try it manually
  return(x - mean(x))
}

PlotGeneNormalized <- function(dat, jtitle){
  if (missing(jtitle)){
    jtitle <- unique(dat$gene)
  }
  dat.norm <- dat %>%
    group_by(tissue) %>%
    mutate(exprs.scaled = Center(exprs))
  p <- ggplot(dat.norm, aes(x = time, y = exprs.scaled, group = tissue, colour = tissue, fill = tissue)) + 
    geom_line() +
    geom_point() + 
    ylab("Centered mRNA expression") +
    xlab("CT")
  return(p)
}

PlotEncodeRnaseq <- function(dat, jtitle, sort.by.tissue = TRUE, by.var = "tpm"){
  source("~/projects/tissue-specificity/scripts/functions/SortByTissue.R")
  dat <- SortByTissue(dat, by.var = by.var)
  
  if (missing(jtitle)){
    jtitle <- unique(dat$gene)
  }
  p <- ggplot(dat, aes(x = tissue, y = tpm)) + geom_bar(stat = "identity") + ggtitle(jtitle) + xlab("") + ylab("Gene expression (TPM)") + 
    theme_bw(24) + theme(axis.text.x=element_text(angle=90,vjust = 0, hjust = 1),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())
  return(p)
}

CalculatePeriodogramLong <- function(dat, jexperiment = "array", remove.inf = TRUE){
  # interval between sampling points, in hours
  if (jexperiment == "array"){
    interval <- 2  # hrs
  } else {
    interval <- 6  # hrs
  }
  exprs <- subset(dat, experiment == jexperiment)$exprs
  p <- CalculatePeriodogram(exprs)
  periods <- signif(interval / p$freq, digits = 3)
  dat.var.s <- data.frame(periodogram = p$p.scaled, period = periods)
  dat.var.s$period <- factor(dat.var.s$period, 
                              levels = sort(unique(dat.var.s$period), decreasing = TRUE))
  if (remove.inf){
    dat.var.s <- subset(dat.var.s, period != Inf)
  }
  return(dat.var.s)
}

PlotPeriodogramLong <- function(dat, jexperiment = "array", jtitle = "Plot title"){
  # expect dat to be by gene
  # Plot periodogram of the gene
  source("~/projects/tissue-specificity/scripts/functions/FourierFunctions.R")
  dat.periodogram <- dat %>%
    group_by(gene, tissue) %>%
    do(CalculatePeriodogramLong(., jexperiment))
  ggplot(dat.periodogram, aes(x = period, y = periodogram)) + facet_wrap(~tissue) +  geom_bar(stat = "identity")
}