library(PhaseHSV)
library(scales)


# Around the Clock Plots --------------------------------------------------

bToKb <- function(){ 
  function(x) x / 1000 
}

scale1decimal <- function(x) sprintf("%.1f", x)
scale2decimal <- function(x) sprintf("%.2f", x)

MakeGenotypeLabeller <- function(counts.long.merged, posrange, fit.first.last.pos){
  # make genotype labeller to edit facet texts. Fits a sinusoid to data and outputs pval, phase, log2fc for each genotype
  jfits <- counts.long.merged %>%
    group_by(genotype, bait) %>%
    do(FormatAndFitRhythmAcrossFrags(., posrange = posrange, log10.to.log2 = TRUE, fit.first.last.pos = fit.first.last.pos))
  geno_lst <- paste(jfits$genotype, paste0("Pval=", signif(jfits$pval, 2)), paste0("Phase=", signif(jfits$phase, 2)), paste0("Log2FC=", signif(jfits$logfc, 2)))
  names(geno_lst) <- jfits$genotype
  return(geno_lst)
}

PlotPvalChiSqr <- function(counts.complex, do.facet=FALSE, pos.max=1e5, jposes=c(0), bsize=11){
  # Plot p-value significance for amplitude (chi-square test)
  jsub <- subset(counts.complex, abs(pos) < pos.max) %>% 
    arrange(pos)
  # plot line OVER point
  m2 <- ggplot(jsub, aes(group = tissue, x = pos, y = -log10(chisqr.pval))) + 
    geom_point(aes(color = hex.col.chisqr, shape = tissue), size = 2.6) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(linetype = tissue), colour = "gray40") +    # two colors
    geom_line(data = subset(jsub, LR == "Right"), aes(linetype = tissue), colour = "gray40") +    # two colors
    scale_color_identity() +  
    geom_vline(xintercept = jposes, linetype = "dotted") +
    theme_bw(bsize) + 
    theme(aspect.ratio = 0.25, legend.position = "bottom") +
    xlab("Position relative to bait [kb]") +
    ylab("-log10(Pvalue)") + 
    scale_x_continuous(labels=bToKb()) + 
    scale_shape_manual(values = c(16,3))
  if (do.facet){
    m2 <- m2 + facet_wrap(~tissue, ncol = 1) + theme(strip.background = element_blank(),strip.text.y = element_blank())
  }
  return(m2)
}

PlotAmpZscore4C <- function(counts.complex, do.facet=FALSE, pos.max=1e5, jposes = c(0)){
  jsub <- subset(counts.complex, abs(pos) < pos.max) %>% 
    arrange(pos)
  
  m2 <- ggplot(jsub, aes(group = tissue, x = pos, y = zscore)) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(linetype = tissue), colour = "gray40") +    # two colors
    geom_line(data = subset(jsub, LR == "Right"), aes(linetype = tissue), colour = "gray40") +    # two colors
    geom_point(aes(color = hex.col2, shape = tissue), size = 2.6) + 
    scale_color_identity() +  
    geom_vline(xintercept = jposes, linetype = "dotted") +
    theme_bw() + 
    theme(aspect.ratio = 0.25, legend.position = "bottom") +
    xlab("Position relative to bait [kb]") +
    ylab("Zscore") + 
    scale_x_continuous(labels=bToKb()) + 
    scale_shape_manual(values = c(16,3))
  if (do.facet){
   m2 <- m2 + facet_wrap(~tissue, ncol = 1) + theme(strip.background = element_blank(),strip.text.y = element_blank())
  }
  return(m2)
}

PlotAmps4C <- function(counts.complex, do.facet=FALSE, pos.max=1e5, jposes = c(0), jshapes = c(16, 3), show.ribbon = FALSE, .log2 = FALSE, bsize = 11){
  # counts.complex from Project4CtoZscore function
  jsub <- subset(counts.complex %>% arrange(pos), abs(pos) < pos.max)
  if (.log2){
    jsub$amp <- jsub$amp / log10(2)
    jsub$amp.se <- jsub$amp.se / log10(2)
    jylab <- "Log2 Fold Change"
  } else {
    eps <- 0.5  # arbitrary
    jsub$amp <- 2^(jsub$amp + eps)
    jsub$amp.se <- 2^(jsub$amp.se + eps)
    jylab <- "Fold Change"
  }
  m1 <- ggplot(jsub, aes(x = pos, y = amp, group = tissue)) + 
    scale_color_identity() + 
    geom_vline(xintercept = jposes, linetype = "dotted") +
    theme_bw(bsize) + 
    theme(aspect.ratio = 0.25, legend.position = "bottom") +
    xlab("Position relative to bait [kb]") +
    ylab(jylab) + 
    scale_shape_manual(values = jshapes) + 
    scale_x_continuous(labels=bToKb())
  # draw lines OVER points
  m1 <- m1 + 
    geom_point(aes(colour = hex.col2, shape = tissue), size=2.6)
  if (show.ribbon){
    m1 <- m1 + 
      geom_ribbon(data = subset(jsub, LR == "Left"), mapping = aes(ymin = amp - amp.se, ymax = amp + amp.se), fill = "gray95") + 
      geom_ribbon(data = subset(jsub, LR == "Right"), mapping = aes(ymin = amp - amp.se, ymax = amp + amp.se), fill = "gray95")
  } else {
    m1 <- m1 + 
      geom_line(data = subset(jsub, LR == "Left"), aes(linetype = tissue), colour = "gray40") +    # two colors
      geom_line(data = subset(jsub, LR == "Right"), aes(linetype = tissue), colour = "gray40")    # two colors
  }
  if (do.facet){
    m1 <- m1 + facet_wrap(~tissue, ncol = 1) + theme(strip.background = element_blank(),strip.text.y = element_blank())
  }
  return(m1)
}

PlotSignalOverTime <- function(counts.long.merged, posrange, do.facet=TRUE, show.first.last.pos=FALSE, show.best=FALSE, best.crit = "pval"){
  
  jtitle <- paste0("Fragment range: ", paste(posrange, collapse = " to "))
  
  counts.long.merged$time <- as.character(counts.long.merged$time)
  jsub <- subset(counts.long.merged, pos > posrange[[1]] & pos < posrange[[2]])
  if (show.first.last.pos){
    jsub <- subset(counts.long.merged, pos %in% range(jsub$pos))
  }
  if (all(startsWith(jsub$time, prefix = "ZT"))){
    jsub$time <- as.numeric(gsub("ZT", "", jsub$time))
  }
  if (show.best){
    # return only a single fragment, choosing one with highest amplitude
    jsub.form <- FormatDat(jsub, posrange, log10.to.log2 = TRUE)
    jfits <- jsub.form %>%
      group_by(genotype, bait, pos) %>%
      do(FitRhythmAcrossFrags(., T.period=24, region.name=paste(posrange, collapse="to"), get.residuals = TRUE, signal.cname="A", label=NULL, window.cname = "pos", time.cname = "time", single.frag=TRUE)) %>%
      arrange(pval)
    print(paste("Selecting best fragment by:", best.crit))
    if (best.crit == "pval"){
      # get most rhythmic fragment by pval
      jfits.sub <- subset(jfits, pval == min(pval))
    } else if (best.crit == "amp"){
      # get most rhythmic fragment by amp
      jfits.sub <- subset(jfits, amp == max(amp))
    }
    pos.best <- jfits.sub$pos[[1]]
    jsub <- subset(jsub, pos == pos.best)
    jstats <- paste0(paste(c("Pval=", "Phase=", "Amp="), signif(c(jfits.sub$pval, jfits.sub$phase, jfits.sub$amp), 2), sep = ""), collapse = ",")
    jtitle <- paste0(jtitle, "\nBest frag:", jfits.sub$pos, " ", jstats)
  }
  m.osc <- ggplot(jsub, aes(x = time, y = A, group = interaction(as.factor(pos), genotype), linetype = genotype)) + 
    geom_point() + 
    geom_line() + 
    # geom_errorbar(mapping = aes(ymin = A - Avar, ymax = A + Avar), width = 0, colour = "gray50", linetype = "solid") + 
    geom_errorbar(mapping = aes(ymin = A - sqrt(Avar * sig2.adj), ymax = A + sqrt(Avar * sig2.adj)), width = 0, colour = "gray50", linetype = "solid") + 
    theme_bw() + 
    ggtitle(jtitle) + 
    theme(legend.position = "bottom", aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("Time [ZT]") + ylab("Log10 Signal")  +
    scale_x_continuous(breaks = seq(0, 20, by = 4))
  if (do.facet){
    m.osc <- m.osc + facet_wrap(~genotype, ncol = 2)
  }
  return(m.osc)
}


# Delta Time Plots --------------------------------------------------------



Log10PosNeg <- Vectorize(function(x){
  # Do log but handle negatives as if it is positive
  if (x > 0){
    return(log10(x))
  } else {
    return(-log10(-x))
  }
}, vectorize.args = "x")

asinh_trans <- function(){
  scales::trans_new(name = 'asinh', transform = function(x) asinh(x), 
                    inverse = function(x) sinh(x))
}

scientific_10_posexpo <- function(x) {
  # http://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
  parse(text=gsub("1e\\+", " 10^", scientific_format()(x)))
}

AnnotateLeftRightTrans <- function(counts){
  counts.cis <- counts[!grepl("^trans", as.character(counts$bin)), ]
  counts.cis$LR <- sapply(counts.cis$pos, AddLR)
  counts.trans <- counts[grepl("^trans", as.character(counts$bin)), ]
  counts.trans$LR <- "trans"
  counts <- rbind(counts.cis, counts.trans)
  return(counts)
}


# define the summary function
f <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

PlotSignalWithHeatmap <- function(counts.delt, counts.long, jbait, jgeno, jtiss, maxdist, binsize=1000, zscorelimits=c(-10, 10)){
  
  # visualize by heatmap
  jsub <- subset(counts.delt, bait == jbait & tissue == jtiss & genotype == jgeno & abs(pos) < maxdist)
  jsubA <- subset(counts.long, bait == jbait & tissue == jtiss & genotype == jgeno & abs(pos) < maxdist)
  
  # do interpolation for heatmap
  x <- seq(-maxdist, maxdist, by = binsize)
  y <- approx(x = jsub$pos, y = -jsub$zscore.row, xout = x)
  jsub.int <- data.frame(pos = x, zscore = y, genotype = jgeno)
  
  plots.signal <- PlotSignalLR(jbait, jsubA, 500, maxdist, jtitle = paste(jbait, jgeno), do.facet = FALSE)
  plots.heatmap <- ggplot(jsub.int, aes(x=pos, y=genotype, fill=zscore.y)) +
    geom_tile() +
    labs(x="",y="") +
    scale_fill_gradientn(name="Zscore", limits=zscorelimits,
                         colours=colorRampPalette(c('red', 'black', 'blue'))(20)) +
    theme(legend.position="bottom", 
          legend.key.width=unit(.1,"npc"),legend.key.height=unit(.05,"npc"),
          axis.ticks=element_blank(), aspect.ratio = 0.125 / 2) + 
    # ggtitle(paste("Zscore", jbait)) + 
    scale_x_continuous(expand = c(0, 0)) +
    geom_vline(xintercept = 0, colour = "white", linetype = "dotted")
  multiplot(plots.signal, plots.heatmap)
}

OrderDecreasing <- function(dat, jfactor, jval){
  # Reorder factors by decreasing value of jval
  # used in conjunction with barplots makes life easier
  dat[[jfactor]] <- factor(dat[[jfactor]], levels = dat[[jfactor]][order(dat[[jval]], decreasing = TRUE)])
  return(dat)
}

PlotVolcano <- function(jbait, counts.delt.sub, jxlab = "", jxlim = c(-0.4, 0.4), by.tissue = TRUE, jtitle = NULL){
  if (is.null(jtitle)){
    jtitle <- jbait
  }
  vol <- ggplot(subset(counts.delt.sub, bait == jbait), aes(x = A.delta, y = -log10(pval.row))) + 
    geom_point(alpha = 0.5) + theme(aspect.ratio = 1) + xlab(jxlab) + xlim(jxlim)
  if (by.tissue){
    vol <- vol + facet_grid(genotype ~ tissue)
  } else {
    vol <- vol + facet_grid(genotype ~ time)
  }
  return(vol)
}

PlotSignal <- function(jbait, counts.long, pseudo.low, mindist, jtitle=NULL){
  if (is.null(jtitle)){
    jtitle <- paste(jbait, "Pseudo:", pseudo.low)
  }
  AvsPosT <- ggplot(subset(counts.long, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait), aes(x = pos, y = A, colour = time)) + 
    geom_line() + facet_grid(genotype ~ tissue) + theme(aspect.ratio = 0.25) + geom_hline(aes(yintercept=0)) + 
    theme_bw() + ggtitle(jtitle) + geom_vline(aes(xintercept=0), linetype="dotted")
  return(AvsPosT)
}

PlotSignalLR <- function(jbait, counts.long, pseudo.low, mindist, jtitle=NULL, do.facet=TRUE, ltype="solid"){
  if (is.null(jtitle)){
    jtitle <- paste(jbait, "Pseudo:", pseudo.low)
  }
  jsub <- subset(counts.long, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  # zts <- c(8, 20) - 8  # rotate
  # phase.palette <- hsv(PhaseToHsv(zts, min.phase = 0, max.phase = 24), 1, 1)
  AvsPosT <- ggplot(jsub, aes(x = pos, y = A, colour = time)) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = A, colour = time), linetype = ltype) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = A, colour = time), linetype = ltype) + 
    theme(aspect.ratio = 0.25) + geom_hline(aes(yintercept=0)) + 
    theme_bw() + ggtitle(jtitle) + geom_vline(aes(xintercept=0), linetype="dotted") + 
    theme(legend.position = "bottom") + 
    xlab("Position relative to bait") + 
    ylab("Log10 Signal")
    # scale_colour_manual(values=phase.palette)
  if (do.facet){
    AvsPosT <- AvsPosT + facet_grid(genotype ~ tissue)
  }
  return(AvsPosT)
}

PlotSignalLR.geno <- function(jbait, counts.long, pseudo.low, mindist, jtitle=NULL, do.facet=TRUE, add.se=TRUE, jgray="gray75", add.vert.line=TRUE){
  jrotate <- -8  # rotate by 8 hrs
  # source("/home/yeung/projects/Rfunctions/RotatePhase.R")
  # source("/home/yeung/projects/Rfunctions/")
  source("scripts/functions/TrackHubFunctions.R")
  
  if (is.null(jtitle)){
    jtitle <- paste(jbait, "Pseudo:", pseudo.low)
  }
  jsub <- subset(counts.long, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  assertthat::assert_that(nrow(jsub) > 0)
  # zts <- c(8, 20) - 8  # rotate
  # phase.palette <- hsv(PhaseToHsv(zts, min.phase = 0, max.phase = 24), 1, 1)
  jtime <- jsub$time
  # jtime <- gsub("ZT12", "ZT012.7", jtime)  # yellow is hard to see on white
  
  zts <- RotatePhase(as.numeric(gsub("ZT", "", jtime)), rotate.hr = jrotate)
  jsub$hsv.col <- hsv(PhaseToHsv(zts, min.phase = 0, max.phase = 24), 1, 1)
  
  # reorder time
  jsub$time <- factor(as.character(jsub$time), levels = paste("ZT", c("04", "08", "12", "16", "20", "00"), sep = ""))
  
  AvsPosT <- ggplot(jsub, aes(x = pos, y = A, linetype = genotype, colour = hsv.col))  # init
  if (add.se){
    AvsPosT <- AvsPosT + 
        geom_ribbon(data = subset(jsub, LR == "Left"), mapping = aes(ymin = A - sqrt(Avar * sig2.adj), ymax = A + sqrt(Avar * sig2.adj)), fill = jgray, colour = jgray) + 
        geom_ribbon(data = subset(jsub, LR == "Right"), mapping = aes(ymin = A - sqrt(Avar * sig2.adj), ymax = A + sqrt(Avar * sig2.adj)), fill = jgray, colour = jgray)
  }
  AvsPosT <- AvsPosT + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = A, linetype = genotype, colour = hsv.col)) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = A, linetype = genotype, colour = hsv.col)) + 
    theme(aspect.ratio = 0.25) + geom_hline(aes(yintercept=0)) + 
    theme_bw() + ggtitle(jtitle) + 
    theme(legend.position = "bottom") + 
    xlab("Position relative to bait") + 
    ylab("Log10 Signal") + 
    scale_color_identity() +
    scale_linetype_manual(values = c("solid", "dashed"))
  # scale_colour_manual(values=phase.palette)
  if (do.facet){
    AvsPosT <- AvsPosT + facet_grid(genotype ~ tissue)
  }
  if (add.vert.line){
    AvsPosT <- AvsPosT + geom_vline(aes(xintercept=0), linetype="dotted")
  }
  return(AvsPosT)
}

PlotPvalLR.geno <- function(jbait, counts.long, pseudo.low, mindist, jtitle=NULL, do.facet=TRUE){
  
  if (is.null(jtitle)){
    jtitle <- paste(jbait, "Pseudo:", pseudo.low)
  }
  jsub <- subset(counts.long, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  # make long
  pstat.cols.i <- grep("Pstat\\.ZT", colnames(jsub))
  jsub.long <- melt(jsub, id.vars = c("pos", "LR"), measure.vars = colnames(jsub)[c(pstat.cols.i)], variable.name = "timeZT", value.name = "pval")
  jsub.long <- jsub.long[!duplicated(jsub.long), ]
  jsub.long$timeZT <- gsub("Pstat\\.", "", jsub.long$timeZT)
  
  AvsPosT <- ggplot(jsub.long, aes(x = pos, y = -log10(pval))) + 
    geom_line(data = subset(jsub.long, LR == "Left"), aes(x = pos, y = -log10(pval))) + 
    geom_line(data = subset(jsub.long, LR == "Right"), aes(x = pos, y = -log10(pval))) + 
    theme(aspect.ratio = 0.25) + geom_hline(aes(yintercept=0)) + 
    theme_bw() + ggtitle(jtitle) + geom_vline(aes(xintercept=0), linetype="dotted") + 
    theme(legend.position = "bottom", aspect.ratio = 1, strip.background = element_blank(),strip.text.y = element_blank()) + 
    scale_x_continuous(labels=bToKb()) + 
    xlab("Position relative to bait") + 
    ylab("-log10(P-value)") + 
    scale_color_identity() 
  # scale_colour_manual(values=phase.palette)
  if (do.facet){
    AvsPosT <- AvsPosT + facet_wrap(~timeZT)
  }
}

PlotYOverTime <- function(jbait, counts.long, pseudo.low, mindist, jtitle=NULL, do.facet=TRUE, cname="Zscore.time"){
  if (is.null(jtitle)){
    jtitle <- paste(jbait, "Pseudo:", pseudo.low)
  }
  jsub <- subset(counts.long, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  # zts <- c(8, 20) - 8  # rotate
  # phase.palette <- hsv(PhaseToHsv(zts, min.phase = 0, max.phase = 24), 1, 1)
  AvsPosT <- ggplot(jsub, aes_string(x = "pos", y = cname, colour = "time")) + 
    geom_line(data = subset(jsub, LR == "Left"), aes_string(x = "pos", y = cname, colour = "time")) + 
    geom_line(data = subset(jsub, LR == "Right"), aes_string(x = "pos", y = cname, colour = "time")) + 
    theme(aspect.ratio = 0.25) + geom_hline(aes(yintercept=0)) + 
    theme_bw() + ggtitle(jtitle) + geom_vline(aes(xintercept=0), linetype="dotted") + 
    theme(legend.position = "bottom") + 
    xlab("Position relative to bait") + 
    ylab("Log10 Signal")
  # scale_colour_manual(values=phase.palette)
  if (do.facet){
    AvsPosT <- AvsPosT + facet_grid(genotype ~ tissue)
  }
  return(AvsPosT)
}

PlotSignalLR.WTKOMerged <- function(jbait, counts.long, pseudo.low, mindist, jtitle=NULL, do.facet=TRUE, jalphas = c(1, 0.5)){
  # works with genotype as dotted
  if (is.null(jtitle)){
    jtitle <- paste(jbait, "Pseudo:", pseudo.low)
  }
  jsub <- subset(counts.long, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  # jsub$genotype <- factor(as.character(jsub$genotype), levels = gene.order)
  # zts <- c(8, 20) - 8  # rotate
  # phase.palette <- hsv(PhaseToHsv(zts, min.phase = 0, max.phase = 24), 1, 1)
  AvsPosT <- ggplot(jsub, aes(x = pos, y = A, colour = time, linetype = genotype, alpha = genotype)) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = A, colour = time)) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = A, colour = time)) + 
    geom_hline(aes(yintercept=0)) + 
    theme_bw() + ggtitle(jtitle) + geom_vline(aes(xintercept=0), linetype="dotted") + 
    theme(legend.position = "bottom") + 
    xlab("Position relative to bait") + 
    ylab("Log10 Signal") + 
    scale_alpha_manual(values = jalphas) +
    theme(aspect.ratio = 0.25)
  # scale_colour_manual(values=phase.palette)
  if (do.facet){
    AvsPosT <- AvsPosT + facet_grid(genotype ~ tissue)
  }
  return(AvsPosT)
}

PlotDelta <- function(jbait, counts.delt, pseudo.low, mindist, jylim = c(-0.2, 0.2), jtitle=NULL){
  if (is.null(jtitle)){
    jtitle <- paste(jbait, "Pseudo:", pseudo.low)
  }
  AdeltvsPosBaitT <- ggplot(subset(counts.delt, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait), aes(x = pos, y = A.delta)) + 
    geom_line() + facet_grid(genotype ~ tissue) + theme(aspect.ratio = 0.25) + geom_vline(aes(xintercept=0), linetype="dotted") + geom_hline(aes(yintercept=0)) + 
    ggtitle(jtitle) + theme_bw() + geom_vline(aes(xintercept=0), linetype="dotted") +
    ylim(jylim)
  return(AdeltvsPosBaitT)
}

PlotDeltaLR <- function(jbait, counts.delt, pseudo.low, mindist, jylim = c(-0.2, 0.2), jtitle=NULL, add.peaks=FALSE){
  if (is.null(jtitle)){
    jtitle <- paste(jbait, "Pseudo:", pseudo.low)
  }
  if (add.peaks != FALSE){
    jsub.peaks <- subset(add.peaks, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  }
  jsub <- subset(counts.delt, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  AdeltvsPosBaitT <- ggplot(jsub, aes(x = pos, y = A.delta)) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = A.delta)) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = A.delta)) + 
    facet_grid(genotype ~ tissue) + theme(aspect.ratio = 0.25) + geom_vline(aes(xintercept=0), linetype="dotted") + geom_hline(aes(yintercept=0)) + 
    ggtitle(jtitle) + theme_bw() + geom_vline(aes(xintercept=0), linetype="dotted") +
    xlab("Position relative to bait") + 
    ylab("Difference in Log10 signal (ZT20 - ZT08)")
  if (add.peaks != FALSE){
    # AdeltvsPosBaitT <-  Adelt+ geom_point(data = subset(counts.delt.sub, bait == jbait), aes(x = pos, y = A.delta))
    AdeltvsPosBaitT <- AdeltvsPosBaitT + geom_point(data = jsub.peaks, aes(x = pos, y = A.delta), colour = "blue")
  }
  if (!is.null(jylim)){
    AdeltvsPosBaitT <- AdeltvsPosBaitT + ylim(jylim)
  }
  return(AdeltvsPosBaitT)
}

PlotPval <- function(jbait, counts.long, pseudo.low, mindist){
  Pvals <- ggplot(subset(counts.delt, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait), aes(x = pos, y = -log10(pval.row))) + 
    geom_line() + facet_grid(genotype ~ tissue) + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0), linetype="dotted") + theme(aspect.ratio = 0.25) +
    theme_bw() + 
    ggtitle(paste(jbait, "Pseudo:", pseudo.low))
  return(Pvals)
}

PlotSignalDelta <- function(jbait, counts.long, counts.delt, pseudo.low, mindist, jylim = c(-0.2, 0.2), jtitle=NULL){
  AvsPosT <- PlotSignal(jbait, counts.long, pseudo.low, mindist, jtitle)
  AdeltvsPosBaitT <- PlotDelta(jbait, counts.delt, pseudo.low, mindist, jylim = jylim, jtitle = "") 
  multiplot(AvsPosT, AdeltvsPosBaitT) 
}

PlotSignalDeltaLR <- function(jbait, counts.long, counts.delt, pseudo.low, mindist, jylim = c(-0.2, 0.2), jtitle=NULL, add.peaks=FALSE){
  AvsPosT <- PlotSignalLR(jbait, counts.long, pseudo.low, mindist, jtitle)
  AdeltvsPosBaitT <- PlotDeltaLR(jbait, counts.delt, pseudo.low, mindist, jylim = jylim, jtitle = "", add.peaks) 
  multiplot(AvsPosT, AdeltvsPosBaitT) 
}

PlotSignalLivVsKid <- function(jbait, counts.long, pseudo.low, mindist){
  jsub <- subset(counts.long, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  AvsPosT <- ggplot(jsub, aes(x = pos, y = A, colour = tissue)) + 
    geom_line() + facet_grid(genotype ~ time) + theme(aspect.ratio = 0.25) + geom_hline(aes(yintercept=0)) + theme_bw() + ggtitle(paste(jbait, "Pseudo:", pseudo.low)) + geom_vline(aes(xintercept=0), linetype="dotted")
  return(AvsPosT)
}

PlotSignalLivVsKidLR <- function(jbait, counts.long, pseudo.low, mindist, jtitle = "", do.facet = TRUE, show.legend = TRUE, in.kb = FALSE){
  jsub <- subset(counts.long, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  if (in.kb){
    jxlab <- "Position Relative to Bait (kb)"
    jsub$pos <- jsub$pos / 1000
  } else {
    jxlab <- "Position Relative to Bait (kb)"
  }
  AvsPosT <- ggplot(jsub, aes(x = pos, y = A, colour = tissue)) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = A)) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = A)) + 
    theme_bw() + 
    theme(aspect.ratio = 0.25, legend.position = "bottom") + geom_hline(aes(yintercept=0)) +  ggtitle(jtitle) + geom_vline(aes(xintercept=0), linetype="dotted") + 
    xlab(jxlab) + ylab("Normalized Reads (Log10)")
  if (do.facet){
    AvsPosT <- AvsPosT + facet_grid(genotype ~ time)
  } 
  if (!show.legend){
    AvsPosT <- AvsPosT + theme(legend.position = "none")
  }
  return(AvsPosT)
}

PlotPvalLK <- function(jbait, counts.delt.lk, pseudo.low, mindist, jtitle = "", jylim = c(0, 10)){
  jdat <- subset(counts.delt.lk, bait == jbait & pseudo == pseudo.low & abs(pos) < mindist)
  PvalsLK <- ggplot(jdat, aes(x = pos, y = -log10(pval.row))) + 
    geom_line() + facet_grid(genotype~time) + theme(aspect.ratio = 0.25) +
    theme_bw() + 
    ggtitle(jtitle) + geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept=0), linetype="dotted") + 
    ylim(jylim) + 
    xlab("Position relative to bait") + 
    ylab("-log10(P-value)")
  return(PvalsLK)
} 

PlotPvalLKLR <- function(jbait, counts.delt.lk, pseudo.low, mindist, add.peaks=FALSE){
  jsub <- subset(counts.delt.lk, bait == jbait & pseudo == pseudo.low & abs(pos) < mindist)
  if (add.peaks != FALSE){
    jsub.peaks <- subset(add.peaks, bait == jbait & pseudo == pseudo.low & abs(pos) < mindist)
  }
  PvalsLK <- ggplot(jsub, aes(x = pos, y = -log10(pval.row))) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = -log10(pval.row))) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = -log10(pval.row))) + 
    facet_grid(genotype~time) + theme(aspect.ratio = 0.25) +
    theme_bw() + 
    ggtitle(paste(jbait)) + geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept=0), linetype="dotted")
  if (add.peaks != FALSE){
    PvalsLK <- PvalsLK + geom_point(data = jsub.peaks, aes(x = pos, y = -log10(pval.row)), colour = "blue")
  }
  return(PvalsLK)
} 

PlotDeltaLKLR <- function(jbait, counts.delt.lk, pseudo.low, mindist, add.peaks=FALSE){
  jsub <- subset(counts.delt.lk, bait == jbait & pseudo == pseudo.low & abs(pos) < mindist)
  if (add.peaks != FALSE){
    jsub.peaks <- subset(add.peaks, bait == jbait & pseudo == pseudo.low & abs(pos) < mindist)
  }
  DeltaLK <- ggplot(jsub, aes(x = pos, y = -log10(pval.row))) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = A.delta)) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = A.delta)) + 
    facet_grid(genotype~time) + theme(aspect.ratio = 0.25) +
    theme_bw() + 
    ggtitle(paste(jbait)) + geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept=0), linetype="dotted")
  if (add.peaks != FALSE){
    DeltaLK <- DeltaLK + geom_point(data = jsub.peaks, aes(x = pos, y = A.delta), colour = "blue")
  }
  return(DeltaLK)
} 

PlotZscorePval <- function(jbait, counts.delt, pseudo.low, mindist, jylim.zscore = c(-5.5, 5.5), jylim.pval = c(0, 10)){
  Zscores <- ggplot(subset(counts.delt, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait), aes(x = pos, y = -zscore.row)) + 
    geom_line() + facet_grid(genotype ~ tissue) + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0), linetype="dotted") + theme(aspect.ratio = 0.25) +
    theme_bw() + 
    ggtitle(paste(jbait, "Pseudo:", pseudo.low)) + ylim(jylim.zscore)
  Pvals <- ggplot(subset(counts.delt, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait), aes(x = pos, y = -log10(pval.row))) + 
    geom_line() + facet_grid(genotype ~ tissue) + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0), linetype="dotted") + theme(aspect.ratio = 0.25) +
    theme_bw() + 
    ggtitle(paste(jbait, "Pseudo:", pseudo.low)) + ylim(jylim.pval)
  multiplot(Zscores, Pvals)
}

PlotPvalLR <- function(jbait, counts.delt, pseudo.low, mindist, jylim, jtitle){
  jsub <- subset(counts.delt, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  Pvals <- ggplot(jsub, aes(x = pos, y = -log10(pval.row))) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = -log10(pval.row))) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = -log10(pval.row))) + 
    facet_grid(genotype ~ tissue) + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0), linetype="dotted") + theme(aspect.ratio = 0.25) +
    theme_bw() + 
    ggtitle(jtitle) + 
    xlab("Position relative to bait") + 
    ylab("-log10(P-value)")
  if (!is.null(jylim)){
    Pvals <- Pvals + ylim(jylim)
  }
  return(Pvals)
}



PlotZscorePvalLR <- function(jbait, counts.delt, pseudo.low, mindist, jylim.zscore = c(-5.5, 5.5), jylim.pval = c(0, 10)){
  jsub <- subset(counts.delt, abs(pos) < mindist & pseudo == pseudo.low & bait == jbait)
  Zscores <- ggplot(jsub, aes(x = pos, y = -zscore.row)) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = -zscore.row)) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = -zscore.row)) + 
    facet_grid(genotype ~ tissue) + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0), linetype="dotted") + theme(aspect.ratio = 0.25) +
    theme_bw() + 
    ggtitle(paste(jbait, "Pseudo:", pseudo.low)) + ylim(jylim.zscore)
  Pvals <- ggplot(jsub, aes(x = pos, y = -log10(pval.row))) + 
    geom_line(data = subset(jsub, LR == "Left"), aes(x = pos, y = -log10(pval.row))) + 
    geom_line(data = subset(jsub, LR == "Right"), aes(x = pos, y = -log10(pval.row))) + 
    facet_grid(genotype ~ tissue) + geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0), linetype="dotted") + theme(aspect.ratio = 0.25) +
    theme_bw() + 
    ggtitle(paste(jbait, "Pseudo:", pseudo.low)) + ylim(jylim.pval)
  multiplot(Zscores, Pvals)
}

SetYXMaxDelta <- function(counts.delt, cname="pval"){
  # calculate xmax and ymax for plotting volcano plots
  xmax <- max(abs(counts.delt$A.delta))
  ymax <- max(abs(-log10(counts.delt[[cname]][which(!is.na(counts.delt[[cname]]))])))
  return(list(xmax=xmax, ymax=ymax))
}

VolcanoPlot <- function(counts.long, adj.pval.cutoff = 0.1){
  # counts.long: object in long format after running model
  # adj.pval.cutoff: after filtering for peaks what we consider to be pval cutoff
  # Return volcano plot 
  # get deltas
  counts.delt.sub <- GetDeltSub(counts.long, adj.pval.cutoff)
  xmax <- max(abs(counts.delt.sub$A.delta))
  ymax <- max(abs(-log10(counts.delt.sub$pval)))
  volcano.plot <- ggplot(counts.delt.sub, aes(x = A.delta, y = -log10(pval), colour = is.sig)) + 
    geom_point(alpha = 0.5) + ylim(0, ymax) + 
    xlim(-xmax, xmax) + facet_wrap(~bait) + theme_bw(24) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))  
  return(volcano.plot)
}

my_line <- function(x,y,...){
  # http://stackoverflow.com/questions/17793690/adding-line-of-identity-to-correlation-plots-using-pairs-command-in-r
  points(x,y, pch=".", ...)
  #text(x, y, labels = dots$jlabels, cex = 1.1, col = "red")
  abline(a = 0, b = 1)
}

scatter_with_xy_line <- function(x, y, ...){
  # if there are horizontal and vertical lines plot those too with options v and h
  vline_hline_coords <- list(...)
  points(x, y, pch = ".")
  abline(a = 0, b = 1)
  # add hline vline
  abline(h = vline_hline_coords$h)
  abline(v = vline_hline_coords$v)
}

my_line_custom <- function(x,y,...){
  # http://stackoverflow.com/questions/17793690/adding-line-of-identity-to-correlation-plots-using-pairs-command-in-r
  dots <- list(...)
  points(x,y, col = dots$col, cex = dots$cex, pch = dots$pch)
  text(x, y, labels = dots$jlabels, cex = 1.1, col = "red")
  abline(a = 0, b = 1)
}

my_line_nolabels <- function(x,y,...){
  points(x,y, pch = ".", ...)
  abline(a = 0, b = 1)
}

panel.hist <- function(x, ...)
{
  # remove minimum value (optional)
  x <- x[which(x > min(x))]
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  # text(0.5, 0.5, txt, cex = cex * abs(r))
  text(0.5, 0.5, txt, cex = 1.5)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # stolen from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  # 
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

PlotTriple <- function(f1, f2, f3, n.boxes = 8){
  require(gridExtra)
  # source("scripts/functions/MergeCounts.R")
  
  jlay <- matrix(c(rep(1, n.boxes - 2), 2, 3), n.boxes, 1, byrow = TRUE)
  n <- table(jlay)
  
  gb1 <- ggplot_build(f1 + theme(legend.position = "none"))
  gb2 <- ggplot_build(f2)
  gb3 <- ggplot_build(f3)
  
  gA <- ggplot_gtable(gb1)
  gB <- ggplot_gtable(gb2)
  gC <- ggplot_gtable(gb3)
  
  g <- rbind(gA, gB, gC)
  
  # locate the panels in the gtable layout
  panels <- g$layout$t[grepl("panel", g$layout$name)]
  # assign new (relative) heights to the panels, based on the number of breaks
  # print(panels)
  # g$heights[panels] <- list(unit(n[[1]],"null"), unit(n[[2]], "null"), unit(n[[3]], "null"))
  
  grid.newpage()
  grid.draw(g)
}

PlotQuad <- function(f1, f2, f3, f4, n.boxes = 8){
  require(gridExtra)
  # source("scripts/functions/MergeCounts.R")  # load this outside of the function so it isnt hardcoded
  
  jlay <- matrix(c(rep(c(1, 2), each = (n.boxes - 2)), 3, 4), n.boxes * 2 - 2, 1, byrow = TRUE)
  n <- table(jlay)
  
  gb1 <- ggplot_build(f1)
  gb2 <- ggplot_build(f2)
  gb3 <- ggplot_build(f3)
  gb4 <- ggplot_build(f4)
  
  gA <- ggplot_gtable(gb1)
  gB <- ggplot_gtable(gb2)
  gC <- ggplot_gtable(gb3)
  gD <- ggplot_gtable(gb4)
  
  g <- rbind(gA, gB, gC, gD)
  
  # locate the panels in the gtable layout
  panels <- g$layout$t[grepl("panel", g$layout$name)]
  # assign new (relative) heights to the panels, based on the number of breaks
  # g$heights[panels] <- list(unit(n[[1]],"null"), unit(n[[2]], "null"), unit(n[[3]], "null"))
  
  grid.newpage()
  grid.draw(g)
}

