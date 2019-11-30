# Jake Yeung
# Date of Creation: 2018-10-28
# File: ~/projects/4c_seq/scripts/primetime_animations/animate_around_the_clock.R
# Animate batch 3 around the clock
# Stolen from: ~/projects/4c_seq/scripts/batch3_4cseq_analysis/load_robj_plot_summary.batch3.WT_vs_Cry1intronKO.R


rm(list=ls())
jstart <- Sys.time()

library(dplyr)
library(ggplot2)
library(PhaseHSV)
library(parallel)

setwd("~/projects/4c_seq")
source("scripts/functions/FitRhythmic4CSmoothed.R")
source("scripts/functions/LoadRobjPlotSummary.R")
source("scripts/functions/PlotFunctions.R")
source("scripts/functions/LoadNormalizationOutputs.R")
source("scripts/functions/TrackHubFunctions.R")
source("scripts/functions/FindPeaks.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/FitRhythmic.R")

library(PhaseHSV)
source("/home/yeung/projects/tissue-specificity/scripts/functions/SvdFunctions.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/CosSineFunctions.R")
source("/home/yeung/projects/tissue-specificity/scripts/functions/PhaseColorFunctions.R")


# Get colors --------------------------------------------------------------

zts <- seq(0, 20, 4)
rotate.hrs <- -8
zts.rotated <- RotatePhase(zts, rotate.hr = rotate.hrs)

hex.cols <- hsv(PhaseToHsv(zts.rotated, 0, 24), 1, 1)

rgb.cols <- HexColorToRgb(hex.cols)

jposes <- c(0)


# Inits -------------------------------------------------------------------

bname <- paste0("2017-09-25-HalfMegabaseFilter-ignore_5_demultiplex_bugfix_mm9_Cry1intronFix")

jbaits <- c("Hoxd4", "Cry1", "Cry1_TSS", "Cry1_UP")
posranges <- list(c(-15000, 15000), c(-29000, -19000), c(-31000, -22200), c(-36000, -29000))
jpval.k <- 0.01  # make more stringent?
jpval.k.chisqr <- 0.01  # colors points for Cry1_UP

jaspect.ratio <- 0.25

lst <- list()
for (i in seq(length(jbaits))){
  lst[[i]] <- list("jbait"=jbaits[[i]], "posrange"=posranges[[i]])
}


indx <- 3
# for (jbait in jbaits){
  
  sublst <- lst[[3]]
  
  
  jbait <- sublst[["jbait"]]
  posrange <- sublst[["posrange"]]
  
  inf <- paste0("/home/yeung/data/4c_seq/", bname, "/", jbait, ".normmax.Inf.log.FALSE.nfrags.5.remove.auto.sig.2500.Robj")
  load(inf, v=T)
  
  jsig <- 2500
  
  if (exists("counts.long")){
    print("Saving counts.long as counts.long.merged")
    counts.long.merged <- counts.long; rm(counts.long)
  }
  
  if (length(jbait) != 1){
    warning(paste("Should be one bait", jbait))
  }
  
  counts.long.merged$LR <- sapply(counts.long.merged$pos, AddLR)
  
  counts.long.merged$genotype <- factor(as.character(counts.long.merged$genotype), levels = c("Liver_WT", "Liver_Cry1intronKO"))
  
  # Fit rhythmic on every position ------------------------------------------
  
  counts.complex <- Project4CtoZscore(subset(counts.long.merged, sig == jsig), jpval.k = jpval.k, jpval.k.chisqr = jpval.k.chisqr)
  counts.complex$tissue <- factor(as.character(counts.complex$tissue), levels = c("Liver_WT", "Liver_Cry1intronKO"))
  
  # Look at Zscore comparisons ----------------------------------------------
  
  counts.long.zscore <- subset(counts.long.merged, time == "ZT20" & sig == jsig)
  counts.long.zscore$A <- counts.long.zscore$Zscore.ZT20
  
  # Find peaks in signal  ---------------------------------------------------
  
  counts.avg <- counts.long.merged %>%
    group_by(bait, genotype, pos, sig, pseudo, LR) %>%
    summarise(A.se = sd(A),
              A.min = min(A),
              A.max = max(A),
              A = mean(A)) %>%
    mutate(time = "ZTmean")
  counts.avg$genotype <- factor(as.character(counts.avg$genotype), levels = c("Liver_WT", "Liver_Cry1intronKO"))
  counts.avg$time <- counts.avg$genotype
  
  # Primetime ---------------------------------------------------------------
  
  
  iris$col <- c('firebrick', 'forestgreen', 'steelblue')[as.integer(iris$Species)]
  iris$size <- 4
  iris$alpha <- 1
  iris <- split(iris, iris$Species)
  
  # Here comes tweenr
  iris_tween <- iris$setosa %>% 
    tween_state(iris$versicolor, ease = 'cubic-in-out', nframes = 30) %>% 
    keep_state(10) %>% 
    tween_state(iris$virginica, ease = 'elastic-out', nframes = 30) %>% 
    keep_state(10) %>% 
    tween_state(iris$setosa, ease = 'quadratic-in', nframes = 30) %>% 
    keep_state(10)
  
  
  
  pos.max <- 1e5
  
  jsub <- subset(counts.long.merged, sig == jsig & genotype == "Liver_WT")
  
  jsub.tween <- tween_elements(jsub %>% mutate(ease = "linear"), "time", "pos", "ease", nframes = 300)
  gapminder_tween <- tween_elements(gapminder_edit,
                                    "time", "id", "ease", nframes = 300)
  
    m1.WT <- PlotSignalLR(jbait, jsub, pseudo.low=500, mindist = pos.max, jtitle=paste0(jbait, " Around the Clock"), ltype = "solid") + 
      facet_wrap(~genotype) + 
      theme(aspect.ratio = jaspect.ratio, strip.background = element_blank(),strip.text.y = element_blank()) + 
      scale_color_manual(values = hex.cols) +
      scale_x_continuous(labels=bToKb()) + 
      labs(x = "Position relative to bait [kb]",
           y = expression('4C signal [log'[10]*']')) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m1.WT)
    
    
    print(m1.WT + coord_cartesian(ylim = c(0, 0.6)))
    print(m1.WT + geom_vline(xintercept = posrange, linetype = "dotted"))
    m1.KO <- PlotSignalLR(jbait, subset(counts.long.merged, sig == jsig & genotype == "Liver_Cry1intronKO"), pseudo.low=500, mindist = pos.max, jtitle=paste0(jbait, " Around the Clock"), ltype = "solid") + 
      facet_wrap(~genotype) + 
      theme(aspect.ratio = jaspect.ratio,strip.background = element_blank(),strip.text.y = element_blank()) + 
      scale_color_manual(values = hex.cols) + 
      labs(x = "Position relative to bait [kb]",
           y = expression('4C signal [log'[10]*']')) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m1.KO)
    print(m1.KO + coord_cartesian(ylim = c(0, 0.6)))
    print(m1.KO + geom_vline(xintercept = posrange, linetype = "dotted"))
    
    multiplot(m1.WT+ geom_vline(xintercept = posrange, linetype = "dotted"), m1.KO+ geom_vline(xintercept = posrange, linetype = "dotted"), cols = 1)
    
    # m2 <- PlotSignalLR(jbait, subset(counts.long.merged, sig == jsig), pseudo.low=500, mindist = pos.max, jtitle=paste0(jbait, " Around the Clock")) + 
    #   facet_wrap(~genotype, ncol = 1) + 
    #   theme(aspect.ratio = jaspect.ratio) + 
    #   scale_color_manual(values = hex.cols) + 
    #   labs(x = "Position relative to bait [kb]",
    #        y = expression('4C signal [log'[10]*']')) + 
    #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # print(m2)
    # print(m2 + coord_cartesian(ylim = c(0, 0.6)))
    # print(m2 + geom_vline(xintercept = posrange, linetype = "dotted"))
    # 
    # for (jdo.facet in c(TRUE, FALSE)){
    #   m1 <- PlotAmps4C(counts.complex, do.facet = jdo.facet, pos.max = pos.max, .log2 = TRUE) + theme(aspect.ratio = jaspect.ratio / 2) + 
    #     labs(x = "Position relative to bait [kb]",
    #          y = expression('Fold change [log'[2]*']')) + 
    #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    #   print(m1)
    #   print(m1 + xlim(-100e3, 50e3) + coord_cartesian(xlim = c(-100e3, 100e3)) + theme(axis.text.x=element_blank()))
    #   # print(m1 + coord_cartesian(ylim = c(0, 1)))
    #   print(m1 + geom_vline(xintercept = posrange, linetype = "dotted"))
    #   
    #   # Plot p-value significance of amplitude
    #   m2 <- PlotPvalChiSqr(counts.complex, do.facet = jdo.facet, pos.max = pos.max) + scale_y_continuous(labels=scale1decimal)  + theme(aspect.ratio = jaspect.ratio / 2) + 
    #     labs(x = "Position relative to bait [kb]",
    #          y = expression('Rhythomicity [log'[10]*'p]')) + 
    #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    #   print(m2)
    #   print(m2 + xlim(-100e3, 50e3) + coord_cartesian(xlim = c(-100e3, 100e3)) + theme(axis.text.x=element_blank()))
    #   print(m2 + coord_cartesian(ylim = c(0, 8.5)))
    #   print(m2 + xlim(-100e3, 50e3) + coord_cartesian(xlim = c(-100e3, 100e3), ylim = c(0, 8.5)) + theme(axis.text.x=element_blank())) 
    #   print(m2 + geom_vline(xintercept = posrange, linetype = "dotted"))
    # }
  



