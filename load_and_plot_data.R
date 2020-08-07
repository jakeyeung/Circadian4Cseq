#' ## Installation: 
# devtools::install_github("jakeyeung/Circadian4Cseq")


# Jake Yeung
# Date of Creation: 2019-11-30
# File: ~/projects/Circadian4Cseq/load_and_plot_data.R
# Load and plot data


library(Circadian4Cseq)

library(dplyr)
library(ggplot2)
library(PhaseHSV)
library(parallel)

# for creating GIFs
library(tweenr)
library(gganimate)


# Get colors --------------------------------------------------------------

zts <- seq(0, 20, 4)
rotate.hrs <- -8
zts.rotated <- RotatePhase(zts, rotate.hr = rotate.hrs)

hex.cols <- hsv(PhaseToHsv(zts.rotated, 0, 24), 1, 1)

rgb.cols <- HexColorToRgb(hex.cols)

jposes <- c(0)




# Load data  --------------------------------------------------------------

jbaits <- c("Hoxd4", "Cry1", "Cry1_TSS", "Cry1_UP")
posranges <- list(c(-15000, 15000), c(-29000, -19000), c(-31000, -22200), c(-36000, -29000))
jpval.k <- 0.01  # make more stringent?
jpval.k.chisqr <- 0.01  # colors points for Cry1_UP

jaspect.ratio <- 0.25

lst <- list()
for (i in seq(length(jbaits))){
  lst[[i]] <- list("jbait"=jbaits[[i]], "posrange"=posranges[[i]])
}

# Do Cry1 as example
indx <- 2

sublst <- lst[[indx]]
jbait <- sublst[["jbait"]]
posrange <- sublst[["posrange"]]


jsig <- 2500

load("data/Cry1.normmax.Inf.log.FALSE.nfrags.5.remove.auto.sig.2500.Robj.RData", v=T)


counts.long.merged <- counts.long; rm(counts.long)

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


#' ## Plot 4C-seq signal over the 24-hour day 

pos.max <- 1e5

jsub <- subset(counts.long.merged, sig == jsig & genotype == "Liver_WT")

jsub$time.numer <- as.numeric(sapply(as.character(jsub$time), GetSampTime, remove.prefix=TRUE))


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



