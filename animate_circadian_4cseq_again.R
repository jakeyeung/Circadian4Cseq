# Jake Yeung
# Date of Creation: 2018-10-28
# File: ~/projects/4c_seq/scripts/primetime_animations/animate_around_the_clock.R
# Animate batch 3 around the clock
# Stolen from: ~/projects/4c_seq/scripts/batch3_4cseq_analysis/load_robj_plot_summary.batch3.WT_vs_Cry1intronKO.R


rm(list=ls())
jstart <- Sys.time()


library(here)
library(ggplot2)
library(gganimate)
library(tweenr)
library(transformr)
library(dplyr)
library(PhaseHSV)
library(Circadian4Cseq)



MakeAnimation <- function(jsub.orig, nf = 20, ks = 0, jfps = 20, vline = -28000,
                          plot.base = "/Users/yeung/projects/4c_seq_analysis/WT_plot"){
  plot.out <- paste0(plot.base, ".fixed.pdf")
  anim.out <- paste0(plot.base, ".fps.", jfps, ".fixed_rotated.11.5.gif")
  jsub <- split(jsub.orig, jsub.orig$time)
  # Here comes tweenr
  jsub_tween <- jsub$ZT00 %>% 
    tween_state(jsub$ZT04, ease = 'linear', nframes = nf) %>% 
    keep_state(ks) %>%
    tween_state(jsub$ZT08, ease = 'linear', nframes = nf) %>% 
    keep_state(ks) %>%
    tween_state(jsub$ZT12, ease = 'linear', nframes = nf) %>% 
    keep_state(ks) %>%
    tween_state(jsub$ZT16, ease = 'linear', nframes = nf) %>% 
    keep_state(ks) %>%
    tween_state(jsub$ZT20, ease = 'linear', nframes = nf) %>%
    keep_state(ks) %>%
    tween_state(jsub$ZT00, ease = 'linear', nframes = nf)
  
  # plot
  ltype <- "solid"; jtitle <- "ZT"
  jaspect.ratio <- 0.25
  zts <- seq(0, 20, 4)
  rotate.hrs <- -8
  zts.rotated <- RotatePhase(zts, rotate.hr = rotate.hrs)
  hex.cols <- hsv(PhaseToHsv(zts.rotated, 0, 24), 1, 1)
  # change yellow to darker yellow
  zt11.rotated <- RotatePhase(11.5, rotate.hr = rotate.hrs)
  zt11.col <- hsv(PhaseToHsv(zt11.rotated, 0, 24), 1, 1)
  # hex.cols[which(hex.cols == "#FFFF00")] <- "#FFD500"
  hex.cols[which(hex.cols == "#FFFF00")] <- zt11.col
  
  p.orig <- ggplot(jsub.orig, aes(x = pos, y = A, colour = time)) + 
    geom_line(data = subset(jsub.orig, LR == "Left"), aes(x = pos, y = A, colour = time), linetype = ltype) + 
    geom_line(data = subset(jsub.orig, LR == "Right"), aes(x = pos, y = A, colour = time), linetype = ltype) + 
    theme(aspect.ratio = 0.25) + geom_hline(aes(yintercept=0)) + 
    theme_bw() + ggtitle(jtitle) + geom_vline(aes(xintercept=vline), linetype="dotted") + 
    theme(legend.position = "bottom") + 
    xlab("Position relative to bait") + 
    ylab("Log10 Signal") + 
    theme(aspect.ratio = jaspect.ratio, strip.background = element_blank(),strip.text.y = element_blank()) + 
    scale_color_manual(values = hex.cols) +
    scale_x_continuous(labels=bToKb()) + 
    labs(x = "Position relative to bait [kb]",
         y = expression('4C signal [log'[10]*']')) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ggsave(filename = plot.out)
  
  jsub_tween <- jsub_tween %>%
    mutate(.ZT = signif(.frame * ((2 * pi) / 120) * (24 / (2 * pi)), digits = 2))
  
  # hex.cols.blue <- rep("#0000ff", 6)
  p <- ggplot(jsub_tween, aes(x = pos, y = A, colour = time)) + 
    geom_line(data = subset(jsub_tween, LR == "Left"), aes(x = pos, y = A, colour = time), linetype = ltype) + 
    geom_line(data = subset(jsub_tween, LR == "Right"), aes(x = pos, y = A, colour = time), linetype = ltype) + 
    theme(aspect.ratio = 0.25) + geom_hline(aes(yintercept=0)) + 
    theme_bw() + ggtitle(jtitle) + geom_vline(aes(xintercept=-28000), linetype="dotted") + 
    theme(legend.position = "bottom") + 
    xlab("Position relative to bait") + 
    ylab("Log10 Signal") + 
    theme(aspect.ratio = jaspect.ratio, strip.background = element_blank(),strip.text.y = element_blank()) + 
    scale_color_manual(values = hex.cols) +
    scale_x_continuous(labels=bToKb()) + 
    labs(x = "Position relative to bait [kb]",
         y = expression('4C signal [log'[10]*']'),
         # title = "ZT: {closest_state}") +
         title = "ZT: {frame_time}") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    transition_time(.ZT) +
    # transition_states(.ZT, transition_length = 1, state_length = 0) + 
    ease_aes('linear')
  animate(p, renderer = gifski_renderer(), fps = jfps)
  anim_save(filename = anim.out)
  return(p)
}



# Try to animate ----------------------------------------------------------


# load("/Users/yeung/projects/4c_seq_analysis/counts.long.merged.Robj", v=T)
data(counts.long.merged.Robj)


jsig <- 2500
pos.max <- 1e5
jsub.orig <- subset(counts.long.merged, sig == jsig & genotype == "Liver_WT" & abs(pos) < pos.max)

nf <- 20
ks <- 0
jfps <- 22

outdir <- "gifs_outputs_test"
dir.create(outdir)

p.WT <- MakeAnimation(subset(counts.long.merged, sig == jsig & genotype == "Liver_WT" & abs(pos) < pos.max),
                      nf = nf, ks = ks, jfps = jfps, 
                      plot.base = file.path(outdir, "WT_plot"))
p.KO <- MakeAnimation(subset(counts.long.merged, sig == jsig & genotype == "Liver_Cry1intronKO" & abs(pos) < pos.max),
                      nf = nf, ks = ks, jfps = jfps, 
                      plot.base = file.path(outdir, "KO_plot"))


print(p.WT)
print(p.KO)

# animate(p.WT, fps = 50, renderer = gifski_renderer())
# animate(p.KO, fps = 50, renderer = gifski_renderer())

