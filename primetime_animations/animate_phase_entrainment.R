# Jake Yeung
# Date of Creation: 2018-10-31
# File: ~/projects/4c_seq/scripts/primetime_animations/animate_phase_entrainment.R
# Animate phase entrainment

rm(list=ls())
library(ggplot2)


phi.deriv <- function(theta, phi, k, Tphi){
  return(2 * pi / Tphi + k * sin(theta - phi))
  # return(2 * pi / Tphi)
}
theta.deriv <- function(Ttheta){
  return(2 * pi / Ttheta)
}

theta0 <- 0
phi0 <- 0
Tphi <- 22
Ttheta <- 24
K <- 1


nsteps <- 50000
tsteps <- seq(nsteps)
tdelt <- 0.01
tvec <- seq(from = 0, by = tdelt, length.out = nsteps)
phi.vec <- rep(NA, length(tsteps))
phi.vec[1] <- phi0

theta.vec <- theta.deriv(Ttheta) * tvec

for (i in tsteps[2:length(tsteps)]){
  phi.vec[i] <- phi.vec[i - 1] + tdelt * phi.deriv(theta.vec[i - 1], phi.vec[i - 1], K, Tphi)
}

plot(phi.vec)

phasediff <- theta.vec - phi.vec

plot(sin(phasediff))

print(paste("Phase angle:", (Ttheta / (2 * pi)) * phasediff[length(phasediff)]))


# Plot on a circle --------------------------------------------------------

anim.color <- "grey50"
sun.color <- "#FFDF00"
dat.out <- data.frame(Time = c(rep(tsteps, 2)), 
                      phase = c(phi.vec, theta.vec) %% 2 * pi,
                      amp = c(rep(1, length(phi.vec)), rep(1.2, length(theta.vec))),
                      oscillator = c(rep("Animal", length(phi.vec)), rep("Sun", length(theta.vec))),
                      col = c(rep(anim.color, length(tsteps)), rep(sun.color, length(tsteps))))

dat.wide <- reshape2::dcast(dat.out, formula = )

indx.filt <- seq(1, nsteps, nsteps / 15)
tsteps.sub <- tsteps[indx.filt]
dat.sub <- subset(dat.out, Time %in% tsteps.sub)
jbreaks <- signif(seq(from = 0, to = 3*pi / 2, length.out = 4), 2)
jlabs <- c(0, expression(over(pi,4)), expression(over(pi,2)), expression(over(3 * pi, 4)))

ggplot() +
  geom_point(data = dat.sub, aes(x = phase, y = amp, fill = col, size = oscillator), shape = 21, color = "black") + 
  geom_segment(data = dat.sub, aes(y=0, x = phase, xend=phase, yend=amp)) + 
  geom_segment(data = dat.sub, aes(y=0, x = phase, xend=phase, yend=0.5)) + 
  coord_polar(theta = "x") + expand_limits(y = 0) + 
  scale_x_continuous(breaks = jbreaks, lim = c(0, 2*pi), labels = jlabs) + 
  scale_size_manual(values = c(2, 5)) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank()) + 
  scale_fill_identity() + 
  ylab("") + scale_y_continuous() + 
  xlab("")





