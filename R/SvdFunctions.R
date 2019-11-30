# SvdFunctions.R
# to be used in svd.on.rhythmic.genes.R and other similar scripts
# March 2 2015

# library(plyr)
library(reshape2)
# library(PhaseHSV)
library(ggplot2)

GetAmpPhaseFromActivities <- function(act.l, mrna.hl, jtiss = "Liver", jgeno = "WT"){
  act.l.complex <- ProjectWithZscore(act.l, omega, n = 4)
  act.l.complex$phase <- AdjustPhase(act.l.complex$phase, half.life = mrna.hl, fw.bw = "bw")
  return(act.l.complex)
}

ProjectWithZscore <- function(act.long, omega, n = 4){
  act.complex <- act.long %>%
    group_by(gene, tissue) %>%
    do(ProjectToFrequency2(., omega, add.tissue=TRUE, propagate.errors = TRUE)) %>%
    # zscore from 0
    mutate(amp = GetAmp(a = Re(exprs.transformed), b = Im(exprs.transformed), n = n),  # peak to trough
           phase = GetPhi(a = Re(exprs.transformed), b = Im(exprs.transformed), omega),  # in hours
           amp.se = GetAmp.se(a = Re(exprs.transformed), b = Im(exprs.transformed), sig.a = se.real, sig.b = se.im, n = n),
           phase.se = GetPhi.se(a = Re(exprs.transformed), b = Im(exprs.transformed), sig.a = se.real, sig.b = se.im, omega),
           zscore = amp / amp.se) %>%
    arrange(desc(zscore))
  return(act.complex)
}

ProjectToFrequency2 <- function(dat, omega, add.tissue=FALSE, propagate.errors=FALSE){
  # simpler than ProjectToFrequency().
  # expect dat to be gene i in tissue c with column names time and exprs
  exprs.transformed <- DoFourier(dat$exprs, dat$time, omega = omega)
  if (add.tissue){    
    tissue <- unique(dat$tissue)
    dat.out <- data.frame(tissue = tissue, exprs.transformed = exprs.transformed)
  } else {
    dat.out <- data.frame(exprs.transformed = exprs.transformed)
  }
  if (propagate.errors){
    se.xy <- DoFourier.se(dat$se, dat$time, omega = omega, normalize = TRUE)
    dat.out$se.real <- se.xy$se.cospart
    dat.out$se.im <- se.xy$se.sinpart
  }
  return(dat.out)
}

ProjectToFrequency <- function(dat, my.omega, normalize = TRUE, rhythmic.only = FALSE, method = "ANOVA", pval.cutoff = 5e-3){
  # Perform fourier transform and normalize across all frequencies (squared and square root)
  # 
  # Input:
  # dat: long dataframe containing expression and time columns for one condition. Expect 'exprs' and 'time' columns.
  # to be transformed to frequency domain.
  # my.omega: the omega of interest.
  # normalize: converts transform into a sort of z-score
  # rhytmic.only: transforms only genes that are rhythmic (by BIC method)
  # 
  # omega in which we are interested.  
  if (rhythmic.only){
    if (IsRhythmic(dat, my.omega, pval.cutoff = pval.cutoff, method = "ANOVA")){
      dat.transformed <- Transform(dat, my.omega, normalize)
    } else {
      dat.transformed <- data.frame(exprs.transformed = 0)
    }
  } else {
    dat.transformed <- Transform(dat, my.omega, normalize)
  }
  return(dat.transformed)
}



Transform <- function(dat, my.omega, normalize = TRUE){
  # Perform fourier transform and normalize across all frequencies (squared and square root)
  # 
  # Input:
  # dat: long dataframe containing expression and time columns for one condition. Expect 'exprs' and 'time' columns.
  # to be transformed to frequency domain.
  # my.omega: the omega of interest.
  # normalize: converts transform into a sort of z-score
  # 
  # omega in which we are interested.
  
  # if we normalize, get list of omegas.
  # otherwise we just use omega
  if (normalize){
    # t <- sort(unique(dat$time))
    # n.timepoints <- length(t)
    # interval <- t[2] - t[1]
    omegas <- GetOmegas(remove.zero = TRUE)
  } else {
    omegas <- my.omega
  }
  
  transforms <- sapply(omegas, DoFourier, exprs = dat$exprs, time = dat$time)
  my.transformed <- transforms[which(omegas == omega)]  # corresponds to omega of interest
  
  if (normalize){
    # Normalize across omega
    # if median is 0, then set factor to 0, otherwise 1
    factor <- 1
    cutoff <- 5
    jmedian <- median(subset(dat, experiment == "rnaseq")$exprs)
    if (jmedian <= cutoff){
      factor <- 0
    }
    
    my.transformed <- (my.transformed / sqrt(sum(Mod(transforms) ^ 2))) * factor
  } 
  return(data.frame(exprs.transformed = my.transformed))
}

GetInterval <- function(time.vector){
  # Given vector of times (equally spaced time points), return the time interval
  time.sort <- sort(unique(dat$time))
  interval <- time.sort[2] - time.sort[1]
  return(interval)
}

IsRhythmic <- function(dat, my.omega, pval.cutoff = 5e-3, method = "ANOVA"){
  # Test if rhythmic by BIC model selection: fit through all omegas.
  # dat: long format, gene and condition
  # my.omega: omega of interest
  # method = "BIC" or "ftest"
  # pval.cutoff: for method = "ftest"
  
  fit.rhyth <- lm(exprs ~ 0 + experiment + sin(my.omega * time) + cos(my.omega * time), data = dat)
  fit.flat <- lm(exprs ~ 0 + experiment, data = dat)  # intercept only
  
  if (method == "BIC"){
    omegas.all <- GetOmegas(remove.zero = TRUE)
    # remove also my.omega to get omegas representing "noise"
    omegas.noise <- omegas.all[which(omegas.all != my.omega)]
    
    bic.test <- BIC(fit.rhyth, fit.flat)
    
    chosen.model <- rownames(bic.test)[which(bic.test$BIC == min(bic.test$BIC))]
    
    if (chosen.model == "fit.flat"){
      # if flat, no need to check for noise.
      return(FALSE)
    }
    # chosen model is fit.rhyth, but is it a noisy gene?
    # Check if noise components have an even better fit than fit.rhyth
    rhyth.bic <- bic.test["fit.rhyth", "BIC"]
    # fit noise
    fits.noise <- lapply(omegas.noise, 
                         function(w) lm(exprs ~ sin(w * time) + cos(w * time), 
                                        data = dat))
    bic.noise.min <- min(sapply(fits.noise, BIC))
    # print(paste(unique(dat$gene), "noise bic:", bic.noise.min, "rhyth bic:", rhyth.bic))
    if (rhyth.bic < bic.noise.min){
      return(TRUE)
    } else {
      return(FALSE)
    } 
  } else if (method == "ANOVA"){
    f.test <- anova(fit.flat, fit.rhyth)
    pval <- f.test[["Pr(>F)"]][[2]]
    if (is.nan(pval)){
      # if gene is all flat, then pval is NaN, force it to be 1 in this case.
      pval <- 1
    }
    if (pval < pval.cutoff){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}


DoFourier <- function(exprs, time, omega, normalize = TRUE){
  p <- exprs %*% exp(1i * omega * time)
  if (normalize){
    p <- p / length(time)  # normalize by number of datapoints 
  }
  return(p)
}

DelFA.DelI <- function(x, omega = 2 * pi / 24){
  # error propagation needs partial derivatives of Fourier with respect to I (A and B)  
  # FA: cospart
  # FB: sinepart
  N <- length(x)
  delFa.delI <- (1 / N) * cos(omega * x)
  return(delFa.delI)
}

DelFB.DelI <- function(x, omega = 2 * pi / 24){
  # error propagation needs partial derivatives of Fourier with respect to I (A and B)  
  # FA: cospart
  # FB: sinepart
  N <- length(x)
  delFb.delI <- (1 / N) * sin(omega * x)
  return(delFb.delI)
}

DelF.DelI <- function(x, omega = 2 * pi / 24){
  # error propagation using exponential form
  N <- length(x)
  delF.delI <- (1 / N) * exp(1i *omega * x)
  return(delF.delI)
}

DoFourier.se <- function(se, jtime, omega = 2 * pi / 24, normalize = TRUE){
  # http://www.tina-vision.net/teaching/cvmsc/docs/fourier.ps
  # Error Propagation and the Fourier Transform􏰀
  N <- length(se)
  se.cospart <- sqrt(t(DelFA.DelI(jtime, omega)) %*% diag(se ^ 2) %*% DelFA.DelI(jtime, omega))
  se.sinpart <- sqrt(t(DelFB.DelI(jtime, omega)) %*% diag(se ^ 2) %*% DelFB.DelI(jtime, omega))
  return(list(se.cospart = se.cospart, se.sinpart = se.sinpart))
}

GetOmegas <- function(n.timepoints=24, interval=2, remove.zero=FALSE){
  # Get omegas from a time series
  # 
  # remove.zero: removes first element (omega = 0).
  # Useful for normalizing Fourier when I don't want omega = 0 case
  
  # Output:
  # list of omegas: can be used for fourier transform
  
  t <- seq(from = 0, by = 1, length.out = n.timepoints)  # will account for interval later
  
  # get harmonic frequencies
  # this creates harmonic frequencies from 0 to .5 in steps of 1/n.timepoints
  
  mid.point <- length(t) / 2 + 1
  freqs <- t[1:mid.point] / n.timepoints
  
  # Convert to period before adjusting for interval
  T.unadj <- 1 / freqs
  
  # Account for interval
  T <- T.unadj * interval
  
  # Calculate omegas
  omegas <- 2 * pi / T
  
  if (remove.zero){
    omegas <- omegas[2:length(omegas)]
  }
  return(omegas)
}

PlotComponentOriginal <- function(orig.dat, s, jgene, component, show.original = TRUE){
  component.mat <- s$d[component] * OuterComplex(s$u[, component, drop = FALSE], t(s$v[, component, drop = FALSE]))
  if (show.original){
    PlotComplex(as.matrix(orig.dat[jgene, ]), labels = colnames(orig.dat), main = paste(jgene, "original"), add.text.plot = FALSE) 
  }
  PlotComplex(as.matrix(component.mat[jgene, ]), labels = colnames(component.mat), main = paste(jgene, "component", component), add.text.plot = FALSE)
}

ProjectOnVector <- function(input.complex, basis.complex){
  # Given an input complex number, project onto a basis complex number using
  # simple trigonometry or dot product
  # 
  # Returns length of vector after projection. >0 indicates along basis vector.
  # <0 indicates anti-correlation with basis vector.
  
  if (is.data.frame(input.complex)){
    input.complex <- as.matrix(input.complex)
  }
  if (is.data.frame(basis.complex)){
    basis.complex <- as.matrix(basis.complex)
  }
  phase.diff <- Arg(input.complex) - Arg(basis.complex)
  projected.mod <- Mod(input.complex) * cos(phase.diff)
  return(projected.mod)
}

TemporalToFrequency <- function(dat, period = 24){
  library(plyr)
  omega <- 2 * pi / period
  
  start.time <- Sys.time()
  dat.complex <- lapply(split(dat, dat$tissue), function(x){
    ddply(x, .(gene), ProjectToFrequency2, omega = omega, add.tissue = TRUE)
  }) %>%
    do.call(rbind, .) %>%
    mutate(magnitude = Mod(exprs.transformed)) %>%
    arrange(desc(magnitude))
  print(Sys.time() - start.time)
  
  detach("package:plyr", unload=TRUE)
  library(dplyr)
  return(dat.complex)
}

GetFourierEntropy <- function(dat, n, interval, method = "direct"){
  # Calculates entropy of all possibl fourier components
  # method either "direct" or "fft". Direct allows uneven sampling.
  # if fft, then it does not care about timing, it's all assumed to be even sampling
  if (method == "direct"){
    # here we are limited to the lowest sampling rate which is rnaseq
    fund.T <- n * interval
    periods <- fund.T / seq(24 / interval)  # for assessing "noise"
    weights <- vector(mode = "numeric", length = length(periods))
    for (i in seq(length(periods))){
      source("~/projects/tissue-specificity/scripts/functions/ShannonEntropy.R")
      p <- periods[i]
      # loop through harmonics and add up weights
      o <- 2 * pi / p
      weights[i] <- Mod(DoFourier(dat$exprs, dat$time, omega = o)) ^ 2
    }
    T.max <- periods[which(weights == max(weights))]
    # normalize by sum then calculate entropy
    weights.norm <- weights / sum(weights)
    H <- ShannonEntropy(weights.norm)
    H.max <- log2(length(weights.norm))
    log.H.weight <- log(H.max / H)
  }
  else if (method == "fft"){
    # use fft instead of direct method because only on array gives us more harmonics maybe
    # gives more accurate results of the "noise"
    source("~/projects/tissue-specificity/scripts/functions/FourierFunctions.R")
    freq.per <- CalculatePeriodogram(dat$exprs)
    # ignore freq = 0
    freq <- freq.and.periodogram$freq[2:length(freq.and.periodogram$freq)]
    periods <- (1 / freq) * interval
    per <- freq.and.periodogram$p.scaled[2:length(freq.and.periodogram$freq)]
    per.norm <- per / sum(per)
    H <- ShannonEntropy(per.norm)
    H.max <- log2(length(per.norm))
    f.max <- FindMaxFreqs(freq, per)[1]  # take top
    T.max <- (1 / f.max) * interval
    log.H.weight <- log(H.max / H)
  }
  return(list(log.H.weight=log.H.weight, T.max=T.max, H.max=H.max, H=H))
}

GetNoiseFactor <- function(dat, period, n, interval, method = "direct"){
  # Calculate weight of other components to give a weight for period of interest
  # method either "direct" or "fft". Direct allows uneven sampling.
  # if fft, then it does not care about timing, it's all assumed to be even sampling
  if (missing(period)){
    period <- 24  # hrs
  }
  if (method == "direct"){
    # here we are limited to the lowest sampling rate which is rnaseq
    fund.T <- n * interval
    periods <- fund.T / seq(24 / interval)  # for assessing "noise"
    weights <- vector(mode = "numeric", length = length(periods))
    for (i in seq(length(periods))){
      source("~/projects/tissue-specificity/scripts/functions/ShannonEntropy.R")
      p <- periods[i]
      # loop through harmonics and add up weights
      o <- 2 * pi / p
      weights[i] <- Mod(DoFourier(dat$exprs, dat$time, omega = o)) ^ 2
    }
    T.max <- periods[which(weights == max(weights))]
    # normalize by sum then calculate entropy
    weights.norm <- weights / sum(weights)
    frac.weight <- weights.norm[which(periods == period)]
    sum.weight <- sum(weights)
  }
  else if (method == "fft"){
    # use fft instead of direct method because only on array gives us more harmonics maybe
    # gives more accurate results of the "noise"
    source("~/projects/tissue-specificity/scripts/functions/FourierFunctions.R")
    freq.per <- CalculatePeriodogram(dat$exprs)
    # ignore freq = 0
    freq <- freq.and.periodogram$freq[2:length(freq.and.periodogram$freq)]
    periods <- (1 / freq) * interval
    per <- freq.and.periodogram$p.scaled[2:length(freq.and.periodogram$freq)]
    per.norm <- per / sum(per)
    f.max <- FindMaxFreqs(freq, per)[1]  # take top
    T.max <- (1 / f.max) * interval
    frac.weight <- per.norm[which(periods == period)]
    sum.weight <- sum(per)
  }
  return(list(frac.weight = frac.weight, T.max=T.max, sum.weight = sum.weight))
}

TemporalToFrequency2 <- function(dat, period = 24, n = 8, interval = 6, add.entropy.method = "array"){
  # use dplyr
  # n: number of non-redundant timepoints over interval hours
  # default n = 8, interval = 6 because we cannot reliable estimate 
  # because we take the lowest sampling rate which is rnaseq
  # add.entropy.method "array" | "both". Array may be better for assessing noise.
  
  # init by checking for exprs = 0, ignore those man
  if (max(dat$exprs) == 0){
    return(data.frame(NULL))
  }
  omega <- 2 * pi / period
  # get the main fourier component
  dat.complex <- ProjectToFrequency2(dat, omega, add.tissue = TRUE)
  # get other fourier components  
  if (add.entropy.method == "direct"){
    H.list <- GetNoiseFactor(dat, period, n, interval, method = "direct")
  }
  else if (add.entropy.method == "array"){
    dat.array <- subset(dat, experiment == "array")
    n.array <- 24  # n.timepoints
    n.interval <- 2  # 2 hrs per timepoint
    H.list <- GetNoiseFactor(dat.array, period, n.array, n.interval, method = "direct")
  }
  else if (add.entropy.method == FALSE){
    return(dat.complex)
  } else {
    warning("Must be 'array' or 'direct' or 'FALSE' for add.entropy.method")
    print(add.entropy.method)
  }
  # using weight approach
  dat.complex$T.max <- H.list$T.max
  dat.complex$frac.weight <- H.list$frac.weight
  dat.complex$sum.weight <- H.list$sum.weight
#   # using entropy approach
#   dat.complex$log.H.weight <- H.list$log.H.weight
#   dat.complex$H.max <- H.list$H.max
#   dat.complex$H <- H.list$H
  return(dat.complex)
}

TemporalToFrequencyDatLong <- function(dat.long, period = 24, n = 8, interval = 6, add.entropy.method = "array"){
  # wrapper to run in parallel
  library(parallel)
  dat.long.by_genetiss <- group_by(dat.long, gene, tissue)
  dat.long.by_genetiss.split <- split(dat.long.by_genetiss, dat.long.by_genetiss$tissue)
  # start <- Sys.time()
  dat.complex.split <- mclapply(dat.long.by_genetiss.split, function(jdf){
    rhyth <- jdf %>%
      group_by(gene) %>%
      do(TemporalToFrequency2(., period, n, interval, add.entropy.method))
  }, mc.cores = 12)
  dat.complex <- do.call(rbind, dat.complex.split)
  # print(Sys.time() - start)
  return(dat.complex)
}


GetEigens <- function(s.complex, period, comp = 1, xlab = "Amp", ylab = "Phase", label.n=30, 
                      eigenval = TRUE, adj.mag = TRUE, pretty.names = FALSE, constant.amp = FALSE, 
                      peak.to.trough = TRUE, jtitle, label.gene=NA, dot.col = "gray85", jsize = 22, dotsize = 1.5,
                      dotshape = 18, disable.text = FALSE, add.arrow = FALSE, disable.repel = FALSE, half.life = 0){
  source("~/projects/tissue-specificity/scripts/functions/PlotFunctions.R")
  if (missing(period)){
    period <- 24
  }
  omega <- 2 * pi / period
  
  var.explained <- s.complex$d ^ 2 / sum(s.complex$d ^ 2)
  eigengene <- s.complex$v[, comp]
  eigensamp <- s.complex$u[, comp]
  # 2015-10-15: Eigengene v is the complex conjugate, unconjugate it to get proper interpretation of the phase
  # otherwise your tissue module will be flipped
  eigengene <- Conj(eigengene)
  if (eigenval){
    eigensamp <- eigensamp * s.complex$d[comp]
  }
  # BEGIN: rotate to phase of largest magnitude in sample of eigengene
  phase.reference <- Arg(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
  rotate.factor <- complex(modulus = 1, argument = phase.reference)
  # rotate eigengene by -phase ref
  eigengene <- eigengene * Conj(rotate.factor)
  # rotate eigensamp by +phase ref
  eigensamp <- eigensamp * rotate.factor
  # END: rotate to phase of largest magnitude in sample of eigengene
  
  if (half.life > 0){
    # if half.life > 0, then adjust eigensamp by rotation corresponding to half-life of mRNA species
    gamma.mrna <- half.life / log(2)
    k.mrna <- 1 / gamma.mrna
    delay.rads <- atan(omega / k.mrna)  # rads
    # subtract delay from eigensamp
    rotate.factor <- complex(modulus = 1, argument = delay.rads)
    eigensamp <- eigensamp * Conj(rotate.factor)
  }
  
  # BEGIN: adjust so largest magnitude in sample of eigengene = 1
  # amp.factor: increases half amplitude to full amplitude. Makes life interpretable.
  if (adj.mag){
    mag.reference <- Mod(eigengene[which(Mod(eigengene) == max(Mod(eigengene)))])
    eigengene <- eigengene / mag.reference
    eigensamp <- eigensamp * mag.reference
  }
  # END: adjust largest magnitude to 1
  
  # BEGIN: adjust mean to peak to peak to trough
  if (peak.to.trough){
    # only on GENE MODULE
    eigengene <- 1 * eigengene
    eigensamp <- 2 * eigensamp
    xlab <- "Log2 Fold Change"
    # ylab <- "Phase (CT)"
  } else {
    xlab <- "Amp (mean to peak)"
  }

  # BEGIN: optinally remove P2 from motif names
  if (pretty.names){
    source("~/projects/tissue-specificity/scripts/functions/RemoveP2Name.R")
    rownames(s.complex$u) <- sapply(rownames(s.complex$u), RemoveP2Name)
  }
  # END
  
  # BEGIN: for gene module: only label the top label.n genes
  names(eigensamp) <- rownames(s.complex$u)
  eigensamp.sorted <- eigensamp[order(Mod(eigensamp), decreasing = TRUE)]
  names.orig <- names(eigensamp.sorted)
  if (label.n < length(eigensamp.sorted)){
    names(eigensamp.sorted)[(label.n + 1):length(eigensamp.sorted)] <- ""
  }
  # END
  
  # BEGIN: for gene module: also include a list of genes if not NA
  if (!is.na(label.gene)){
    names(eigensamp.sorted)[which(names.orig %in% label.gene)] <- names.orig[which(names.orig %in% label.gene)]
  }
  # END
  
  gene.labels <- rownames(s.complex$u)
  if (missing(jtitle)){
    # jtitle1 <- paste0("Tissue Module ", comp, " (", length(gene.labels), " genes)\n", signif(var.explained[comp], 2), " of total circadian variance)")  
    jtitle1 <- "Tissue Module"
    jtitle2 <- paste0("Gene Module ", comp, " (", length(gene.labels), " genes)\n", signif(var.explained[comp], 2), " of total circadian variance)")
  } else {
    jtitle1 <- jtitle
    jtitle2 <- jtitle
  }
  v.plot <- PlotComplex2(eigengene, labels = rownames(s.complex$v), omega = omega, 
                         title = jtitle1, 
                         xlab = "Tissue Weights", ylab = ylab, ampscale = 1, constant.amp = constant.amp,
                         dot.col = "black", jsize = jsize, add.arrow = add.arrow, disable.repel = FALSE)
    
  # u.plot <- PlotComplex2(eigensamp, labels = rownames(s.complex$u), omega = omega, title = paste0("Gene Module ", comp, " (", signif(var.explained[comp], 2), " of total circadian variance)"), xlab = xlab, ylab = ylab)
  u.plot <- PlotComplex2(eigensamp.sorted, labels = names(eigensamp.sorted), omega = omega, 
                         title = jtitle2, 
                         xlab = xlab, ylab = ylab, ampscale = 2, constant.amp = constant.amp, dot.col = dot.col, jsize = jsize, 
                         dotsize = dotsize, dotshape = dotshape, disable.text = disable.text,
                         add.arrow = add.arrow, disable.repel = disable.repel)
  return(list(v.plot = v.plot, u.plot = u.plot, eigengene = eigengene, eigensamp = eigensamp))
}

SvdOnComplex <- function(dat.complex, value.var = "exprs.transformed"){
  source("~/projects/tissue-specificity/scripts/functions/LongToMat.R")
  # used in heatmap_tissue_specific_rhythms.R
  M.complex <- LongToMat(dat.complex, value.var = value.var)
  s <- svd(M.complex)
  # add row and colnames
  rownames(s$u) <- rownames(M.complex)
  rownames(s$v) <- colnames(M.complex)
  # screeplot
  # plot(s$d ^ 2 / sum(s$d ^ 2), type = 'o')  # eigenvalues
  return(s)
}

NormalizeComplexMat <- function(dat, dic){
  key.tissue <- as.character(dat$tissue)[1]
  norm.factor <- dic[[key.tissue]]  # magnitude of reference gene (e.g. Arntl)
  dat$exprs.norm <- dat$exprs.transformed / norm.factor
  return(dat)
}
