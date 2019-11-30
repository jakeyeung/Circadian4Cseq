# Jake Yeung
# Date of Creation: 2017-09-21
# File: ~/projects/4c_seq/scripts/functions/FitRhythmic4CSmoothed.R
# Fit rhythmic 4C but smooth of n.fragments

# Fit Fragment by Fragment

Project4CtoZscore <- function(counts.long.copy, jpval.k = 0.005, jpval.k.chisqr = 1e-5, match.to.chisqr.color = TRUE){
  counts.long.copy$gene <- as.character(counts.long.copy$pos)
  counts.long.copy$experiment <- "4CSeq"
  counts.long.copy$exprs <- counts.long.copy$A
  counts.long.copy$zttime <- counts.long.copy$time
  counts.long.copy$time <- as.numeric(gsub("ZT", "", counts.long.copy$time))
  counts.long.copy$tissue <- counts.long.copy$genotype
  counts.long.copy$se <- sqrt(counts.long.copy$Avar * counts.long.copy$sig2.adj)
  # counts.long.copy$se <- counts.long.copy$sig2.adj
  
  # Find rhythmic peaks by propagating the errorrs
  counts.complex <- ProjectWithZscore(counts.long.copy, omega = 2 * pi / 24, n=4)
  counts.complex$pos <- as.numeric(counts.complex$gene)
  
  # Handle LR in case LR is not defined by 0 due to shifting of coordinates
  counts.complex <- dplyr::inner_join(counts.complex, subset(counts.long.copy, select = c(pos, tissue, LR)))
  # remove duplicate rows
  counts.complex <- counts.complex[!duplicated(counts.complex), ]
  # counts.complex <- dplyr::left_join(counts.complex, subset(counts.long.copy, select = c(pos, tissue, LR)))
  # counts.complex$LR <- sapply(counts.complex$pos, function(p) ifelse(p < 0, "Left", "Right"))
  counts.complex$pval <- 2*pnorm(-abs(counts.complex$zscore))
  
  # set when dots go gray
  counts.complex$hex.col <- PhaseAmpPvalToColor(counts.complex$phase, counts.complex$amp, counts.complex$pval, rotate.hr = -8, amp.k = 0.01, pval.k = jpval.k, method = "cutoff")
  counts.complex$hex.col2 <- sapply(counts.complex$hex.col, function(h) ifelse(h == "#FFFFFF", "gray85", h))
  counts.complex$amp.zscore <- mapply(function(amp, amp.se) amp / amp.se, counts.complex$amp, counts.complex$amp.se)
  
  # add chisquare test for significance of amplitude (real and imaginary parts are gaussian variables, so squared sums should be chi-square distributed)
  counts.complex$real.z <- Re(counts.complex$exprs.transformed) / counts.complex$se.real
  counts.complex$im.z <- Im(counts.complex$exprs.transformed) / counts.complex$se.im
  counts.complex$ssz <- counts.complex$real.z ^ 2 + counts.complex$im.z ^ 2
  counts.complex$chisqr.pval <- 1 - pchisq(counts.complex$ssz, 2) 
  
  # redo color
  hex.col.chisqr <- PhaseAmpPvalToColor(counts.complex$phase, counts.complex$amp, counts.complex$chisqr.pval, rotate.hr = -8, amp.k = 0.01, pval.k = jpval.k.chisqr, method = "cutoff")
  counts.complex$hex.col.chisqr <- sapply(hex.col.chisqr, function(h) ifelse(h == "#FFFFFF", "gray85", h))
  
  if (match.to.chisqr.color){
    counts.complex$hex.col2 <- counts.complex$hex.col.chisqr
  }
  return(counts.complex)
}


# modified for 4CSeq
# 2017-09-21

FormatDat <- function(counts.sub, posrange, log10.to.log2){
  pos.start <- posrange[[1]]; pos.end <- posrange[[2]]
  region <- paste0(pos.start, pos.end, collapse = "to")
  counts.sub$pos <- as.character(counts.sub$pos)
  counts.sub$time <- as.numeric(sapply(counts.sub$time, function(x) gsub("ZT", "", x)))
  if (log10.to.log2){
    counts.sub$A <- counts.sub$A / (log10(2))
  }
  return(counts.sub)
}

FormatAndFitRhythmAcrossFrags <- function(counts.long, posrange, label=NULL, log10.to.log2=FALSE, fit.first.last.pos=FALSE, fit.every.frag=FALSE){
  # position should be in character, not numeric
  # make region.name that is reasonable
  pos.start <- posrange[[1]]; pos.end <- posrange[[2]]
  counts.sub <- subset(counts.long, pos > pos.start & pos < pos.end)
  
  if (fit.first.last.pos){
    # Because positions nearby are "correlated" due to sigma smoothing, fit on 
    # first and last position of position to have a more fair estimate of p-value
    counts.sub <- subset(counts.sub, pos %in% range(counts.sub$pos))
  }
  
  counts.sub <- FormatDat(counts.sub, posrange, log10.to.log2)
  region <- paste0(pos.start, pos.end, collapse = "to")
  # counts.sub$pos <- as.character(counts.sub$pos)
  # counts.sub$time <- as.numeric(sapply(counts.sub$time, function(x) gsub("ZT", "", x)))
  # if (log10.to.log2){
  #   counts.sub$A <- counts.sub$A / (log10(2))
  # }
  jfit <- FitRhythmAcrossFrags(counts.sub, T.period=24, region.name=region, get.residuals = TRUE, signal.cname="A", label=label, window.cname = "pos", time.cname = "time")
  if (log10.to.log2){
    # show also log2FC, which is 2 times the amplitude
    jfit$logfc <- 2 * jfit$amp
  }
  return(jfit)
}

FitRhythmAcrossFrags <- function(dat, T.period, region.name, get.residuals = TRUE, signal.cname="A", label=NULL, window.cname = "pos", time.cname = "time", single.frag=FALSE){
  # Fit rhythmic model to long dat. Expect dat to be per tissue per gene.
  # use lapply and ddply to vectorize this code.
  # dat needs columns: tissue time experiment exprs
  
  
  # Get parameters from complex model: used later
  GetParamsRhythModel <- function(myfit){
    model.params <- coef(myfit)
    params.b <- model.params[length(model.params)]
    params.a <- model.params[(length(model.params) - 1)]
    params.int <- model.params[1:(length(coef(rhyth.fit)) - 2)]
    windows <- names(model.params[1:(length(coef(rhyth.fit)) - 2)])
    windows <- sapply(windows, function(w) gsub("window", "", w))
    return(list(intercept = params.int,
                a = params.a,
                b = params.b,
                windows = windows))
  }
  if (!single.frag){
    form.rhyth <- as.formula(paste0(signal.cname, " ~ 0 + ", window.cname, " + cos(w * ", time.cname, ") + sin(w * ", time.cname, ")"))
    form.flat <- as.formula(paste0(signal.cname, " ~ 0 + ", window.cname))
  } else {
    form.rhyth <- as.formula(paste0(signal.cname, " ~ 1 + cos(w * ", time.cname, ") + sin(w * ", time.cname, ")"))
    form.flat <- as.formula(paste0(signal.cname, " ~ 1"))
  }
  
  w = 2 * pi / T.period
  rhyth.fit <- lm(form.rhyth, data = dat)
  flat.fit <- lm(form.flat, data = dat)
  compare.fit <- anova(flat.fit, rhyth.fit)
  pval <- compare.fit["Pr(>F)"][[1]][2]
  model.params <- GetParamsRhythModel(rhyth.fit)  # y = experimentarray + experimentrnaseq + a cos(wt) + b sin(wt)
  amp <- sqrt(model.params$a ^ 2 + model.params$b ^ 2)
  phase.rad <- atan2(model.params$b, model.params$a)
  phase.time <- (phase.rad / w) %% T.period
  if (get.residuals){
    # only care about residuals: do this if you want to see how the fit varies by changing period.
    ssq.residuals <- anova(rhyth.fit)["Residuals", "Sum Sq"]
    ssq.residuals.flat <- anova(flat.fit)["Residuals", "Sum Sq"]
    variance <- ssq.residuals / nrow(dat) # both rhyth and flat the same
    #     chi_sqr.rhyth <- ssq.residuals / variance
    #     chi_sqr.flat <- ssq.residuals.flat / variance
    
    #     dat.out <- data.frame(ssq.residuals = ssq.residuals, variance = variance)
    dat.out <- data.frame(region = region.name,
                          cos.part = model.params$a, sin.part = model.params$b,
                          amp = amp, phase = phase.time, pval = pval,
                          ints = as.character(paste0(model.params$intercept, collapse=",")),
                          ssq.residuals = ssq.residuals,
                          ssq.residuals.flat = ssq.residuals.flat,
                          variance = variance,
                          ints.mean = mean(model.params$intercept),
                          stringsAsFactors = FALSE)
  } else {
    dat.out <- data.frame(tissue = tissue,
                          cos.part = model.params$a, sin.part = model.params$b,
                          amp = amp, phase = phase.time, pval = pval,
                          ints = model.params$intercept,
                          ints.mean = mean(model.params$intercept),
                          stringsAsFactors = FALSE)
  }
  # add label optionally
  if (!is.null(label)){
    dat.out$label <- label
  }
  return(dat.out)
}

FitRhythmRegion <- function(jsub, n.frags, ncores=10){
  # fit 4cseq signal across entire "region". Expect jsub to be single genotype, Left/Right of bait with ideally 6 timepoints
  jsub <- jsub %>%
    arrange(pos)
  n.time <- length(unique(jsub$time))
  
  jsub.wide <- dcast(jsub, "pos ~ time", value.var = "A")
  jenv <- new.env()
  # chunk matrix into environment and then run fits
  for (i in seq(1, nrow(jsub.wide) - n.frags)){
    jlab <- as.character(jsub.wide$pos[i + n.frags / 2])
    jsub.sub <- melt(jsub.wide[i:(i+n.frags), ], id.vars = "pos", variable.name = "time", value.name = "A")
    jsub.sub$time <- as.numeric(sapply(jsub.sub$time, function(x) gsub("ZT", "", x)))  # convert ZT00 to 0
    jsub.sub$pos <- as.character(jsub.sub$pos)
    jenv[[jlab]] <- jsub.sub
  }
  fits <- mclapply(ls(jenv), function(region.label){
    jsub <- get(region.label, jenv)
    fit <- FitRhythmAcrossFrags(jsub, T.period = 24, get.residuals = TRUE, region.name = as.numeric(region.label), signal.cname = "A", window.cname = "pos", time.cname = "time")
  }, mc.cores = ncores)
  fits <- do.call(rbind, fits)
  return(fits)
}
