FitRhythmicAssays <- function(dat, T.period, get.residuals = TRUE, signal.cname="signal.log10", label=NULL, window.cname = "window", time.cname = "time"){
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
  
  form.rhyth <- as.formula(paste0(signal.cname, " ~ 0 + ", window.cname, " + cos(w * ", time.cname, ") + sin(w * ", time.cname, ")"))
  form.flat <- as.formula(paste0(signal.cname, " ~ 0 + ", window.cname))

  region <- paste0(as.character(unique(dat$chromo)), ":", min(dat$start), "-", max(dat$end), collapse="")
  windows <- paste0(as.character())
  w = 2 * pi / T.period
  # Expect columns exprs and time
  # rhyth.fit <- lm(signal.log10 ~ 0 + window + cos(w * time) + sin(w * time), data = dat)
  # flat.fit <- lm(signal.log10 ~ 0 + window, data = dat)
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
    dat.out <- data.frame(region = region, 
                          cos.part = model.params$a, sin.part = model.params$b, 
                          amp = amp, phase = phase.time, pval = pval,
                          ints = as.character(paste0(model.params$intercept, collapse=",")),
                          windows = as.character(paste0(model.params$windows, collapse=",")),
                          ssq.residuals = ssq.residuals,
                          ssq.residuals.flat = ssq.residuals.flat,
                          variance = variance)
  } else {
    dat.out <- data.frame(tissue = tissue, 
                          cos.part = model.params$a, sin.part = model.params$b, 
                          amp = amp, phase = phase.time, pval = pval,
                          ints = model.params$intercept)
  }
  # add label optionally
  if (!is.null(label)){
    dat.out$label <- label
  }
  return(dat.out)
}


