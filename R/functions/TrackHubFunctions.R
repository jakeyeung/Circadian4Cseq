library(PhaseHSV)
WriteBed <- function(dat.bed, outdir, nwindows, jbait, jsep = "_"){
  if (missing(jbait)){
    baitfix <- ""
  } else {
    baitfix <- paste0("_", jbait)
  }
  # output to name: assay_genotype_nwindows.bed
  assay <- unique(dat.bed$assay)
  genotype <- unique(dat.bed$genotype)
  if (length(assay) != 1 & length(genotype) != 1){
    print("Assay and geno should be UNIQUE")
    return(NULL)
  }
  outfname <- paste0(assay, jsep, genotype, jsep, nwindows, baitfix, ".bed")
  outfname.pval <- paste0(assay, jsep, genotype, jsep, nwindows, baitfix, ".pval.bed")
  trackname <- paste0(assay, jsep, genotype, jsep, nwindows)
  outf <- file.path(outdir, outfname)
  outf.pval <- file.path(outdir, outfname.pval)
  
  outbed <- subset(dat.bed, select = c(chromo, start, end, pval, score.bed))
  outbed$pval <- format(outbed$pval, format="scientific", digits=2)
  outbed$strand <- "+"
  outbed$thickStart <- dat.bed$start
  outbed$thickEnd <- dat.bed$end
  outbed$itemRgb <- dat.bed$rgb
  head(outbed)
  # write header
  
  # write phase amp track
  sink(file = outf, append = FALSE)
  cat("browser position chr6:142371133-142421629")
  cat("\n")
  cat(paste0('track name=AmpPhase-"', trackname, '" description="AmpPhase-', trackname, '" visibility=dense itemRgb="On" useScore=1'))
  cat("\n")
  sink()
  write.table(outbed, file = outf, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  # write again but no itemRgb, just useScore
  sink(file = outf.pval, append = FALSE)
  cat("browser position chr6:142371133-142421629")
  cat("\n")
  cat(paste0('track name=Pval-"', trackname, '" description="Pvalue-', trackname, '" visibility=dense itemRgb="Off" useScore=1'))
  cat("\n")
  sink()
  write.table(outbed, file = outf.pval, append = TRUE, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  print(paste("Output files:", outf, "\n", outf.pval))
  return(data.frame(NULL))
}


HexColorToRgb <- function(hex){
  # HEX input: #FFC8DD
  R <- as.integer(paste0("0x", substr(hex, 2, 3)))
  G <- as.integer(paste0("0x", substr(hex, 4, 5)))
  B <- as.integer(paste0("0x", substr(hex, 6, 7)))
  return(paste(R, G, B, sep = ","))
}

SaturationCurve <- function(x, Vmax, k, x0, expo = 1){
  # exponent for making things curvier
  if (expo == 1){
    warning("Exponential 1 is Sigmoidal function, not Michaelis")
    y <- (Vmax * 2 / (1 + exp(-k * (x - x0)))) - 1
    # y <- Vmax * x^expo / (k^expo + x^expo)
  } else {
    y <- Vmax * x^expo / (k^expo + x^expo)
  }
}

RotatePhase <- Vectorize(function(phase, rotate.hr=0){
  # rotate phase, in hours. If negative add 24. If > 24 minus 24.
  if (is.na(phase)) return(phase)
  phase <- phase + rotate.hr
  if (phase < 0) phase <- phase + 24
  if (phase > 24) phase <- phase - 24
  return(phase)
}, "phase")

StepCurve <- function(x, x.step){
  # Map x to step function from 0 to 1, 1 if x >= x.step
  return(ifelse(x >= x.step, 1, 0))
}

CutoffMap <- function(amp, pval, amp.cutoff, pval.cutoff){
  if (any(is.na(c(amp, pval)))){
    return(0)
  }
  if (amp > amp.cutoff & pval < pval.cutoff){
    return(1)
  } else {
    return(0)
  }
}

HexToHsv <- Vectorize(function(hex){
  R <- as.integer(paste0("0x", substr(hex, 2, 3)))
  G <- as.integer(paste0("0x", substr(hex, 4, 5)))
  B <- as.integer(paste0("0x", substr(hex, 6, 7)))
  rgb2hsv(R, G, B)  
}, "hex")

PhaseAmpPvalToColor <- function(phase, amp, pval, rotate.hr=-8, amp.k = 2, pval.k = 0.25, method = "smooth"){
  # method: "smooth" or "step"
  # if smoothi -> uses SaturationCurve with a power exponent of i. i must be integer
  # rotate.phase: rotate phase by some hours
  phase <- RotatePhase(phase, rotate.hr)
  phase.col <- PhaseToHsv(phase, min.phase = 0, max.phase = 24)
  if (method == "smooth"){
    amp.col <- SaturationCurve(amp, Vmax = 1, k = amp.k, x0 = 0)
    pval.col <- SaturationCurve(-log10(pval), Vmax = 1, k = pval.k, x0 = 0)
  } else if (method == "step"){
    warning("Function not optimized. Some white spots may occur")
    amp.col <- StepCurve(amp, x.step = amp.k)
    pval.col <- StepCurve(-log10(pval), x.step = pval.k)
  } else if (method == "cutoff"){
    amp.col <- mapply(function(amp, pval) CutoffMap(amp, pval, amp.k, pval.k), amp, pval)
    pval.col <- amp.col
  } else {
    # assume it is smooth with a power
    jexpo <- as.integer(strsplit(method, "smooth")[[1]][[2]])
    # amp.col <- SaturationCurve(amp, Vmax = 1, k = amp.k, x0 = 0, expo = jexpo)
    amp.col <- rep(1, length(pval))
    pval.col1 <- SaturationCurve(-log10(pval), Vmax = 1, k = pval.k, x0 = 0, expo = jexpo)
    pval.col2 <- SaturationCurve(amp, Vmax = 1, k = amp.k, x0 = 0, expo = jexpo)
    # take min of pval and amp
    pval.df <- cbind(pval.col1, pval.col2)
    pval.col <- apply(pval.df, 1, function(row) min(row))
  }
  # convert bad phase amp pvals to 0,0,0
  bad.i <- which(is.na(phase.col) | is.na(amp.col) | is.na(pval.col))
  # make it 0 0 0
  if (length(bad.i) > 0){
    warning(paste0("Converting to 0,0,0 for:", paste0(bad.i, collapse = ",")))
    phase.col[bad.i] <- 0; amp.col[bad.i] <- 0; pval.col[bad.i] <- 0  
  }
  hsv.col <- hsv(phase.col, amp.col, pval.col)
  return(hsv.col)
}

