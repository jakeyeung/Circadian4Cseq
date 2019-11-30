LoadRobjPlotSummary <- function(inf, outdir, sigma=2500, pseudo=500){
  # inf: Robj merged, run from vitalit
  
  source("scripts/functions/PlotFunctions.R")
  source("scripts/functions/FindPeaks.R")
  source("scripts/functions/RawDataFunctions.R")
  source("scripts/functions/LoadNormalizationOutputs.R")
  source("scripts/functions/PlotGeneAcrossTissues.R")
  
  load("Robjs/dat.long.fixed_rik_genes.Robj")
  outf <- strsplit(basename(inf), "\\.")[[1]][[1]]
  print(paste("Outfile:", outf))
  pseudo.low <- 500
  sigma <- 2500
  
  dir.create(outdir, showWarnings = FALSE)
  load(inf, v=T)
  
  counts.long <- counts.long.merged; rm(counts.long.merged)
  
  counts.deltsub.lst <- GetDeltaSubLst(counts.long, 0.05)
  counts.delt <- counts.deltsub.lst$counts.delt; counts.delt.sub <- counts.deltsub.lst$counts.delt.sub; counts.long <- counts.deltsub.lst$counts.long; rm(counts.deltsub.lst)
  # filter close to bait
  # counts.long <- subset(counts.long, abs(pos) > 8000)
  
  pseudo.high <- max(unique(counts.long$pseudo))
  pseudo.low <- min(unique(counts.long$pseudo))
  mindist <- 0.5e6
  maxdist <- 0
  jbaits <- as.character(unique(counts.long$bait))
  
  
  # compare Liver and Kidney
  baits.lk <- as.character(unique(subset(counts.long, tissue == "Kidney")$bait))
  baits.ko <- as.character(unique(subset(counts.long, genotype == "KO")$bait))
  
  counts.delt.lk <- subset(counts.long, bait %in% baits.lk & genotype != "KO") %>%
    group_by(pseudo) %>%
    do(GetDeltaLK(.))
  
  # jsub.delt <- subset(counts.delt, bait==jbait, pos <= mindist)
  # add LR
  
  counts.long$LR <- sapply(counts.long$pos, AddLR)
  counts.delt$LR <- sapply(counts.delt$pos, AddLR)
  counts.delt.lk$LR <- sapply(counts.delt.lk$pos, AddLR)
  
  # set up the volcanos
  counts.delt.sub <- counts.delt %>%
    group_by(bait, sig, pseudo, genotype, tissue) %>%
    do(FilterForPeaks(.))
  counts.delt.sub <- counts.delt.sub %>%
    group_by(bait, sig, pseudo, genotype, tissue) %>%
    mutate(pval.adj = p.adjust(pval.row, method = "BH"))
  
  counts.delt.lk.sub <- counts.delt.lk %>%
    group_by(bait, sig, pseudo, genotype, time) %>%
    do(FilterForPeaks(.))
  counts.delt.lk.sub <- counts.delt.lk.sub %>%
    group_by(bait, sig, pseudo, genotype, time) %>%
    mutate(pval.adj = p.adjust(pval.row, method = "BH"))
  
  pdf(paste0(outdir, "/", outf, ".", pseudo.low, ".sigma.", sigma, ".pdf"))
  for (jbait in jbaits){
    print(paste("Printing bait:", jbait))
    # PlotSignalDeltaLR(jbait, counts.long, counts.delt, pseudo.low, mindist, add.peaks = counts.delt.sub)
    PlotSignalDeltaLR(jbait, counts.long, counts.delt, pseudo.low, mindist, add.peaks = FALSE)
    PlotZscorePvalLR(jbait, counts.delt, pseudo.low, mindist)
    v <- PlotVolcano(jbait, counts.delt.sub, jxlab = "Liver log10 delta signal (ZT20 - ZT08)", by.tissue = TRUE)
    print(v)
    if (jbait %in% baits.lk){
      m1 <- PlotSignalLivVsKidLR(jbait, subset(counts.long, genotype != "KO"), pseudo.low, mindist)
      m2 <- PlotPvalLKLR(jbait, counts.delt.lk, pseudo.low, mindist)
      #     m2 <- PlotDeltaLKLR(jbait, counts.delt.lk, pseudo.low, mindist, add.peaks = counts.delt.lk.sub)  # test it worsk
      multiplot(m1, m2, col=1)
      v.lk <- PlotVolcano(jbait, counts.delt.lk.sub, jxlab = "WT log10 delta signal (Kidney - Liver)", by.tissue = FALSE)
      print(v.lk)
    }
  }
  dev.off()
  return(list("counts.long"=counts.long, "counts.delt"=counts.delt, "counts.delt.lk"=counts.delt.lk, "counts.delt.sub"=counts.delt.sub, "counts.delt.lk.sub"=counts.delt.lk.sub))
}