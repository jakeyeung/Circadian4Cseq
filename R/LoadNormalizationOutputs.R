AddLR <- function(pos){
  if (pos < 0){
    return("Left")
  } else {
    return("Right")
  }  
}

GetDeltaSubLst <- function(counts.long, adj.pval.cutoff = 0.05){
  source("scripts/functions/FindPeaks.R")
  # reorder factors
  counts.long$genotype <- factor(as.character(counts.long$genotype), levels = c("WT", "KO"))
  
  start <- Sys.time()
  counts.delt <- counts.long %>%
    group_by(pseudo, sig) %>%
    do(GetDelta(., filter.peaks = FALSE))
  
  # adj.pval.cutoff <- 0.1
  #   adj.pval.cutoff <- 0.05
  
  counts.delt.sub <- counts.delt %>%
    group_by(bait, sig, pseudo, genotype, tissue) %>%
    do(FilterForPeaks(.))
  
  counts.delt.sub <- counts.delt.sub %>%
    group_by(bait, genotype, tissue, sig, pseudo) %>%
    mutate(pval.adj = p.adjust(pval),
           pval.adj.KO = p.adjust(pval.KO),
           pval.adj.kid = p.adjust(pval.kid))
  
  counts.delt.sub$is.sig <- as.factor(sapply(counts.delt.sub$pval.adj, IsSignif, adj.pval.cutoff))
  counts.delt.sub$is.sig.KO <- as.factor(sapply(counts.delt.sub$pval.adj.KO, IsSignif, adj.pval.cutoff))
  counts.delt.sub$is.sig.kid <- as.factor(sapply(counts.delt.sub$pval.adj.kid, IsSignif, adj.pval.cutoff))
  
  counts.delt$pval.row <- apply(counts.delt, 1, function(row) AssignPvalToCondition(row))
  counts.delt$zscore.row <- apply(counts.delt, 1, function(row) AssignZscoreToCondition(row))
  counts.delt.sub$pval.row <- apply(counts.delt.sub, 1, function(row) AssignPvalToCondition(row))
  counts.delt.sub$signif.row <- apply(counts.delt.sub, 1, function(row) AssignSignifToCondition(row))
  
  return(list(counts.delt = counts.delt, counts.delt.sub = counts.delt.sub, counts.long = counts.long))
}

IsSignif <- function(p, adj.pval.cutoff){
    if (is.na(p)) return("not.signif")
    if (p < adj.pval.cutoff){
      return("is.signif")
    } else {
      return("not.signif")
    }
}

AssignZscoreToCondition <- function(row, genotype.i=3, tissue.i=4, zscore.i=10, zscore.ko.i=11, zscore.kid.i=12){
  if (row[genotype.i] == "WT" & row[tissue.i] == "Liver"){
    zscore.row = as.numeric(row[zscore.i])
  } else if (row[genotype.i] == "KO" & row[tissue.i] == "Liver"){
    zscore.row = as.numeric(row[zscore.ko.i])
  } else if (row[genotype.i] == "WT" & row[tissue.i] == "Kidney"){
    zscore.row = as.numeric(row[zscore.kid.i])
  } else {
    stop("genotype must be KO or WT, tissue must be Liver or Kidney")
  }
  return(zscore.row)
}

AssignZscoreToConditionLK <- function(row, genotype.i=3, time.i=4, zscore.08.i=10, zscore.20.i=12){
  if (row[genotype.i] == "WT" & row[time.i] == "ZT08"){
    zscore.row = as.numeric(row[zscore.08.i])
  } else if (row[genotype.i] == "WT" & row[time.i] == "ZT20"){
    zscore.row = as.numeric(row[zscore.20.i])
  } else {
    stop("genotype must be KO or WT, tissue must be Liver or Kidney")
  }
  return(zscore.row)
}

AssignPvalToCondition <- function(row, genotype.i=3, tissue.i=4, pval.i=13, pval.ko.i=14, pval.kid.i=15){
  if (row[genotype.i] == "WT" & row[tissue.i] == "Liver"){
    pval.row = as.numeric(row[pval.i])
  } else if (row[genotype.i] == "KO" & row[tissue.i] == "Liver"){
    pval.row = as.numeric(row[pval.ko.i])
  } else if (row[genotype.i] == "WT" & row[tissue.i] == "Kidney"){
    pval.row = as.numeric(row[pval.kid.i])
  } else {
    stop("genotype must be KO or WT, tissue must be Liver or Kidney")
  }
  return(pval.row)
}

AssignPvalToConditionLK <- function(row, genotype.i=3, time.i=4, pval.08.i=11, pval.20.i=13){
  if (row[genotype.i] == "WT" & row[time.i] == "ZT08"){
    pval.row = as.numeric(row[pval.08.i])
  } else if (row[genotype.i] == "WT" & row[time.i] == "ZT20"){
    pval.row = as.numeric(row[pval.20.i])
  } else {
    stop("genotype must be KO or WT, tissue must be Liver or Kidney")
  }
  return(pval.row)
}

AssignSignifToCondition <- function(row, genotype.i=3, tissue.i=4, signif.i=19, signif.ko.i=20, signif.kid.i=21){
  if (row[genotype.i] == "WT" & row[tissue.i] == "Liver"){
    signif.row = row[signif.i]
  } else if (row[genotype.i] == "KO" & row[tissue.i] == "Liver"){
    signif.row = row[signif.ko.i]
  } else if (row[genotype.i] == "WT" & row[tissue.i] == "Kidney"){
    signif.row = row[signif.kid.i]
  } else {
    stop("genotype must be KO or WT, tissue must be Liver or Kidney")
  }
  return(signif.row)
}

GetDeltaLK <- function(counts.long){
  # Get Delta for Liver and Kidney comparison
  counts.delt <- counts.long %>%
    group_by(pos, bait, genotype, time, sig, pseudo) %>%
    summarise(A.mean = mean(A), A.var = mean(Avar), A.delta = diff(A), 
              ZscoreLK08 = unique(Zscore.lk08), pvalLK08 = unique(Pstat.lk08), ZscoreLK20 = unique(Zscore.lk20), pvalLK20 = unique(Pstat.lk20))  # ZT20 - ZT08 WT
  counts.delt$pval.row <- apply(counts.delt, 1, function(row) AssignPvalToConditionLK(row))
  counts.delt$zscore.row <- apply(counts.delt, 1, function(row) AssignZscoreToConditionLK(row))
  return(counts.delt)
}

GetDelta <- function(counts.long, filter.peaks = FALSE, adj.pval.cutoff = 0.1){
  # Get delta from counts.long, assumes counts.long is properly subsetted
  # filter.peaks: returns delt which are filtered by peaks 
  counts.delt <- counts.long %>%
    group_by(pos, bait, genotype, tissue, sig, pseudo) %>%
    summarise(A.mean = mean(A), A.var = mean(Avar), A.delta = diff(A), Zscore = unique(Zscore), Zscore.KO = unique(Zscore.KO), Zscore.kid = unique(Zscore.kid), 
              pval = unique(Pstat), pval.KO = unique(Pstat.KO), pval.kid = unique(Pstat.kid))  # ZT20 - ZT08 WT
  if (filter.peaks){
    counts.delt <- counts.delt %>%
      group_by(bait, sig, pseudo, genotype, tissue) %>%
      do(FilterForPeaks(.))
    counts.delt <- counts.delt %>%
      group_by(bait, sig, pseudo, genotype, tissue) %>%
      mutate(pval.adj = p.adjust(pval, method = "BH"),
             pval.adj.KO = p.adjust(pval.KO, method = "BH"),
             pval.adj.kid = p.adjust(pval.kid, method = "BH"))
    counts.delt$is.sig <- as.factor(sapply(counts.delt$pval.adj, function(p){
      if (is.na(p)) return("not.signif")
      if (p < adj.pval.cutoff){
        return("is.signif")
      } else {
        return("not.signif")
      }
    }))
    counts.delt$is.sig.KO <- as.factor(sapply(counts.delt$pval.adj.KO, function(p){
      # need to handle NAs
      if (is.na(p)) return("not.signif")
      if (p < adj.pval.cutoff){
        return("is.signif")
      } else {
        return("not.signif")
      }
    }))
    counts.delt$is.sig.kid <- as.factor(sapply(counts.delt$pval.adj.kid, function(p){
      # need to handle NAs
      if (is.na(p)) return("not.signif")  
      if (p < adj.pval.cutoff){
        return("is.signif")
      } else {
        return("not.signif")
      }
    }))
  }
  return(counts.delt)
}

GetDeltaSub <- function(counts.delt, adj.pval.cutoff = 0.1){
  out.delt <- counts.delt %>%
    group_by(bait, sig, pseudo, genotype, tissue) %>%
    do(FilterForPeaks(.)) %>%
    mutate(pval.adj = p.adjust(pval, method = "BH"),
           pval.adj.KO = p.adjust(pval.KO, method = "BH"),
           pval.adj.kid = p.adjust(pval.kid, method = "BH"))
#   out.delt <- out.delt %>%
#     group_by(bait, sig, pseudo, genotype, tissue) %>%

  out.delt$is.sig <- as.factor(sapply(out.delt$pval.adj, function(p){
    if (is.na(p)) return("not.signif")
    if (p < adj.pval.cutoff){
      return("is.signif")
    } else {
      return("not.signif")
    }
  }))
  out.delt$is.sig.KO <- as.factor(sapply(out.delt$pval.adj.KO, function(p){
    # need to handle NAs
    if (is.na(p)) return("not.signif")
    if (p < adj.pval.cutoff){
      return("is.signif")
    } else {
      return("not.signif")
    }
  }))
  out.delt$is.sig.kid <- as.factor(sapply(out.delt$pval.adj.kid, function(p){
    # need to handle NAs
    if (is.na(p)) return("not.signif")  
    if (p < adj.pval.cutoff){
      return("is.signif")
    } else {
      return("not.signif")
    }
  }))
  return(out.delt)
}

MergeFixBatch1and2 <- function(counts.long.batch1, counts.long.batch2, add.label = NA){
  # Merge counts batch1 and batch2 together from Jerome
  # batch1 is from old
  # batch2 is new batch
  
  # label batch1 and batch2 just in case
  counts.long.batch1$batch <- 1
  counts.long.batch2$batch <- 2
  counts.long <- rbind(counts.long.batch1, counts.long.batch2)
  rm(counts.long.batch1, counts.long.batch2)
  
  if (!is.na(add.label)){
    counts.long$label <- add.label
  }
  
  # rename genotype Liver to WT (I messed up with labels in my batch 2)
  counts.long$genotype <- sapply(as.character(counts.long$genotype), function(g){
    if (g == "Liver"){
      return("WT")
    } else {
      return(g)
    }
  })
  
  # rearrange factor so it is WT KO, Liver Kidney
  counts.long$genotype <- factor(as.character(counts.long$genotype), levels = c("WT", "KO"))
  
  # rename Gys2 to GYS2
  counts.long$bait <- as.factor(sapply(as.character(counts.long$bait), function(b){
    if (b == "Gys2"){
      return("GYS2") 
    } else {
      return(b)
    }
  }))
  
  # rename Pi3kap1 -> Pik3ap1
  counts.long$bait <- as.factor(sapply(as.character(counts.long$bait), function(b){
    if (b == "Pi3kap1"){
      return("Pik3ap1") 
    } else {
      return(b)
    }
  }))
  return(counts.long)
}

AddLogPosLR <- function(counts.sub, lr.only = FALSE, make.neg = FALSE, pos.name = "pos", out.name = "pos.log"){
  counts.sub$lr <- sapply(counts.sub[[pos.name]], function(p){
    if (p > 0){
      return("right")
    } else {
      return("left")
    }
  })
  
  if (lr.only) return(counts.sub)
  
  counts.sub[[out.name]] <- sapply(counts.sub[[pos.name]], function(p){
    if (p > 0){
      return(log10(p))
    } else if (p < 0){
      if (make.neg == FALSE) return(log10(-p))
      if (make.neg == TRUE) return(-log10(-p))
    } else {
      warning("p shouldnt be 0")
    }
  })
  return(counts.sub)
}

GetExcludedRange <- function(pos, jrange){
  # jrange if in cis contains a "gap" which should be region taht is 
  # too close to bait and was excluded
  # extract this exlcuded range as a genomic pos
  full.range.i <- min(jrange):max(jrange)
  excluded.range.i <- full.range.i[which(!full.range %in% range)]
  excluded.pos <- pos[excluded.range.i]
  return(excluded.pos)
}

GetSampTime <- function(jname, remove.prefix=FALSE){
  # jname: ZT08 WT -> 8 or ZT08 depending on remove.prefix
  time <- strsplit(jname, split = " ")[[1]][[1]]  # ZT08
  if (remove.prefix){
    time <- strsplit(time, split = "ZT")[[1]][[2]]
  }
  return(time)
}

GetGenotype <- function(jname, handle.tissue = TRUE){
  # jname: ZT08 WT -> WT
  # exception: if Kidney then return WT
  genotype <- strsplit(jname, split = " ")[[1]][[2]]  # WT or KO or Kidney
  if (genotype == "Kidney"){
    genotype <- "WT"
  }
  if (genotype != "WT" & genotype != "KO"){
    warning(paste("Genotype neither WT or KO", genotype))
  }
  return(genotype)
}

GetTissue <- function(jname, default.tissue = "Liver"){
  # if tissue extracted matches does not match either "Kidney" or "Liver", then set to default.tissue
  tissue <- strsplit(jname, split = " ")[[1]][[2]]  # Kidney, WT or KO
  if (tissue != "Liver" & tissue != "Kidney"){
    tissue <- default.tissue
  }
  return(tissue)
}

MakeLong <- function(pos, jrange, pstats, zstats, A, Avar, bait.lst, relative=TRUE, pstats.KO = NA, zstats.KO = NA, jchromo = NULL){
  # After loading .RData (output from norm_fit_all.R), put it
  # into a long data format
  # 
  # pos (from clog list): 
  #   position which will be used to translate the
  #   jrange into a real genomic position
  # jrange (numeric):
  #   index of positions of the genomic positions considered 
  #    (need pos[jrange[i]] to get real genomic position)
  # pstats:
  #  list of 4 (only first one matters)
  #  contains p-values from statistics
  # zstats:
  #   lke pstats, but with z-scores
  # A:
  #   like pstats, but contains list of 4 should correspond to the 4 samples (ZTs WT-KO). This is the "position effect"
  # Avar:
  #   like pstats, but conatins list of 4 should correspond to the 4 samples (ZTs WT-KO). This is the variance of position effect
  # bait.lst:
  #  contains useful things like name, chr, pos
  # relative:
  #  make positions relative to bait position
  # pstats.KO and zstats.KO
  #  pstats and zstats for KO is optional
  
  
  pos.range <- pos[jrange]  # position within range
  
  if (relative){
    pos.range <- pos.range - bait.lst$pos
  }
  
  # do for all
  N <- length(A[[1]])
  nsamp <- length(names(bait.lst$cols))
  A.all <- unlist(A)
  Avar.all <- unlist(Avar)
  times <- sapply(names(bait.lst$cols), GetSampTime)
  # handle kidney
  genotypes <- sapply(names(bait.lst$cols), GetGenotype)
  tissues <- sapply(names(bait.lst$cols), GetTissue)
  # genotype_time <- sapply(names(bait.lst$cols), function(s) gsub(" ", "_", s))
  genotype_time = paste(tissues, genotypes, times, sep = "_")
  
  # handle chromosomes: if NULL (unspecified), use bait chromo, otherwise use current.chromo specified in .RData
  if (is.null(jchromo)){
    # take bait region if chromo not specifieid 
    chromo.vec <- rep(bait.lst$chr, N * nsamp)
  } else {
    chromo.vec <- rep(jchromo, N * nsamp)
  }
  
  if (is.null(pstats.KO) & is.null(zstats.KO)){
    print("Getting WT pstat and zstat only...")
    dlong <- data.frame(bait = rep(bait.lst$name, N * nsamp), 
                        genotype_time = rep(genotype_time, each = N), 
                        tissue = rep(tissues, each = N), genotype = rep(genotypes, each = N), 
                        time = rep(times, each = N),
                        chromo = chromo.vec, 
                        pos = rep(pos.range, nsamp), 
                        A = A.all, 
                        Avar = Avar.all, 
                        Zscore = rep(zstats, nsamp), 
                        Pstat = rep(pstats, nsamp))
  } else {
    print("Getting WT and KO pstat and zstat...")
#     print(length(zstats))
#     print(length(zstats.KO))
    dlong <- data.frame(bait = rep(bait.lst$name, N * nsamp), 
                        genotype_time = rep(genotype_time, each = N), 
                        tissue = rep(tissues, each = N), genotype = rep(genotypes, each = N), 
                        time = rep(times, each = N),
                        chromo = chromo.vec, 
                        pos = rep(pos.range, nsamp), 
                        A = A.all, 
                        Avar = Avar.all, 
                        Zscore = rep(zstats, nsamp), 
                        Pstat = rep(pstats, nsamp),
                        Zscore.KO = rep(zstats.KO, nsamp),
                        Pstat.KO = rep(pstats.KO, nsamp))
  }
  return(dlong)
}

GetFnamePairs <- function(datdir){
  # Search data dir for pairs of .RData (F_ and N_ for each sample)
  # return list of pairs of filenames we can just load and feed into MakeLong
  files.all <- list.files(datdir)
  samps <- sapply(files.all, function(s) strsplit(s, "_")[[1]][[3]])
  
  fname.pairs <- list()
  for (samp in samps){
    grepstr <- paste(samp, "_", sep = "")
    fname.pair <- files.all[grepl(grepstr, files.all)]
    if (length(fname.pair) != 2) warning("Unexpected number of files. Should be 2")
    # with dir
    fname.pair <- sapply(fname.pair, function(f) file.path(datdir, f))
    fname.pairs[[samp]] <- fname.pair
  }
  return(fname.pairs)
}

MakeLongAllSamps <- function(datdir, sig, pseudo){
  # create counts.long from .RData F_ and N_ samplesA
  # sig and pseudo are the sigmas and pseudocounts used in the analysis
  jenv <- new.env()
  fname.pairs <- GetFnamePairs(datdir)
  
  counts.long <- data.frame()
  for (fpair in fname.pairs){
    # need to load both N_ and F_ of a sample, then run MakeLong
    d.temp <- GetDatLongFromFpair(fpair, jenv)
    #     for (f in fpair){
    #       print(paste("Loading file:", f))
    #       load(f, verbose=F, envir = jenv) 
    #     }
    #     #  print(ls(envir = jenv))
    #     d.temp <- MakeLong(pos = jenv$clog$pos, jrange = jenv$range, pstats = jenv$OUTPSTAT[[1]], zstats = jenv$OUTZSTAT[[1]], A = jenv$OUT, Avar = jenv$OUTVAR, bait.lst = jenv$bait)
    counts.long <- rbind(counts.long, d.temp)
  }
  counts.long$sig = sig
  counts.long$pseudo = pseudo
  return(counts.long)
}

GetDatLongFromFpair <- function(fpair, jenv){
  # .GlobalEnv is global environment for jenv
  for (f in fpair){
    print(paste("Loading file:", f))
    load(f, verbose=F, envir = jenv) 
  }
  #  print(ls(envir = jenv))
  # chec, if OUTPSTAT and OUTZSTAT are null, if null, then take only 1st column, otherwise take 2 columns
  if (is.null(jenv$OUTPSTAT[[2]]) & is.null(jenv$OUTZSTAT[[2]])){
    dlong <- MakeLong(pos = jenv$clog$pos, 
                      jrange = jenv$range, 
                      pstats = jenv$OUTPSTAT[[1]], 
                      zstats = jenv$OUTZSTAT[[1]], 
                      A = jenv$OUT, 
                      Avar = jenv$OUTVAR, 
                      bait.lst = jenv$bait,
                      chromo = jenv$current.chr)
  } else if (!is.null(jenv$OUTPSTAT[[2]]) & !is.null(jenv$OUTZSTAT[[2]])){
    dlong <- MakeLong(pos = jenv$clog$pos, 
                      jrange = jenv$range, 
                      pstats = jenv$OUTPSTAT[[1]], 
                      zstats = jenv$OUTZSTAT[[1]], 
                      A = jenv$OUT, 
                      Avar = jenv$OUTVAR, 
                      bait.lst = jenv$bait,
                      pstats.KO = jenv$OUTPSTAT[[2]],
                      zstats.KO = jenv$OUTZSTAT[[2]],
                      chromo = jenv$current.chr)
  } else {
    warning("OUTPSTAT and OUTZSTAT neither both null or both not null...")
  }
  return(dlong)
}