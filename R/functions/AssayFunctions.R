AssayToBed <- function(fits.zscore, outdir, nwindows, mybait, jamp.k, jpval.k, jmethod){
  fits.zscore <- AnnotateDat(fits.zscore, cname = "label")
  fits.zscore$nwindows <- sapply(as.character(fits.zscore$ints), function(s) length(strsplit(s, ",")[[1]]))
  print("Fits zscore:")
  print(head(fits.zscore))
  fits.zscore <- subset(fits.zscore, select = c(-ints, -windows))
  
  # Get bed files from rhythmicity ------------------------------------------
  
  fits.zscore$phase <- sapply(fits.zscore$phase, HandleNAs, set.val = 0)
  fits.zscore$amp <- sapply(fits.zscore$amp, HandleNAs, set.val = 0)
  fits.zscore$pval <- sapply(fits.zscore$pval, HandleNAs, set.val = 1)
  
  max.score <- max(-log10(fits.zscore$pval))  # -log10pval
  score.bed <- LinearInterpolate(min(-log10(fits.zscore$pval)), max.score, 0, 1000, x = -log10(fits.zscore$pval))
  score.bed <- sapply(score.bed, function(s){
    if (is.na(s)){
      return(0)
    }
    if (s > 1000){
      s <- 1000
    } 
    return(as.integer(s))
  })
  fits.zscore$score.bed <- score.bed
  
  fits.zscore$hex.col <- PhaseAmpPvalToColor(fits.zscore$phase, fits.zscore$amp, fits.zscore$pval, rotate.hr = -8, amp.k = jamp.k, pval.k = jpval.k, method = jmethod)
  fits.zscore$rgb <- sapply(fits.zscore$hex.col, HexColorToRgb)
  
  fits.zscore %>%
    group_by(assay, genotype) %>%
    do(WriteBed(., outdir, nwindows, mybait))
}

LoadAssayZscores <- function(robjpath = "/home/yeung/projects/4c_seq/data/Robjs/cyclix/all_baits/fits.zscore.jenv.nwin.6.Robj"){
  load(robjpath, v=T)
  fits.zscore <- subset(fits.zscore, select = c(-windows, -ints))
  
  fits.zscore$bait <- sapply(fits.zscore$label, function(l) strsplit(l, ";")[[1]][[4]])
  fits.zscore$genotype <- sapply(fits.zscore$label, function(l) strsplit(l, ";")[[1]][[3]])
  fits.zscore$assay <- sapply(fits.zscore$label, function(l) strsplit(l, ";")[[1]][[2]])
  fits.zscore$chromo <- sapply(as.character(fits.zscore$region), function(l) strsplit(l, ":")[[1]][[1]])
  fits.zscore$start <- sapply(as.character(fits.zscore$region), function(l) strsplit(strsplit(l, ":")[[1]][[2]], "-")[[1]][[1]])
  fits.zscore$end <- sapply(as.character(fits.zscore$region), function(l) strsplit(strsplit(l, ":")[[1]][[2]], "-")[[1]][[2]])
  
  return(fits.zscore)
}

Annotate4CSeqWithCyclix <- function(counts.delt.subsub, fits.zscore, bait.locs, jassay, amp.k = 2, pval.k = 0.25, mapping.method = "smooth"){
  if (missing(fits.zscore)){
    print("Fits zscore not found, loading from functions...")
    fits.zscore <- LoadAssayZscores()
  }
  
  # if (USE.FRAGMENTS){
  #   # counts.delt.subsub <- subset(counts.delt, bait %in% jbaits & tissue == "Liver")
  #   # counts.delt.subsub <- counts.delt
  # } else {
  #   # max/min peaks only , counts.delt.sub should not be null 
  #   # counts.delt.subsub <- subset(counts.delt.sub, bait %in% jbaits & tissue == "Liver")
  # }
  
  baits.locs.hash <- hash(as.character(bait.locs$bait), bait.locs$pos)
  counts.delt.subsub$frags.abs <- apply(counts.delt.subsub, 1, function(row){
    pos <- as.numeric(row[1]); bait <- as.character(row[2])
    return(pos + baits.locs.hash[[bait]])
  })
  
  # annotate to region
  fits.zscore.sub <- subset(fits.zscore, assay == jassay)
  counts.delt.subsub <- counts.delt.subsub %>%
    group_by(bait, genotype, tissue) %>%
    mutate(regions = GetCyclixRegions(frags.abs, genotype, fits.zscore = fits.zscore.sub, nwindow = 6, winsize = 500))
  # head(counts.delt.subsub$regions)
  
  # annotate
  # get H3K27ac signal from region
  pval.hash <- hash(as.character(fits.zscore.sub$label), as.numeric(fits.zscore.sub$pval))
  amp.hash <- hash(as.character(fits.zscore.sub$label), as.numeric(fits.zscore.sub$amp))
  phase.hash <- hash(as.character(fits.zscore.sub$label), as.numeric(fits.zscore.sub$phase))
  
  # add pval amp phase to counts.delt.sub
  counts.delt.subsub[[paste0("pval.", jassay)]] <- sapply(as.character(counts.delt.subsub[["regions"]]), function(r){
    pval <- pval.hash[[r]]
    if (is.null(pval)){
      return(NA)
    }
    return(pval)
  })
  counts.delt.subsub[[paste0("amp.", jassay)]] <- sapply(as.character(counts.delt.subsub[["regions"]]), function(r){
    amp <- amp.hash[[r]]
    if (is.null(amp)){
      return(NA)
    }
    return(amp)
  }) 
  counts.delt.subsub[[paste0("phase.", jassay)]] <- sapply(as.character(counts.delt.subsub[["regions"]]), function(r){
    phase <- phase.hash[[r]]
    if (is.null(phase)){
      return(NA)
    } else {
      return(phase)
    }
  })
  
  # plot result
  counts.delt.subsub$colour <- PhaseAmpPvalToColor(counts.delt.subsub[[paste0("phase.", jassay)]], 
                                                   counts.delt.subsub[[paste0("amp.", jassay)]], 
                                                   counts.delt.subsub[[paste0("pval.", jassay)]], 
                                                   rotate.hr = -8, 
                                                   amp.k = amp.k,
                                                   pval.k = pval.k,
                                                   method = mapping.method)
  counts.delt.subsub$alpha <- SaturationCurve(-log10(counts.delt.subsub[[paste0("pval.", jassay)]]), Vmax = 1, k = pval.k, x0 = 0)
  return(counts.delt.subsub)
}


GetCyclixRegions <- function(frags.abs, genotype, fits.zscore, nwindow, winsize){
  geno <- unique(genotype)
  if (length(geno) != 1) stop("Geno must be length 1")
  fits.zscore.sub <- subset(fits.zscore, genotype == geno)
  
  # assign frags.abs to a region
  starts <- as.numeric(fits.zscore.sub$start)
  ends <- as.numeric(fits.zscore.sub$end)
  # nwindow <- 6  # 3 to left, 3 to right
  # winsize <- 500
  starts <- starts + winsize * nwindow / 2
  ends <- ends - winsize * nwindow / 2
  regions <- as.character(fits.zscore.sub$label)
  
  frags.abs.regions <- sapply(frags.abs, function(f){
    i <- which(f > starts & f < ends)
    if (length(i) > 1){
      warning(paste0(f, " assigned to first element of i:", paste0(regions[i], collapse = ","), "\n"))
      return(regions[i][[1]])  # take first but warn
    } else if (length(i) == 0){
      return(NA)  # no match
    }
    return(regions[i])
  }, USE.NAMES = FALSE)
  # print(frags.abs.regions)
  # Replace 0s with NAs
  # frags.abs.regions[which(sapply(frags.abs.regions, function(f) length) == 0)] <- NA
  return(frags.abs.regions)
}

AddCyclixTo4CSeq <- function(counts.delt.subsub, jbaits, fits.zscore, jassay, baits.locs){
  baits.locs.hash <- hash(as.character(bait.locs$bait), bait.locs$pos)
  counts.delt.subsub$frags.abs <- apply(counts.delt.subsub, 1, function(row){
    pos <- as.numeric(row[1]); bait <- as.character(row[2])
    return(pos + baits.locs.hash[[bait]])
  })
  
  # annotate to region
  fits.zscore.sub <- subset(fits.zscore, assay == jassay)
  counts.delt.subsub <- counts.delt.subsub %>%
    group_by(bait, genotype, tissue) %>%
    mutate(regions = GetCyclixRegions(frags.abs, genotype, fits.zscore = fits.zscore.sub, nwindow = 6, winsize = 500))
  head(counts.delt.subsub$regions)
  
  # annotate
  # get H3K27ac signal from region
  pval.hash <- hash(as.character(fits.zscore.sub$label), as.numeric(fits.zscore.sub$pval))
  amp.hash <- hash(as.character(fits.zscore.sub$label), as.numeric(fits.zscore.sub$amp))
  phase.hash <- hash(as.character(fits.zscore.sub$label), as.numeric(fits.zscore.sub$phase))
  
  # add pval amp phase to counts.delt.sub
  counts.delt.subsub[[paste0("pval.", jassay)]] <- sapply(as.character(counts.delt.subsub[["regions"]]), function(r){
    pval <- pval.hash[[r]]
    if (is.null(pval)){
      return(NA)
    }
    return(pval)
  })
  counts.delt.subsub[[paste0("amp.", jassay)]] <- sapply(as.character(counts.delt.subsub[["regions"]]), function(r){
    amp <- amp.hash[[r]]
    if (is.null(amp)){
      return(NA)
    }
    return(amp)
  }) 
  counts.delt.subsub[[paste0("phase.", jassay)]] <- sapply(as.character(counts.delt.subsub[["regions"]]), function(r){
    phase <- phase.hash[[r]]
    if (is.null(phase)){
      return(NA)
    } else {
      return(phase)
    }
  })  
  counts.delt.subsub$colour <- PhaseAmpPvalToColor(counts.delt.subsub[[paste0("phase.", jassay)]], 
                                                   counts.delt.subsub[[paste0("amp.", jassay)]], 
                                                   counts.delt.subsub[[paste0("pval.", jassay)]])
  counts.delt.subsub$alpha <- SaturationCurve(-log10(counts.delt.subsub[[paste0("pval.", jassay)]]), Vmax = 1, k = 0.25, x0 = 0)
  return(counts.delt.subsub)
}

HandleNAs <- function(x, set.val = 0){
  # convert NAs to set.val
  if (is.na(x)){
    x <- set.val
  }
  return(x)
}

GetRegionFromLab <- function(lab, bed, baits.long, jrange, bait.range, nwindows, nstep, scale.signal){
  # get assay, geno, bait from lab: chr6:142441602-142443602;DNAse;WT;GYS2
  jassay <- strsplit(lab, ";")[[1]][[2]]
  jgeno <- strsplit(lab, ";")[[1]][[3]]
  jbait <- strsplit(lab, ";")[[1]][[4]]
  return(GetRegionsMatchBait(bed, baits.long, jassay, jgeno, jbait, jrange, bait.range, nwindows, nstep, scale.signal = jscale))
}

DiffPhase <- function(phase1, phase2, period=24){
  # Get phase1 and phase2 difference, modulo period
  # phase1 is reference
  jdiff <- abs(diff(c(phase1, phase2)))
  jdiff.min <- min(period - jdiff, jdiff)
  return(jdiff.min)
}

GetRegions <- function(bed, jchromo, jstart, jend, jassay, jgeno, nwindows, nstep, jbait = FALSE, scale.signal = FALSE, bind.bait = FALSE){
  bed.sub <- subset(bed, chromo == jchromo & start >= jstart & end <= jend & assay == jassay & genotype == jgeno)
  if (nrow(bed.sub) == 0) return(NULL)
  if (bind.bait){
    bed.bait <- subset(bed, chromo == jchromo & start >= (baits.sub$pos - bait.range) & end <= (baits.sub$pos + bait.range) & genotype == "WT" & assay == jassay)
    bed.bait$region.type <- "bait"
  }
  
  windows <- sort(as.character(unique(bed.sub$window)))
  if (nstep == 1){
    # in nstep 1, we do a "sliding window", windows must be odd number
    # windows must be odd number
    if (nwindows %% 2 == 1){
      nwindows <- nwindows + 1
      print(paste0("Changing window to even: New window=", nwindows))
    }
  }
  
  i <- 1
  
  jlst <- expandingList()
  # init
  while (TRUE){
    window.left.i <- i
    window.right.i <- i + nwindows
    if (window.right.i > length(windows)) break
    window.vec <- windows[window.left.i:window.right.i]
    # get window.label
    # print(window.vec)
    i <- i + nstep
    bed.window <- subset(bed.sub, window %in% window.vec)
    if (nstep > 1){
      region <- paste0(bed.window$chromo[[1]], ":", min(bed.window$start), "-", max(bed.window$end))
    } else if (nstep == 1){
      # take middle window as region
      if (length(window.vec) %% 2 == 1){
        region <- window.vec[(length(window.vec) + 1) / 2]  # assumes window.vec is odd
      } else {
        print(window.vec)
        print(length(window.vec))
        stop("window length must be odd")
      }
    }
    if (jbait != FALSE){
      region <- paste0(region, ";", jassay, ";", jgeno, ";", jbait)
    } else {
      region <- paste0(region, ";", jassay, ";", jgeno)
    }
    bed.window$region.type <- "distal"
    if (bind.bait){
      bed.baitwin <- rbind(bed.bait, bed.window)
    } else {
      bed.baitwin <- bed.window
    }
    bed.baitwin$region.label <- region
    # scale to handle different amplitudes being assigned different models
    if (scale.signal){
      bed.baitwin <- bed.baitwin %>%
        group_by(region.type, window) %>%
        mutate(signal.log10 = scale(signal.log10, center=TRUE, scale=TRUE))
      if (any(is.nan(range(bed.baitwin$signal.log10)))){
        bed.baitwin <- data.frame(NULL)
      }
      # jlst.flat$signal.log10 <- scale(jlst.flat$signal.log10, center = TRUE, scale = TRUE)
    }
    
    jlst$add(bed.baitwin)
  }
  jlst.flat <- jlst$as.list()
  jlst.flat <- do.call(rbind, jlst.flat)
  return(jlst.flat)
}

GetRegionsMatchBait <- function(bed, baits.long, jassay, jgeno, jbait, jrange, bait.range, nwindows, nstep, scale.signal = FALSE, bind.bait = TRUE){
  # Output dataframe of region of nwindows and its bait region.
  # This df can be used to fit rhythmic parameters.
  # Push these into an environment and it should be ready to go
  baits.sub <- subset(baits.long, bait == jbait)
  jchromo <- as.character(baits.sub$chromo[[1]])
  jstart <- baits.sub$pos - jrange
  jend <- baits.sub$pos + jrange
  
  jlst.flat <- GetRegions(bed, jchromo, jstart, jend, jassay, jgeno, nwindows, nstep, jbait = jbait, scale.signal = scale.signal, bind.bait = bind.bait)
  # bed.sub <- subset(bed, chromo == jchromo & start >= jstart & end <= jend & assay == jassay & genotype == jgeno)
  # if (nrow(bed.sub) == 0) return(NULL)
  # if (bind.bait){
  #   bed.bait <- subset(bed, chromo == jchromo & start >= (baits.sub$pos - bait.range) & end <= (baits.sub$pos + bait.range) & genotype == "WT" & assay == jassay)
  #   bed.bait$region.type <- "bait"
  # }
  # 
  # windows <- sort(as.character(unique(bed.sub$window)))
  # if (nstep == 1){
  #   # in nstep 1, we do a "sliding window", windows must be odd number
  #   # windows must be odd number
  #   if (nwindows %% 2 == 1){
  #     nwindows <- nwindows + 1
  #     print(paste0("Changing window to even: New window=", nwindows))
  #   }
  # }
  # 
  # i <- 1
  # 
  # jlst <- expandingList()
  # # init
  # while (TRUE){
  #   window.left.i <- i
  #   window.right.i <- i + nwindows
  #   if (window.right.i > length(windows)) break
  #   window.vec <- windows[window.left.i:window.right.i]
  #   # get window.label
  #   # print(window.vec)
  #   i <- i + nstep
  #   bed.window <- subset(bed.sub, window %in% window.vec)
  #   if (nstep > 1){
  #     region <- paste0(bed.window$chromo[[1]], ":", min(bed.window$start), "-", max(bed.window$end))
  #   } else if (nstep == 1){
  #     # take middle window as region
  #     if (length(window.vec) %% 2 == 1){
  #       region <- window.vec[(length(window.vec) + 1) / 2]  # assumes window.vec is odd
  #     } else {
  #       print(window.vec)
  #       print(length(window.vec))
  #       stop("window length must be odd")
  #     }
  #   }
  #   region <- paste0(region, ";", jassay, ";", jgeno, ";", jbait)
  #   bed.window$region.type <- "distal"
  #   if (bind.bait){
  #     bed.baitwin <- rbind(bed.bait, bed.window)
  #   } else {
  #     bed.baitwin <- bed.window
  #   }
  #   bed.baitwin$region.label <- region
  #   # scale to handle different amplitudes being assigned different models
  #   if (scale.signal){
  #     bed.baitwin <- bed.baitwin %>%
  #       group_by(region.type, window) %>%
  #       mutate(signal.log10 = scale(signal.log10, center=TRUE, scale=TRUE))
  #     if (any(is.nan(range(bed.baitwin$signal.log10)))){
  #       bed.baitwin <- data.frame(NULL)
  #     }
  #     # jlst.flat$signal.log10 <- scale(jlst.flat$signal.log10, center = TRUE, scale = TRUE)
  #   }
  #   
  #   jlst$add(bed.baitwin)
  # }
  # jlst.flat <- jlst$as.list()
  # jlst.flat <- do.call(rbind, jlst.flat)
  return(jlst.flat)
}

AnnotateDat <- function(fits.all.long, cname = "gene"){
  fits.all.long$chromo <- sapply(as.character(fits.all.long[[cname]]), function(g){
    return(strsplit(g, ":")[[1]][[1]])
  })
  fits.all.long$start <- sapply(as.character(fits.all.long[[cname]]), function(g){
    coord <- strsplit(g, ";")[[1]][[1]]
    startend <- strsplit(coord, ":")[[1]][[2]]
    start <- strsplit(startend, "-")[[1]][[1]]
    return(as.numeric(start))
  })
  fits.all.long$end <- sapply(as.character(fits.all.long[[cname]]), function(g){
    coord <- strsplit(g, ";")[[1]][[1]]
    startend <- strsplit(coord, ":")[[1]][[2]]
    end <- strsplit(startend, "-")[[1]][[2]]
    return(as.numeric(end))
  })
  fits.all.long$center <- mapply(function(s, e) mean(c(s, e)), fits.all.long$start, fits.all.long$end)
  
  # assays and geno
  fits.all.long$assay <- sapply(as.character(fits.all.long[[cname]]), function(g){
    return(strsplit(g, ";")[[1]][[2]])
  })
  fits.all.long$genotype <- sapply(as.character(fits.all.long[[cname]]), function(g){
    return(strsplit(g, ";")[[1]][[3]])
  })
  # bait
  fits.all.long$bait <- sapply(as.character(fits.all.long[[cname]]), function(g){
    return(strsplit(g, ";")[[1]][[4]])
  })
  return(fits.all.long)
}