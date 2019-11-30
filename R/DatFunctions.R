SumAroundFragsLong <- function(pos.dat, all.dat, frag.dist = 20000, normalize = TRUE){
  # position of fragment near significant peak
  # all.dat: A.delta not filtered for peaks 
  # frag.dist: how far from fragment you care about
  # normalize: normalize by number of fragments you summed up 
  pos <- unique(pos.dat$pos)
  print(pos.dat)
  if (length(pos) != 1) warning(paste("Expected position length to be 1, found", length(pos)))
  
  dat.sub <- subset(all.dat, abs(pos - top.hit) < frag.dist)
  # sum and normalize by number of bases
  print(pos)
  print(dat.sub)
  integral <- sum(dat.sub$A.delta)
  if (normalize){
    integral <- integral / nrow(dat.sub)
  }
  pos.dat$integral <- integral
  return(pos.dat)
}

SumAroundFrags <- function(top.hit, jbait, all.dat, frag.dist = 20000, normalize = TRUE){
  # position of fragment near significant peak
  # all.dat: A.delta not filtered for peaks 
  # frag.dist: how far from fragment you care about
  # normalize: normalize by number of fragments you summed up 
  dat.sub <- subset(all.dat, bait == jbait & abs(pos - top.hit) < frag.dist)
  # sum and normalize by number of bases
  integral <- sum(dat.sub$A.delta)
  if (normalize){
    integral <- integral / nrow(dat.sub)
  }
  return(integral)
}

SetWidthLimit <- function(frag.pos, jfrag, pos.limit = 3000){
  # frag.pos: where the local min/max is
  # jfrag: where we think the width limit is
  # if jfrag is past pos.limit then set fragment position at position limit (closest to frag.pos)  
  if (abs(frag.pos - jfrag) > abs(frag.pos - sign(frag.pos) * pos.limit)){
    # jfrag is beyond the position limit (closest to frag.pos) therefore we cut to pos limit
    jfrag <- sign(frag.pos) * pos.limit
  }
  return(jfrag)
}


GetPeakWidth <- function(frag.pos, jbait, all.dat, min.delt = 0.05, pos.limit = 3000){
  # get width of peak by finding the distance between fragment to nearest
  # get right limit
  # min.delt: where you want to call beginning of peak
  # pos.limit: if edge fragment is past the pos.limit, set to nearest limit (+/- pos limit)
  # you add pos limit to take into account that if peak is near bait, then you may jump
  # too far from -3000 to 3000
  
  # order data from smallest position to largest, so the fragment closest
  # to frag.pos can naturally be extracted as statistic i from dat.right or dat.left
  all.dat <- all.dat[order(all.dat$pos), ]
  dat.right <- subset(all.dat, bait == jbait & pos - frag.pos > 0)  # [1] is closest to frag.pos
  dat.left <- subset(all.dat, bait == jbait & pos - frag.pos < 0)  # [nrow] is closest to frag.pos
  
  dat.right.filt <- dat.right[which(abs(dat.right$A.delta) < min.delt), ]
  dat.right.frag <- dat.right.filt[1, "pos"]
  dat.right.frag.pos <- SetWidthLimit(frag.pos, dat.right.frag$pos, pos.limit)
    
  dat.left.filt <- dat.left[which(abs(dat.left$A.delta) < min.delt), ]
  dat.left.frag <- dat.left.filt[nrow(dat.left.filt), "pos"] 
  dat.left.frag.pos <- SetWidthLimit(frag.pos, dat.left.frag$pos, pos.limit)
  
  peakwidth <- dat.right.frag.pos - dat.left.frag.pos
  return(peakwidth)
}

GetNextFrag <- function(top.hit, all.dat, direction="right", stat.i = 1){
  # Given peak, get fragment right of peak
  # direction: right | left
  # stat.i, take the ith statistic (can take a vector of values as well, like seq(5) for top 5)
  if (direction == "right"){
    jsub <- subset(all.dat, (pos - top.hit) > 0)
  } else if (direction == "left"){
    jsub <- subset(all.dat, (pos - top.hit) < 0)
  } else {
    warning("Direction must be right or left")
    return(NA)
  }
  # take fragment closest to bait 
  jsub <- jsub[order(abs(jsub$pos - top.hit)), ]
  # take ith row
  return(jsub[stat.i, ])
}

GetNextFragRange <- function(top.hit, all.dat, direction="right", jrange = 20000){
  # Given peak, get fragment right of peak
  # direction: right | left
  # jrange: how far from top.hit (left or right) do you want to take fragments?
  if (direction == "right"){
    jsub <- subset(all.dat, (pos - top.hit) > 0 & abs(pos - top.hit) < jrange)
  } else if (direction == "left"){
    jsub <- subset(all.dat, (pos - top.hit) < 0 & abs(pos - tophit) < jrange)
  } else {
    warning("Direction must be right or left")
    return(NA)
  }
  return(jsub)
}

GetAvgFirstDeriv <- function(top.hit, all.dat, jrange = 20000){
  # Given peak, get second derivative (approximate the step size as average between +/- fragment)
  # alert user if the two step sizes are too large?
  # jrange: left and right of bait how many nucleotides to take average?
  frag.left <- GetNextFragRange(top.hit, all.dat, direction = "left", jrange = jrange)
  frag.right <- GetNextFrag(top.hit, all.dat, direction = "right", jrange = jrange)  
  # TODO: actually avg first derivative might not be the best option  
  warning("Not finished here")
}


LabelSignifHits <- function(dat){
  # check if is.sig group is is.signif or not.signif 
  sig.status <- as.character(unique(dat$is.sig))
  if (length(sig.status) != 1){
    warning("Length should be 1")
  } else {
    if (sig.status == "not.signif"){
      dat$label <- ""
    } else if (sig.status == "is.signif"){
      dat$label <- as.character(seq(nrow(dat)))
    } else {
      warning("Sig status should be not.signif or is.signif")
    }
  }
  return(dat)
}

MergeZscoreAndZscoreKO <- function(dat){
  # sort by genotype
  dat.sorted <- dat[order(dat$genotype, decreasing = TRUE), ]
  
  # check if identical
  nn <- nrow(dat.sorted)
  if (identical(dat.sorted$Zscore[1:(nn/2)], dat.sorted$Zscore[(nn/2 + 1) : nn])){
    dat.sorted$Zscore[(nn/2 + 1) : nn] <- dat.sorted$Zscore.KO[1:(nn/2)]
    dat.sorted$Pstat[(nn/2 + 1) : nn] <- dat.sorted$Pstat.KO[1:(nn/2)]
  } else {
    warning("Not identical, doing nothing...")
    return(dat.sorted)
  }
  return(dat.sorted)
}

IsSignif <- function(genotype, pval, pval.KO, threshold){
  if (genotype == "WT"){
    if (pval < threshold){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if (genotype == "KO"){
    if (pval.KO < threshold){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    warning("Genotype must be WT or KO")
    print(genotype)
  }
}

IsSignif2 <- function(dat){
  # check if length 1 or 2, if length 1 check if it is significant, then report its genotype
  # if length 2, check which is significant, report its genotype, if both report both
  if (nrow(dat) == 1){
    if (dat$is.sig == TRUE){
      is.sig <- dat$genotype[[1]]
    } else {
      is.sig <- "None"
    }
  } else if (nrow(dat) == 2){
    if (length(which(dat$is.sig == TRUE)) == 0){
      is.sig <- "None"
    } else if (length(which(dat$is.sig == TRUE)) == 1){
      is.sig <- dat$genotype[which(dat$is.sig == TRUE)]
    } else if (length(which(dat$is.sig == TRUE)) == 2){
      is.sig <- "KO,WT"
    } else {
      warning("Not yet implemented this case")
    }
  } else {
    warning("is.sig can only be length 1 or 2")
    is.sig <- NA
  }
  return(data.frame(is.sig = is.sig))
}