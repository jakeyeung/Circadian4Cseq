GetSeriesFromSampi <- function(s){
  s.numb <- as.numeric(gsub("S", "", s))
  series.numb <- (floor((s.numb - 1) / 12)) + 1
  return(series.numb)
}

CumSumTrans <- function(jsub){
  jsub <- subset(jsub, cistrans == "trans") %>% arrange(start)
  jsub$counts.cumsum <- cumsum(jsub$counts)
  jsub$counts.cumsum.norm <- jsub$counts.cumsum / sum(jsub$counts)
  return(jsub)
}

FitLinear <- function(jdat, slope.only=FALSE){
  jfit <- lm(counts.cumsum ~ rel.pos, jdat)
  m <- jfit$coefficients[["rel.pos"]]
  if (is.na(m)){
    pval <- NA
  } else {
    pval <- summary(jfit)$coefficients["rel.pos", "Pr(>|t|)"]
  }
  if (slope.only){
    return(m)
  } else{
    return(data.frame(m = m, pval = pval))
  }
}

SlidingWindowLinReg <- function(dat, spots.start, spots.end, jwin.step = FALSE, pos.cname = "jpos"){
  # http://stats.stackexchange.com/questions/3051/mean-of-a-sliding-window-in-r
  # result <- vector(length = length(spots.start))
  result <- rep(NA, length(spots.start))
  
  # prev.pos <- dat[[pos.cname]][1]  # init
  for(i in 1:length(spots.start)){
    curr.pos <- dat[[pos.cname]][i]
    if (is.na(spots.start[i]) | is.na(spots.end[i])) next
    
    # if (jwin.step != FALSE & i > 1){
    #   # skip current fragment unless frag dist is larger than jwin.step away from previous frag dist
    #   if (abs(curr.pos) - abs(prev.pos) < jwin.step){
    #     next
    #   }
    # }
    result[i] <- FitLinear(dat[spots.start[i]:spots.end[i], ], slope.only=TRUE)
    # # update prev.pos
    # prev.pos <- curr.pos
  }
  return(result)
}


FitExpDecay <- function(dat, pseudo=1e-5){
  # handle slopes of zeros with a pseudocount
  fits <- lm(log10(slopes.roll + pseudo) ~ jpos, dat)
  yint <- coef(fits)[[1]]
  slope <- coef(fits)[[2]]
  return(data.frame(yint = yint, slope = slope))
}

GetSlopes.trans <- function(dat.sub, jstep = 1e6, nfrags.min = 10, jwin.step = 1e6){
  # dat.sub, single bait in trans
  # lin.step how many steps by distance
  jposs <- dat.sub$rel.pos
  spots.start <- rep(NA, length(jposs))
  spots.end <- rep(NA, length(jposs))
  
  for (i in 1:length(jposs)){
    spots.start[i] <- i
    # only consider psots that are within jwin.step away
    curr.pos <- jposs[i]
    if (i > 1){
      if (abs(curr.pos) - abs(prev.pos) < jwin.step){
        next
      }
    }
    # print(spots.start[i])
    # do max of nfrags.min fragments or jposs + jstep
    fragsi.jstep <- which(abs(jposs) > abs(jposs[i]) & abs(jposs) < abs(jposs[i]) + jstep)
    if (length(fragsi.jstep) > nfrags.min){
      spots.end[i] <- max(fragsi.jstep)
    } else {
      spots.end[i] <- min((i + nfrags.min), length(jposs))  # dont go past length of jposs??
    }
    prev.pos <- curr.pos
  }
  
  max.i <- which.max(spots.end)
  dat.sub <- dat.sub[1:max.i, ]
  dat.sub.fits.rolling <- SlidingWindowLinReg(dat.sub, spots.start[1:max.i], spots.end[1:max.i], jwin.step = jwin.step, pos.cname = "rel.pos")
  dat.sub$slopes.roll <- abs(dat.sub.fits.rolling)  # take absolute
  return(dat.sub)
}

GetSlopes <- function(dat.sub, jstep = 0.2, nfrags.min = 10, clip.frags = 100, jwin.step = FALSE, .log = TRUE){
  # dat.sub: single bait, single sample, condition, left OR right. .
  # jstep: smooth over what distance?
  # nfrags.min: minimum number of frags over which to smooth
  # clip.frags: clip how many frags at the end? 
  
  # check signs match
  jfactor <- unique(sign(dat.sub$rel.pos))
  if (length(jfactor) != 1){
    warning("Expect all relative positions to be same sign.")
  }
  
  if (.log){
    if (any(dat.sub$rel.pos < 0)){
      dat.sub$jpos <- -log10(-dat.sub$rel.pos)
    } else {
      dat.sub$jpos <- log10(dat.sub$rel.pos)
    }
  } else {
    dat.sub$jpos <- dat.sub$rel.pos
  }
  
  # BEGIN: get spots corresponding to indexes start and stops 
  # for linear regression across a dynamic sliding window
  # take next 10 fragments or jpos of 0.5, whichever is larger 
  # create "spots.start and sports.end"
  jposs <- dat.sub$jpos
  spots.start <- rep(NA, length(jposs))
  spots.end <- rep(NA, length(jposs))
  
  # jstep <- 0.2
  # nfrags.min <- 10
  for (i in 1:length(jposs)){
    spots.start[i] <- i
    
    curr.pos <- jposs[i]
    if (i > 1){
      if (abs(curr.pos) - abs(prev.pos) < jwin.step){
        next
      }
    }
    
    # do max of nfrags.min fragments or jposs + jstep
    fragsi.jstep <- which(abs(jposs) > abs(jposs[i]) & abs(jposs) < abs(jposs[i]) + jstep)
    if (length(fragsi.jstep) > nfrags.min){
      spots.end[i] <- max(fragsi.jstep)
    } else {
      spots.end[i] <- min((i + nfrags.min), length(jposs))  # dont go past length of jposs??
    }
    prev.pos <- curr.pos
  }
  
  max.i <- which.max(spots.end)
  dat.sub <- dat.sub[1:max.i, ]
  dat.sub.fits.rolling <- SlidingWindowLinReg(dat.sub, spots.start[1:max.i], spots.end[1:max.i], jwin.step = jwin.step)
  dat.sub$slopes.roll <- abs(dat.sub.fits.rolling)  # take absolute
  
  # remove last x fragments if last x fragments
  # clip.frags <- 100
  max.i.slopes <- max.i - clip.frags
  dat.sub <- dat.sub[1:max.i.slopes, ]
  # add log10pos if log, do nothing in linear
  if (.log){
    dat.sub$log10pos <- dat.sub$jpos
  }
  return(dat.sub)
}

RemoveFragments <- function(jsub, frags.to.remove = 5){
  warning("Potential bug where frag.i = 0 is not defined. ")
  # remove n fragments left and right of bait.
  jsub.left <- subset(jsub, rel.pos < 0)
  jsub.right <- subset(jsub, rel.pos > 0)
  left.pos <- jsub.left[order(jsub.left$rel.pos, decreasing = TRUE), ]$rel.pos[frags.to.remove]
  right.pos <- jsub.right[order(jsub.right$rel.pos, decreasing=FALSE), ]$rel.pos[frags.to.remove]
  # remove frags (less than or greater than, no equal)
  jsub <- subset(jsub, rel.pos < left.pos | rel.pos > right.pos)
  return(jsub)
}

CumSumLeftRight <- function(jsub){
  jsub$LR <- sapply(jsub$rel.pos, function(pos) ifelse(pos < 0, "left", "right"))
  i <- 1
  jsub.lst <- list()
  for (left.right in c("left", "right")){
    jjsub <- subset(jsub, LR == left.right)
    jjsub <- jjsub[order(abs(jjsub$rel.pos), decreasing = FALSE), ]
    jjsub$counts.cumsum <- cumsum(jjsub$counts)
    jjsub$counts.cumsum.norm <- jjsub$counts.cumsum / sum(jjsub$counts)
    i <- i + 1
    jsub.lst[[i]] <- jjsub
  }
  jsub2 <- do.call(rbind, jsub.lst)
  return(jsub2)
}

# used in vitalit

OrderAroundZero <- function(x){
  # c(-10, -9, ,,, 100, 101) -> c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)'  
  # zero means closest fragment, ties go by earliest index
  # x <- seq(-10, by = 2, length.out = 20)
  if (is.unsorted(x)){
    warning("Vector x must be sorted, returning NA")
    return(NA)
  }
  x.right <- sort.int(x[x >= 0], decreasing=FALSE, index.return = TRUE)$ix
  x.left <- sort.int(x[x < 0], decreasing=TRUE, index.return = TRUE)$ix  # adjust sign at the end
  # set closest fragment to be 0, so +1 and -1 means 1 fragment from closest fragment
  if (x[which.min(abs(x))] >= 0){
    # shift x.rights so it begins with 0
    x.right <- x.right - 1
  } else if (x[which.min(abs(x))] < 0){
    x.left <- x.left - 1
  }
  x.i <- c(-1 * x.left, x.right)
  return(x.i)
}

OrderCisFragments <- function(x, cistrans){
  if (all(cistrans == "trans")){
    # return(rep("trans", length(x)))
    return(rep(Inf, length(x)))
  } else if (all(cistrans == "cis")){
    return(OrderAroundZero(x))
  } else{
    warning(paste0("cistrans must be cis or trans: ", cistrans))
  }
}

GetBadSamples <- function(inmain, bait, prefix="removed_samps.", suffix=".txt"){
  if (missing(inmain)){
    inmain <- "/home/shared/4c_seq/removed_samps"
  }
  fname <- paste0(prefix, bait, suffix)
  samps <- readLines(file.path(inmain, fname))
  return(samps)
}
