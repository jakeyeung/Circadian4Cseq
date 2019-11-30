FindPeaks <- function(vec,bw=1,x.coo=c(1:length(vec)))
{
  # https://rtricks.wordpress.com/2009/05/03/an-algorithm-to-find-local-extrema-in-a-vector/
  # bw controls how noisy you tolerate your signal
  pos.x.max <- NULL
  pos.y.max <- NULL
  pos.x.min <- NULL
  pos.y.min <- NULL 	
  for(i in 1:(length(vec)-1)) 	{ 		if((i+1+bw)>length(vec)){
    sup.stop <- length(vec)}else{sup.stop <- i+1+bw
    }
    if((i-bw)<1){inf.stop <- 1}else{inf.stop <- i-bw}
    subset.sup <- vec[(i+1):sup.stop]
    subset.inf <- vec[inf.stop:(i-1)]
    
    is.max   <- sum(subset.inf > vec[i]) == 0
    is.nomin <- sum(subset.sup > vec[i]) == 0
    
    no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
    no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
    
    if(is.max & is.nomin){
      pos.x.max <- c(pos.x.max,x.coo[i])
      pos.y.max <- c(pos.y.max,vec[i])
    }
    if(no.max & no.nomin){
      pos.x.min <- c(pos.x.min,x.coo[i])
      pos.y.min <- c(pos.y.min,vec[i])
    }
  }
  return(list(pos.x.max, pos.y.max, pos.x.min, pos.y.min))
}

FilterForPeaks <- function(dat){
  # filter for peaks given a counts.long subsetted data frame
  peaks <- FindPeaks(dat$A.delta, bw = 1, x.coo = dat$pos)
  peaks.pos <- c(peaks[[1]], # max
                 peaks[[3]])  # min
  return(subset(dat, pos %in% peaks.pos))
}