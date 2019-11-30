FilterTADs <- function(bait.sub, dat.tads, dist.cutoff = 1000000){
  dat.sub <- subset(dat.tads, chromo == as.character(bait.sub$chromo[[1]]) &
                      (abs(dat.tads$start - bait.sub$pos) < dist.cutoff | 
                         abs(dat.tads$end - bait.sub$pos) < dist.cutoff))
  # min(abs(start - bait.sub$pos), abs(end - bait.sub$pos)) < dist.cutoff)
  #                       abs(center - bait.sub$pos) < dist.cutoff)
  # label with bait name
  # print(abs(dat.tads$start - bait.sub$pos) < dist.cutoff | abs(dat.tads$end - bait.sub$pos) < dist.cutoff)
  if(nrow(dat.sub) == 0)
    return(data.frame())
  dat.sub$bait <- bait.sub$bait
  # add distance relative 
  dat.sub$start.rel <- dat.sub$start - bait.sub$pos
  dat.sub$end.rel <- dat.sub$end - bait.sub$pos
  return(dat.sub)
}
