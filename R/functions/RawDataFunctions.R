# 2016-02-03
# Handle raw data functions here 
library(hash)

GetSampleNumber <- function(counts.long.raw.sub){
  counts.long.raw.sub <- counts.long.raw.sub %>%
    group_by(bait, time, bkgr, rel.pos) %>%
    mutate(sample.i = seq(length(counts)),
           samp.name = paste(time, bkgr, sample.i, sep = ";"))
  return(counts.long.raw.sub)
}

GetRegion <- function(avg.scores, frag.pos, bait.lst, norm.meth, rmin, rmax, strict.min=0,
                      take.abs.range=TRUE){
  # filter row index based on normalizatoin method, min, max
  # strict.min is a numeric if we want to set non-absolute range but want to
  # restrict from regions too close to bait
  
  if (take.abs.range){
    if (norm.meth == "region"){
      I = which(avg.scores[,1] == bait.lst$chr
                & abs(frag.pos - bait.lst$pos) > strict.min
                & abs(frag.pos - bait.lst$pos) > rmin
                & abs(frag.pos - bait.lst$pos) < rmax)
    } else if (norm.meth == "notregion"){
      I = which(avg.scores[,1] == bait.lst$chr
                & abs(frag.pos - bait.lst$pos) > strict.min
                & abs(frag.pos - bait.lst$pos) > rmax)  # assumes rmin contained in rmax
    } else {
      I = which(avg.scores[,1] == bait.lst$chr  # chromosome of avg score matches bait chromo
                & abs(frag.pos - bait.lst$pos) > strict.min
                & abs(frag.pos - bait.lst$pos) > rmin)
    }
  } else{
    if (norm.meth == "region"){
      I = which(avg.scores[,1] == bait.lst$chr
                & abs(frag.pos - bait.lst$pos) > strict.min
                & (frag.pos - bait.lst$pos) > rmin
                & (frag.pos - bait.lst$pos) < rmax)
    } else if (norm.meth == "notregion"){
      I = which(avg.scores[,1] == bait.lst$chr
                & abs(frag.pos - bait.lst$pos) > strict.min
                & (frag.pos - bait.lst$pos) > rmax)  # assumes rmin contained in rmax
    } else {
      I = which(avg.scores[,1] == bait.lst$chr  # chromosome of avg score matches bait chromo
                & abs(frag.pos - bait.lst$pos) > strict.min
                & (frag.pos - bait.lst$pos) > rmin)
    }
  }
  return(I)
}

GetIndexOfNoisyFrags <- function(mat, jbait, med.cutoff){
  # Given bait remove all fragments that do not have a median of less than med.cutoff across
  # all samples and replicates
  
  # filter columns to look at only one bait
  mat.filt <- mat[, grepl(jbait, colnames(mat))]
  
  # get index of rows where median is less than cutoff
  frags.i <- which(apply(mat.filt, 1, median) > med.cutoff)
  return(frags.i)
}

MakeSampNames <- function(baits, zts, bkgrs, jsep = ";"){
  # Make sampnames from baits, zts, bkgrs
  # e.g., "HOXD4;ZT08;WT"
  # can easily add sample names afterwards
  samp.names <- paste(baits, zts, bkgrs, sep = jsep)
  
  # create hash to track sample numbers
  samp.hash <- hash(unique(samp.names), rep(1, length(unique(samp.names))))
  
  samp.names.withnumb <- samp.names  # init
  # add sample i 
  for (i in seq(length(samp.names))){
    samp.numb <- samp.hash[[samp.names[i]]]
    # update new sampnames
    samp.names.withnumb[i] <- paste(samp.names[[i]], samp.numb, sep = ";")
    # update hash using original sampnames
    samp.hash[[samp.names[i]]] <- samp.hash[[samp.names[[i]]]] + 1
  }
  return(samp.names.withnumb)
}

StrToCond <- function(s){
  # String of form: ZT08;WT;1, get zt, bkgr, samp.i from it
  zt <- strsplit(s, ";")[[1]][[1]]
  bkgr <- strsplit(s, ";")[[1]][[2]]
  samp.i <- strsplit(s, ";")[[1]][[3]]
  return(list(zt = zt, bkgr = bkgr, samp.i = samp.i))
}

GetSampIndex <- function(bait, zt, bkgr, samp.i, baits, zts, bkgrs){
  # get sample index given baits, zts, bkgrs
  samp.names <- MakeSampNames(baits, zts, bkgrs, jsep = ";")
  samp.to.get <- paste(bait, zt, bkgr, samp.i, sep = ";")
  samp.i <- which(samp.names == samp.to.get)
  return(samp.i)
}

RemoveSample <- function(mat, bait, zt, bkgr, samp.i, baits, zts, bkgrs, keep.n.colnames=3){
  # Remove a sample from avg.scores
  # bait: bait you want to remove
  # time for bait you wnat to remove
  # background (KO, WT) fo time, bait you want to remove
  # sample for bait;time;bkgr to remove
  # baits: vector of baits corresponding to colnames in bait
  # zts: vector of times 
  # bkgrs: vector of backgrounds
  # keep.n.colnames: if your mat has chromosome, start, end then keep by setting to 3, otherwise set to 0
  
  # return list:
  # mat.filt, baits.filt, zts.filt, bkgrs.filt: elements filtered out
  # 
  samp.names <- MakeSampNames(baits, zts, bkgrs, jsep = ";")
  
  # grepstr <- paste0("^", bait, ".*;", as.character(samp.i), "$")
  grepstr <- paste(bait, zt, bkgr, samp.i, sep = ";")
  # remove samples and bait
  keep.cols <- rep(TRUE, keep.n.colnames)
  
  mat.filt <- mat[, c(keep.cols, !grepl(grepstr, samp.names))]
  return(mat.filt)
}

FilterSingleBait <- function(bait, avg.scores, baits, zts, bkgrs, metadata.columns = 3){
  # Filter avg scores for a single bait e.g., "HOXD4"
  # need baits, zts, bkgrs to be updated as well
  # allows you to save this all into an object later
  # metadata.columns: index to specify the last column name containing meta data 
  # (e.g., 3 means first 3 cols are metadata eg chromosome, start, end)
  Ibaits <- which(baits == bait)
  # filter baits, zts, bkgrs they contain no metadata cols
  baits <- baits[Ibaits]
  zts <- zts[Ibaits]
  bkgrs <- bkgrs[Ibaits]
  avg.scores <- avg.scores[, c(1:metadata.columns, Ibaits + metadata.columns)]
  return(list(avg.scores = avg.scores, baits = baits, zts = zts, bkgrs = bkgrs))
}