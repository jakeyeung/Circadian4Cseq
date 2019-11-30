LoadBed <- function(inf, jassay, jtime, jgeno){
  print(paste("Loading", inf))
  bed <- read.table(inf, header=FALSE)
  bed$assay <- jassay
  bed$time <- jtime
  bed$genotype <- jgeno
  return(bed)
}

SplitPosIntoRegions <- function(pos.vec, distmax){
  # pos.vec must be ordered
  pos.starts <- c()
  pos.ends <- c()
  pos.starts <- c(pos.starts, pos.vec[1])
  
  for (i in 2:length(pos.vec)){
    pos.prev <- pos.vec[i - 1]
    pos.now <- pos.vec[i]
    if ((pos.now - pos.prev) > distmax){
      pos.ends <- c(pos.ends, pos.prev)
      pos.starts <- c(pos.starts, pos.now)
    }
  }
  if ((length(pos.starts) - length(pos.ends)) == 1){
    # add pos.ends
    pos.ends <- c(pos.ends, pos.prev)
  }
  return(data.frame(starts=pos.starts, ends=pos.ends))
}


GetRegions <- function(counts.hits.sub, distmax){
  if (nrow(counts.hits.sub) < 2) return(data.frame(NULL))
  #   jbait <- as.character(counts.hits.sub$bait[[1]])
  regions.df <- SplitPosIntoRegions(counts.hits.sub$pos, distmax)
  return(regions.df)
}


GetAdeltSigDelt <- function(jsub, bed.27ac, counts.delt, minsize=2000){
  if (jsub$region.size < minsize) return(data.frame(NULL))
  bed.sub <- subset(bed.k27ac, chromo == jsub$chromo & start > jsub$starts.abs & end < jsub$ends.abs & genotype == "WT")
  delt.sub <- subset(counts.delt, bait == as.character(jsub$bait) & pos > jsub$starts & pos < jsub$ends & tissue == "Liver" & genotype == "WT") %>%
    group_by(bait, genotype, tissue) %>%
    summarise(A.delta = mean(A.delta))
  
  bed.zt08 <- subset(bed.sub, time %in% c("ZT06", "ZT10")) %>%
    group_by(assay, genotype) %>%
    summarise(signal = mean(signal)) %>%
    mutate(time = "ZT08")
  
  bed.zt20 <- subset(bed.sub, time %in% c("ZT18", "ZT22")) %>%
    group_by(assay, genotype) %>%
    summarise(signal = mean(signal)) %>%
    mutate(time = "ZT20")
  
  # delt.both <- data.frame(bait = as.character(jsub$bait), genotkpe = "WT", tissue = "Liver", sig.delta = bed.zt20$signal - bed.zt08$signal, A.delta = delt.sub$A.delta)
  delt.both <- data.frame(sig.delta = bed.zt20$signal - bed.zt08$signal, A.delta = delt.sub$A.delta)
  return(delt.both)
}
