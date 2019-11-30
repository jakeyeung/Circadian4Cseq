GetBaitLocations <- function(datmain, quiet=FALSE){
  # given directory containing N_*.RData, retrieve a data frame of bait names and their locations
  if (missing(datmain)){
    #     datmain <- "/home/yeung/data/4c_seq/region_normalized_pseudo_sigma_1000000/sigma_2500.psuedo_300"
    datmain <- "/home/yeung/data/4c_seq/bait_locations"
  }
  datfiles <- list.files(datmain)
  baitfiles <- sapply(datfiles, function(f){
    if (substring(f, first = 1, last = 2) == "N_"){
      return(f)
    } else {
      return(NULL)
    }
  })
  baitfiles <- unlist(baitfiles, use.names = F)
  # Load up each baitfile and concatenate into a dataframe
  bait.long <- lapply(baitfiles, function(f){
    
    if (!quiet) print(paste("Loading: ", f))
    load(file.path(datmain, f))  # contains object bait
    return(data.frame(bait = bait$name, chromo = bait$chr, pos = bait$pos))
  })
  bait.long <- do.call(rbind, bait.long)
  return(bait.long)
}

GetBaitLocations2 <- function(inf){
  # use proper bait locations from metadata
  if (missing(inf)){
    inf <- "/home/shared/4c_seq/metadata/batch1.batch2.baits.bed"
  }
  bait.locs <- read.table(inf, header=FALSE, sep = "\t", col.names = c("bait", "chromo", "start", "end"))
  bait.locs <- bait.locs %>%
    mutate(pos = (start + end) / 2) %>%
    arrange(bait)
  return(subset(bait.locs, select = c(bait, chromo, pos)))
}

Vectorize(GetGeneFromBaitName <- function(jgene, striptext = "auto", get.upper = "auto"){
  # Cry1intron -> CRY1 (if get.upper == TRUE)
  # striptext can be auto, then it searches for intron, long, short as striptext
  if (striptext == "auto" & get.upper == "auto"){
    if (endsWith(jgene, suffix = "intron")){
      striptext <- "intron"
      get.upper <- TRUE
    } else if (endsWith(jgene, suffix = "long")){
      striptext <- "long"
      get.upper <- FALSE
    } else if (endsWith(jgene, suffix = "short")){
      striptext <- "short"
      get.upper <- FALSE
    } else {
      warning(paste("unknowng striptext for gene", jgene))
      return(jgene)
    }
  }
  jgene <- strsplit(jgene, striptext)[[1]][[1]]
  if (get.upper == TRUE){
    jgene <- toupper(jgene)
  }
  return(jgene)
}, vectorize.args = "jgene")