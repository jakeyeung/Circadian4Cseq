# Jake Yeung
# BedGraphFunctions.R
#  
# 2019-07-15

GetRightLeftLimits <- function(pos, end.length = 0, get.right=TRUE){
  # pos are positions that should be sorted!
  # end.length: the right and left ends add an arbitrary distance
  pos <- sort(pos)
  pos.left <- rep(NA, length(pos))
  pos.right <- rep(NA, length(pos))
  for (i in seq(length(pos))){
    p <- pos[i]
    p1r <- pos[i + 1]  # 1 to the right
    p1l <- pos[i - 1]  # 1 to the left

    if (is.na(p1r)) p1r <- p + end.length  # if index is beyond length, it returns NA
    if (length(p1l) == 0) p1l <- p - end.length  # if index is 0, it returns numeric(0)

    pos.right[i] <- (p + p1r) / 2
    pos.left[i] <- (p + p1l) / 2
  }
  if (get.right){
    return(round(pos.right))
  } else {
    return(round(pos.left))
  }
}

WriteToFile <- function(outfile, chromo, start, end, value){
  sink(file = outfile)
  for (i in seq(length(chromo))){
    chr <- chromo[i]
    s <- start[i]
    e <- end[i]
    v <- value[i]
    cat(paste0(chr, "\t", s, "\t", e, "\t", v, "\n"))
  }
  sink()
}

WriteBedGraph <- function(chromo, start, end, value, outfile, overwrite=TRUE){

  if (length(chromo) != length(start)) stop("Unequal lengths")
  if (length(chromo) != length(end)) stop("Unequal lengths")
  if (length(chromo) != length(value)) stop("Unequal lengths")

  if (!file.exists(outfile)){
    WriteToFile(outfile, chromo, start, end, value)
  } else {
    if (overwrite){
      print("File exists, overwriting...")
      WriteToFile(outfile, chromo, start, end, value)
    } else {
      warning("File exists, not overwriting")
    }
  }
  return(NULL)
}

