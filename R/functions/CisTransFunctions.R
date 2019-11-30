PutPosInBin <- function(pos){
  # Possible bins
  # nearbin1 (<50kb)
  # nearbin2 (50-100kb)
  # nearbin3 (100-500kb)
  # distal1 (500kb-1000kb)
  # distal2 (1000kb-2000kb)
  # distal3 (2000kb to Inf)
  # trans (trans)
  if (length(pos) > 1){
    pos <- pos[1]
  }
  if (abs(pos) < 50000){
    return("50kb-")
  } else if (abs(pos) > 50000 & abs(pos) < 100000){
    return("50kb-100kb")
  } else if (abs(pos) > 100000 & abs(pos) < 500000){
    return("100kb-500kb")
  } else if (abs(pos) > 500000 & abs(pos) < 1000000){
    return("500kb-1mb")
  } else if (abs(pos) > 1000000 & abs(pos) < 2000000){
    return("1mb-2mb")
  } else if (abs(pos) > 2000000){
    return("2mb+")
  }
}

PutPosInBin.cistrans <- function(pos){
  # Possible bins
  # near (<2000kb)
  # distal3 (2000kb to Inf)
  # trans (trans)
  if (length(pos) > 1){
    pos <- pos[1]
  }
  if (abs(pos) <= 2000000){
    return("<2mb")
  } else if (abs(pos) > 2000000){
    return("2mb+")
  }
}
