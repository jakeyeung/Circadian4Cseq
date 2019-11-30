LoadExonIntronCounts <- function(inf){
  # inf <- "/home/shared/4c_seq/bowtie2_WT_Cry1Intron/Liver_Liv06_S6_R1_filtered.counts_mat"
  bname <- basename(inf)
  jtiss <- strsplit(bname, split = "_")[[1]][[1]]
  if (jtiss == "Kidney"){
    # Kidney_ZT8_k_42_filtered.counts_mat
    samp.numb <- as.numeric(gsub("[^\\d]+", "", strsplit(bname, split = "_")[[1]][[4]], perl=TRUE))
    jtime <- as.numeric(gsub("ZT", "", strsplit(bname, split = "_")[[1]][[2]]))
    jgeno.sym <- strsplit(bname, split = "_")[[1]][[3]]  # k or w
    if (jgeno.sym == "k"){
      jgeno <- "Cry1IntronKO"
    } else if (jgeno.sym == "w"){
      jgeno <- "WT"
    } else {
      stop(paste("jgeno must be k or w", jgeno.sym))
    }
  } else if (jtiss == "Liver") {
    # Liver_Liv09_S9_R1_filtered.counts_mat
    samp.numb <- as.numeric(gsub("[^\\d]+", "", strsplit(bname, split = "_")[[1]][[2]], perl=TRUE))
    jtime <- SeriesToTime(samp.numb)
    jgeno <- SeriesToGeno(samp.numb)
  } else {
    stop(paste("jtiss must be Liver or Kidney", jtiss))
  }
  dat <- read.table(inf, header = TRUE, sep = "\t", stringsAsFactors=FALSE)
  dat$samp.numb <- samp.numb
  dat$tissue <- jtiss
  dat$geno <- jgeno
  dat$time <- jtime
  return(dat)
}

PhaseDelta <- function(p.vec){
  if (abs(diff(p.vec)) > 12){
    # take smallest p and add 24
    p.vec[which.min(p.vec)] <- p.vec[which.min(p.vec)] + 24
  }
  # take 2nd minus 1st, as in diff function
  return(diff(p.vec))
}

SeriesToTime <- function(samp.numb){
   samp.numb.by.12 <- (samp.numb - 1) %% 12
  # odds are WT, evens are KO, timepoints 0, 4, 8 ... 
  jtimes <- ifelse(samp.numb %% 2, (floor(samp.numb.by.12/2)) * 4, (floor(samp.numb.by.12/2)) * 4)  # 0,0,4,4,8,8 ... 68,68
}

SeriesToGeno <- function(samp.numb){
  jgenos <- ifelse(samp.numb %% 2 == 1, "WT", "Cry1IntronKO")
}


LoadDataset <- function(jtissue, jeps = 1e-3){
  inf <- paste0("/home/shared/4c_seq/kallisto_WT_Cry1Intron/Mermet_", jtissue, "_Series_1_WT_Cry1IntronKO.txt")
  dat <- read.table(inf, header = TRUE, sep = "\t")
  
  
  
  
  # Ensembl to gene ID ------------------------------------------------------
  
  dat$gene <- Transcript2Gene(dat$target_id, return.original = TRUE)
  
  # If gene name all numbers, then use transcript name
  dat$gene[grepl("^[1-9]*$", dat$gene)] <- dat$target_id[grepl("^[1-9]*$", dat$gene)]
  
  colnames(dat) <- gsub("Liv", "Liver", colnames(dat))
  
  dat.exprs <- dat[, grepl(paste0("^", jtissue),  colnames(dat))]
  sampnames <- colnames(dat)[grepl(paste0("^", jtissue), colnames(dat))]
  
  if (jtissue == "Liver"){
    # colnames need to be readjusted to be "Tissue_Time_Geno"
    samp.numb <- as.numeric(sapply(sampnames, function(s) substr(s, start=6, stop = 7), USE.NAMES = TRUE))
    # odds are WT, evens are KO, timepoints 0, 4, 7 ... 
    
    # times
    samp.numb.by.12 <- (samp.numb - 1) %% 12
    jtimes <- ifelse(samp.numb %% 2, (floor(samp.numb.by.12/2)) * 4, (floor(samp.numb.by.12/2)) * 4)  # 0,0,4,4,8,8 ... 68,68
    jgenos <- ifelse(samp.numb %% 2 == 1, "WT", "Cry1IntronKO")
    jtissues <- rep("Liver", length(jtimes))
    sampnames <- make.names(paste(jtissues, jtimes, jgenos, sep = "_"), unique=TRUE)
  }
  
  
  # fix sampnames if necessary
  if (any(grepl("\\.", sampnames))){
    psplit <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[3]], simplify = TRUE, USE.NAMES = FALSE)
    # set biorep to 1 if no period, otherwise take number after "." and add 1
    biorep <- as.numeric(sapply(psplit, function(p){
      tryCatch(as.numeric(strsplit(p, "\\.")[[1]][[2]]) + 1, error = function(e) 1)
    }, simplify = TRUE, USE.NAMES = FALSE))
    # fix sampnames
    samp.fix.i <- grepl("\\.[1-9]$", sampnames)
    sampnames[samp.fix.i] <- substr(sampnames[samp.fix.i], start = 1, stop = nchar(sampnames[samp.fix.i]) - 2)
  } else {
    biorep <- rep(1, length(sampnames))
  }
  
  
  
  tissues <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[1]], USE.NAMES = FALSE)
  times <- sapply(sampnames, function(s) as.numeric(gsub("ZT", "", strsplit(s, "_")[[1]][[2]])), USE.NAMES = FALSE)
  genos <- sapply(sampnames, function(s) strsplit(s, "_")[[1]][[3]], USE.NAMES = FALSE)
  
  dat.long <- data.frame(gene = dat$gene, 
                         tissue = rep(tissues, each = nrow(dat.exprs)),
                         time = rep(times, each = nrow(dat.exprs)),
                         geno = rep(genos, each = nrow(dat.exprs)),
                         biorep = rep(biorep, each = nrow(dat.exprs)),
                         sampnames = rep(sampnames, each = nrow(dat.exprs)), 
                         exprs = unlist(dat.exprs))
  dat.long$gene <- as.character(dat.long$gene)
  dat.long$geno <- factor(dat.long$geno, levels = c("WT", "Cry1IntronKO"))
  
  dat.long$tissgeno <- paste(dat.long$tissue, dat.long$geno, sep = "_")
  dat.long$tissue <- dat.long$tissgeno
  
  # fix time to plot over entire series
  dat.long$time <- mapply(function(jtime, jbiorep) jtime + (jbiorep - 1) * 24, dat.long$time, dat.long$biorep)
  
  dat.long <- dat.long %>%
    group_by(gene, tissue, time, geno) %>% 
    summarise(exprs = sum(exprs))
  
  # jeps <- 1e-3
  jcutoff <- -2
  ggplot(dat.long, aes(x = log2(exprs + jeps))) + geom_density() + geom_vline(xintercept = jcutoff)
  
  # log2 and tpm separate
  dat.long$tpm <- dat.long$exprs
  dat.long$exprs <- log2(dat.long$exprs + jeps)
  dat.long$tissue.orig <- dat.long$tissue
  dat.long$tissue <- paste(dat.long$tissue.orig, dat.long$geno, sep = "_")
  return(dat.long)
}

PlotRnaseq <- function(dat.long.merged.sub, single.day=FALSE, time.space = 8, jtitle="", linear.eps = FALSE){
  if (linear.eps != FALSE){
    if (is.na(as.numeric(linear.eps))) stop("linear.eps must be numeric if not FALSE")
    dat.long.merged.sub$exprs <- dat.long.merged.sub$tpm
    jylab <- "TPM"
  } else {
    jylab <- "Log2(TPM)"
  }
  if (jtitle == ""){
    jtitle <- dat.long.merged.sub$gene[[1]]
  }
  if (single.day){
    source("/home/yeung/projects/tissue-specificity/scripts/functions/PlotGeneAcrossTissues.R")
    dat.long.merged.sub <- dat.long.merged.sub %>% 
      ungroup() %>% 
      mutate(experiment = "RNASeq")
    dat.long.merged.sub <- ConvertToSingleDay(dat.long.merged.sub)
  }
  m <- ggplot(dat.long.merged.sub, aes(x = time, y = exprs, linetype = geno, colour = tissue)) + geom_point() + geom_line() +
    theme_bw(24) + xlab("ZT") + ylab(jylab) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle) + 
    scale_x_continuous(breaks = seq(0, max(dat.long.merged.sub$time), time.space))
  if (single.day){
    m <- m + geom_errorbar(data = dat.long.merged.sub, aes(ymin=exprs-sem, ymax=exprs+sem, linetype = geno, colour = tissue), linetype = "solid", size=0.5, width=0.5) 
  }
  return(m)
}
