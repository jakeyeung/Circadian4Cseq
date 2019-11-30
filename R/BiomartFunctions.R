# functions from UTR scripts
library(methods)
library(IRanges)

CapitalizeFirstLetter <- function(gene){
  # CRY1 -> Cry1
  return(paste(toupper(substr(gene, 1, 1)), tolower(substr(gene, 2, nchar(gene))), sep = ""))
}

Gene2StartEnd <- function(gene.list, return.original=TRUE, mm9=TRUE) {
  # make gene.list begin with Upper, end with Lower
  gene.list <- CapitalizeFirstLetter(gene.list)	
  # library("biomaRt")
  # https://support.bioconductor.org/p/74322/  need to use host and biomart
  if (mm9){
    jhost="may2012.archive.ensembl.org"
    gene.list <- Gene2Ensembl(gene.list, return.original=TRUE)
    gene.attr <- "ensembl_gene_id"
    gene.filtr <- "ensembl_gene_id"
  } else {
    jhost="www.ensembl.org"
    gene.attr <- "external_gene_name"
    gene.filtr <- gene.attr
  }
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host=jhost)
    gos <- getBM(gene.list,
                 attributes=c(gene.attr, "start_position", "end_position"),
                 filters=c(gene.filtr),
                 mart=mart.obj)
    gl <- gos[match(gene.list, gos[,1]), c(2, 3)]
      ## if not found, then keep ENSEMBL name
      print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
      if (return.original){
	          gl[is.na(gl)] <- gene.list[is.na(gl)]
        }
        return(gl)
}

setMethod("$", "GRanges", function(x, name) { # {{{
  # http://grokbase.com/t/r/bioconductor/122mwahdhw/bioc-add-extra-columns-to-granges-metadata
  elementMetadata(x)[, name]
}) # }}}

geneRanges <- 
  function(db, column="ENTREZID")
  {
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementLengths(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }

splitColumnByOverlap <-
  function(query, subject, column="ENTREZID", ...)
  {
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
  }

geneRanges <- 
  function(db, column="ENTREZID")
  {
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementLengths(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }

splitColumnByOverlap <-
  function(query, subject, column="ENTREZID", ...)
  {
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
  }

genomicRangesToBed <- function(gr, gene.name.col = "gene", strand.col = "strand"){
  bed <- data.frame(seqnames=gr$seqnames,
                    starts=gr$start-1,
                    ends=gr$end,
                    genename=unlist(gr[[gene.name.col]]),
                    strand=gr[[strand.col]])
  return(bed)
}

Vectorize(AnnotateFromHash <- function(x, hash.tbl){
  x <- as.character(x)
  out <- hash.tbl[[x]]
  if (is.null(out)){
    out <- NA
  }
  return(out)
}, vectorize.args = "x")


# Functions from promoter and enhancers scripts -----------------------------------------



Transcript2Gene <- function(gene.list, return.original=TRUE) {
  library("biomaRt")
  # https://support.bioconductor.org/p/74322/  need to use host and biomart
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
  gos <- getBM(gene.list, attributes=c("ensembl_transcript_id", "external_gene_name"),
               filters=c("ensembl_transcript_id"),
               mart=mart.obj)
  gl <- gos[match(gene.list, gos[,1]), 2]
  ## if not found, then keep ENSEMBL name
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original){
    gl[is.na(gl)] <- gene.list[is.na(gl)] 
  }
  return(gl)
}

Transcript2Strand <- function(gene.list, return.original=TRUE) {
  library("biomaRt")
  # https://support.bioconductor.org/p/74322/  need to use host and biomart
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
  gos <- getBM(gene.list, attributes=c("ensembl_transcript_id", "strand"),
               filters=c("ensembl_transcript_id"),
               mart=mart.obj)
  gl <- gos[match(gene.list, gos[,1]), 2]
  ## if not found, then keep ENSEMBL name
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original){
    gl[is.na(gl)] <- gene.list[is.na(gl)] 
  }
  return(gl)
}

Gene2Strand <- function(gene.list, return.original=TRUE) {
  library("biomaRt")
  # https://support.bioconductor.org/p/74322/  need to use host and biomart
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
  gos <- getBM(gene.list, attributes=c("external_gene_name", "strand"),
               filters=c("external_gene_name"),
               mart=mart.obj)
  gl <- gos[match(gene.list, gos[,1]), 2]
  ## if not found, then keep ENSEMBL name
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original){
    gl[is.na(gl)] <- gene.list[is.na(gl)] 
  }
  return(gl)
}



EnsemblGene2Gene <- function(gene.list, return.original=TRUE) {
  library("biomaRt")
  # https://support.bioconductor.org/p/74322/  need to use host and biomart
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
  gos <- getBM(gene.list,attributes=c("ensembl_gene_id", "external_gene_name"),
               filters=c("ensembl_gene_id"),
               mart=mart.obj)
  gl <- gos[match(gene.list, gos[,1]), 2]
  ## if not found, then keep ENSEMBL name
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original){
    gl[is.na(gl)] <- gene.list[is.na(gl)] 
  }
  return(gl)
}

Gene2Ensembl <- function(gene.list, return.original=TRUE) {
  library("biomaRt")
  # https://support.bioconductor.org/p/74322/  need to use host and biomart
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
  gos <- getBM(gene.list,attributes=c("external_gene_name", "ensembl_gene_id"),
               filters=c("external_gene_name"),
               mart=mart.obj)
  gl <- gos[match(gene.list, gos[,1]), 2]
  ## if not found, then keep ENSEMBL name
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original){
    gl[is.na(gl)] <- gene.list[is.na(gl)] 
  }
  return(gl)
}

Attribute2Gene <- function(gene.list, jattribute, return.original=TRUE) {
  library("biomaRt")
  # https://support.bioconductor.org/p/74322/  need to use host and biomart
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
  gos <- getBM(gene.list,attributes=c(jattribute, "external_gene_name"),
               filters=c("ensembl_gene_id"),
               mart=mart.obj)
  gl <- gos[match(gene.list, gos[,1]), 2]
  ## if not found, then keep ENSEMBL name
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original){
    gl[is.na(gl)] <- gene.list[is.na(gl)] 
  }
  return(gl)
}

AppendGeneID <- function(dat){
  genes.tens <- rownames(dat)
  Gene.ID <- Transcript2Gene(genes.tens)
  dat <- cbind(Gene.ID, dat)
  return(dat)
}

AnnotatePseudogenes <- function(gene.list, return.original=TRUE){
  library("biomaRt")
  # https://support.bioconductor.org/p/74322/  need to use host and biomart
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
  gos <- getBM(gene.list,attributes=c("external_gene_name", "gene_biotype"),
               filters=c("external_gene_name"),
               mart=mart.obj)
  gl <- gos[match(gene.list, gos[,1]), 2]
  ## if not found, then keep ENSEMBL name
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original){
    gl[is.na(gl)] <- gene.list[is.na(gl)] 
  }
  return(gl)
}

AnnotateTranscriptLength <- function(gene.list, return.original=TRUE, min.length=1000){
  library("biomaRt")
  # https://support.bioconductor.org/p/74322/  need to use host and biomart
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl', host="www.ensembl.org")
  gos <- getBM(gene.list,attributes=c("external_gene_name", "transcript_length"),
               filters=c("external_gene_name"),
               mart=mart.obj)
  print(gos[match(gene.list, gos[, 1]), ])
  gl <- gos[match(gene.list, gos[,1]), 2]
  ## if not found, then keep ENSEMBL name
  print(paste0("Could not match ", length(gl[is.na(gl)]), " genes."))
  if (return.original){
    gl[is.na(gl)] <- gene.list[is.na(gl)] 
  }
  return(gl)
}
