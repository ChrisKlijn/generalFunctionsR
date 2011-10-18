# chris_cghsegmentation_functions.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: CGH segmentation
# Description: Functions to segment your CGH data
# -------------------------------------------------------------------

doCghSeg <- function(chromNum, allKC, chrom, maploc) {
  
  require(cghseg) 
  require(GenomicRanges)
  
  tumnames <- colnames(allKC)
  profiles <- as.matrix(allKC[chrom == chromNum,-c(1,2)])

  # Segmentation
  n <- nrow(profiles)
  CGHd <- new("CGHdata",Y=profiles)
  CGHo  <- new("CGHoptions")
  alpha(CGHo) <- 500 / n    ## look for 1 to a maximum of 500 segments per tumor, 
  select(CGHo) <- "mBIC"   ## this define the way the number of segment is selected ( in between 1 and Kmax)
  wavenorm(CGHo)<-'position' ## Wave correction using position

  CGHr  <- multiseg(CGHd,CGHo)   ##The segments with their mean can be found in CGHr@mu.

  allstart <- unlist(lapply(CGHr@mu, function (x) {return(x$begin)}))
  allend <- unlist(lapply(CGHr@mu, function (x) {return(x$end)}))
  allmeans <- unlist(lapply(CGHr@mu, function (x) {return(x$mean)}))

  allID <-rep(gsub('X','', names(CGHr@mu)), 
    unlist(lapply(CGHr@mu, nrow)))

  segments <- GRanges(
    seqnames = Rle(paste('chr', chromNum, sep=''), length(allstart)),
    ranges=IRanges(start=maploc[allstart], end=maploc[allend]), 
    seg.mean=allmeans, 
    ID=allID)

  # Make a range object to contain the start and the end probe for each segment

  return(segments)
}

