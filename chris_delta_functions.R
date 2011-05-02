# chris_delta_functions.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: General - work on delta profiles
# Description: Functions to work with delta profiles, pairwise
#              analyses of CGH profiles
# -------------------------------------------------------------------

deltaMakeDiff <- function(KC, sampleCombs) {
  
  # This function returns delta profiles
  # First it does quantile normalization on the KC data frame 
  # Then it calculates all differences in the variable sampleCombs
  # for now sampleCombs is sample numbers, so c(1,2) means columns 3 and 4, it disregards the
  # chrom and maploc columns
  
  # Args: KC - KC data frame with 1st col = chrom, 2nd col = maploc, rest of the cols are log2 sampledata
  #       sampleCombs - integer matrix with two colums, each row is a comparison to be done 
  
  # To Do:
  # - input checking
  
  require(preprocessCore)
  
  # Quantile normalization
  
  dataMatrix <- as.matrix(KC[,3:ncol(KC)])
  dataMatrix <- normalize.quantiles(dataMatrix)
  KCnorm <- KC
  KCnorm[,3:ncol(KCnorm)] <- dataMatrix
  
  for (i in 1:nrow(sampleCombs)) {
  
  
  
  }
  
  


}
