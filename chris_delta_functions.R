# chris_delta_functions.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: General - work on delta profiles
# Description: Functions to work with delta profiles, pairwise
#              analyses of CGH profiles
# -------------------------------------------------------------------

deltaQuant <- function(KC, sampleCombs) {
  
  # This function returns delta profiles
  # First it does quantile normalization on the KC data frame 
  # Then it calculates all differences in the variable sampleCombs
  # for now sampleCombs is sample numbers, so c(1,2) means columns 3 and 4, it disregards the
  # chrom and maploc columns
  
  # Args: KC - KC data frame with 1st col = chrom, 2nd col = maploc, rest of the cols are log2 sampledata
  #       sampleCombs - matrix/data.frame of column numbers or names. First entry 
  #                     is subtracted from the second.
  
  # To Do:
  # - input checking
  
  require(preprocessCore)
  
  # Quantile normalization
  
  dataMatrix <- as.matrix(KC[,3:ncol(KC)])
  dataMatrix <- normalize.quantiles(dataMatrix)
  KCnorm <- KC
  KCnorm[,3:ncol(KCnorm)] <- dataMatrix
  
  KCdiff <- KCnorm[, c(1,2)]

  for (i in 1:nrow(sampleCombs)) {
    KCdiff$temp <- KCnorm[,sampleCombs[i, 1]] - 
      KCnorm[,sampleCombs[i, 2]]
    newColName <- paste(sampleCombs[i,], collapse='-')
    colnames(KCdiff) <- gsub('temp', newColName, colnames(KCdiff))
  }

  return(KCdiff)

}

deltaLinear <- function(comb, KC, KCseg, thres=.2) {
  
  # This function calculates delta profiles based on robust 
  # linear regression. Inputs are:
  # comb - vector of two column names to compare
  # KC - KC data frame of the data
  # KCseg - segmented version of the KC data frame
  # thres - threshold on which to select the segments to use
  #         for linear regression (over/under values)
  #
  # important: column 4 will be subtracted from column 3,
  # so the metastasis sample will have to be in column 4

  require(robustbase)
  require(DNAcopy)
    
  smallKC <- KC[,c('chrom', 'maploc', comb)]
  smallSeg <- subset(KCseg, samplelist=comb)
  smallFreq <- glFrequency(smallSeg, thres)
  
  # Select only those probes that exceed the threshold for either gains or losses
  ind <- smallFreq$gain == 1 | smallFreq$loss == -1
  
  fitrob <- lmrob(smallKC[ind,4] ~ smallKC[ind,3])
  
  diffNorm <- (smallKC[,3] + coef(fitrob)[[1]]) * coef(fitrob)[[2]] - smallKC[,4]
  
  return(diffNorm)
  
}

deltaLinearSeg <- function(comb, KC, KCseg, thres=.2) {
  
  # This function calculates delta profiles based on robust 
  # linear regression. 
  # This variant first segments the data and sets probevalues
  # to their segment means. These probes are then used to
  # correct for tumor content and to cacluate the difference
  # The delta profile is immediately a segmented profile.
  #
  # Inputs are:
  # comb - vector of two column names to compare
  # KC - KC data frame of the data
  # KCseg - segmented version of the KC data frame
  # thres - threshold on which to select the segments to use
  #         for linear regression (over/under values)
  #
  # important: column 4 will be subtracted from column 3,
  # so the metastasis sample will have to be in column 4

  require(robustbase)
  require(DNAcopy)
  
  smallKC <- KC[,c('chrom', 'maploc', comb)]
  smallSeg <- subset(KCseg, samplelist=comb)

  segOutput <- smallSeg$output
  segKC <- setProbeToSeg(smallKC, smallSeg)

  # Select only those probes that exceed the threshold for either gains or losses
  ind <- abs(segKC[,3]) > .2 & abs(segKC[,4]) > .2
  lmfit <- lm(segKC[ind,4] ~ segKC[ind,3])
  
  diffNorm <- (segKC[,3] + coef(lmfit)[[1]]) * coef(lmfit)[[2]] - segKC[,4]
  
  return(diffNorm)
  
}

setProbeToSeg <- function (KC, KCseg) {
  
  # Function to set probe values to their calculated segments
  # ID in KCseg$output must be equal to the colnames in KC
    
  segOutput <- KCseg$output
  segProbeKC <- KC

  for (i in 1:nrow(segOutput)) {
    probesInSeg <- with(segOutput, segProbeKC$chrom == chrom[i] & 
      segProbeKC$maploc >= loc.start[i] &
      segProbeKC$maploc < loc.end[i])

    segProbeKC[probesInSeg, segOutput$ID[i]] <- 
      segOutput$seg.mean[i]
  }

  return(segProbeKC)

}