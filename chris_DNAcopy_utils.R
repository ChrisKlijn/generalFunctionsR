# chris_DNAcopy_utils.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Description: Functions to facilitate work with DNAcopy data
# -------------------------------------------------------------------


extractSeg <- function (segResult, minMark, cutoff=NULL, higher=c('higher', 'lower', 'both')) { 
  
  # Extract the output segments, for which there are at least minMark markers
  # Use cutoff as a cutoff and higher as a bool to indicate higher or lower (higher=F) than the
  # threshold
  
  if (!is.null(cutoff)) {
  
    outputFrame <- switch(higher,
      higher=segResult$output[which(segResult$output$seg.mean >= cutoff),],
      lower=segResult$output[which(segResult$output$seg.mean <= cutoff),],
      both=segResult$output[which(abs(segResult$output$seg.mean) >= abs(cutoff)),]
    )
  
  }
  else {
    outputFrame <- segResult$output
  }
  
  outputFrame <- outputFrame[outputFrame$num.mark >= minMark,]
  
  return(outputFrame)

}

