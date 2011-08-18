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

getSegCount <- function (segFrame, originalFrame) {
  
  # Count segments per tumor, separated by groups in the group variable in the dataframes
  
  segCount <- as.data.frame(table(segFrame$ID, segFrame$group), stringsAsFactors=F)
  segCount <- segCount[-which(segCount$Freq == 0),]
  colnames(segCount) <- c('TumorID','group', 'NumSeg')
  
  # Add missing tumor IDs
  
  AllID <- unique(originalFrame$ID)
  
  missingID <- which(!(AllID %in% segFrame$ID))
  
  for (m in missingID) {
      rowInsert <- data.frame(TumorID=AllID[m], 
        group=originalFrame$group[originalFrame$ID==AllID[m]][1], 
        NumSeg=0)
      segCount <- rbind(segCount, rowInsert)    
  }
  
  
  return(segCount)

}
