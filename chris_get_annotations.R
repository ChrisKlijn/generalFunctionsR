# chris_get_annotations.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: General R functions
# Description: Scripts to get and update annotation files
# -------------------------------------------------------------------

### MOUSE GENES ###

# Function

getMouseGenesUCSC <- function () {
  
  require(rtracklayer)
  
  session <- browserSession("UCSC")
  genome(session) <- "mm9"
  query <- ucscTableQuery(session, "refGene")
  mouseGenes <- getTable(query)

  for (i in 1:ncol(mouseGenes)) {
    if (is.factor(mouseGenes[,i])) {
      mouseGenes[,i] <- as.character(mouseGenes[,i])
    }
  }
  
  return(mouseGenes)
  
}

# Run and save on medoid

mouseGenes <- getMouseGenesUCSC()
save(file='~/data/annData/mouseGenesUCSC.Rda', list='mouseGenes')

