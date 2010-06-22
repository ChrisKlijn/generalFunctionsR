# File to parse GEO matrix files.
# Written because GEOquery does not work
#
# Make sure you preprocessed the series file unsing the bash_GEO.sh script

chris_import_GEO <- function(filenameData, filenamePlatform) {

tempData <- read.delim(filenameData, header=T, stringsAsFactors=F)
tempPlatform <- read.delim(filenamePlatform, header=T, stringsAsFactors=F)
tempPlatform$ID <- gsub('-1', '', tempPlatform$ID)
outData <- merge(x=tempData, y=tempPlatform, by.x="ID_REF", by.y="ID")
return(outData)

}
