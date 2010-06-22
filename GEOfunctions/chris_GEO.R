chrisGEODataCorrMetYap <- function(DataSet, normalString=NULL, nperm=3000) {
	
	outputList <- list(realCor=NULL, randCor=NULL)
	
	# Make sure to select only the cancer samples and not normal samples.
	# normalString can be a regular expression
	
	if (is.null(normalString)) {
		datacols <- grep('GSM', colnames(DataSet))
	}
	else {
		datacols <- setdiff(grep('GSM', colnames(DataSet)), grep(normalString, DataSet[1,]))
	}
	
	YAPexprs <- sapply(DataSet[grep('YAP', DataSet$Gene.symbol), datacols], as.numeric)
	METexprs <- sapply(DataSet[grep('MET', DataSet$Gene.symbol), datacols], as.numeric)
	
	# Get means if more than one probe is found - Any NAs are inherited for now
	if (!is.null(nrow(YAPexprs))) {
		YAPexprs <- colMeans(YAPexprs, na.rm=F)
	}
	if (!is.null(nrow(METexprs))) {
		METexprs <- colMeans(METexprs, na.rm=F)
	}
	
	# Remove NA containing samples for both genes	
	naIndex <- which(!is.na(YAPexprs) & !is.na(METexprs))
	METexprs <- METexprs[naIndex]
	YAPexprs <- YAPexprs[naIndex]
	
	# Calculate the correlation between the top of the gene expression
	YAPtop <- which(YAPexprs > mean(YAPexprs, na.rm=T) + sd(YAPexprs, na.rm=T))
	METtop <- which(METexprs > mean(METexprs, na.rm=T) + sd(METexprs, na.rm=T))
	alltop <- unique(c(METtop,YAPtop))
	plot(YAPexprs[alltop], METexprs[alltop])
	outputList$realCor <- cor(YAPexprs[alltop], METexprs[alltop])
	
	# Do a permutation to determine the significance of the corralation value
	
	randCor <- vector(mode='numeric', length=nperm)
	for (i in 1:nperm) {
		YAPexprsRand <- YAPexprs[sample(1:length(YAPexprs))]
		METexprsRand <- METexprs[sample(1:length(METexprs))]
		YAPtopRand <- which(YAPexprsRand > mean(YAPexprsRand, na.rm=T) + sd(YAPexprsRand, na.rm=T))
		METtopRand <- which(METexprsRand > mean(METexprsRand, na.rm=T) + sd(METexprsRand, na.rm=T))
		randTop <- unique(c(METtopRand, YAPtopRand))
		randCor[i] <- cor(YAPexprsRand[randTop], METexprsRand[randTop])
	}
	
	outputList$randCor <- randCor	
	
	return(outputList)
	
}

chrisGEODataTtestMetYap <- function(DataSet, normalString=NULL, nperm=3000) {
	
	outputList <- list(realCor=NULL, randCor=NULL)
	
	# Make sure to select only the cancer samples and not normal samples.
	# normalString can be a regular expression
	
	if (is.null(normalString)) {
		datacols <- grep('GSM', colnames(DataSet))
	}
	else {
		datacols <- setdiff(grep('GSM', colnames(DataSet)), grep(normalString, DataSet[1,]))
	}
	
	YAPexprs <- sapply(DataSet[grep('YAP', DataSet$Gene.symbol), datacols], as.numeric)
	METexprs <- sapply(DataSet[grep('MET', DataSet$Gene.symbol), datacols], as.numeric)
	
	# Get means if more than one probe is found - Any NAs are inherited for now
	if (!is.null(nrow(YAPexprs))) {
		YAPexprs <- colMeans(YAPexprs, na.rm=F)
	}
	if (!is.null(nrow(METexprs))) {
		METexprs <- colMeans(METexprs, na.rm=F)
	}
	
	# Remove NA containing samples for both genes	
	naIndex <- which(!is.na(YAPexprs) & !is.na(METexprs))
	METexprs <- METexprs[naIndex]
	YAPexprs <- YAPexprs[naIndex]
	
	# Calculate the correlation between the top of the gene expression
	YAPtop <- which(YAPexprs > mean(YAPexprs, na.rm=T) + sd(YAPexprs, na.rm=T))
	METtop <- which(METexprs > mean(METexprs, na.rm=T) + sd(METexprs, na.rm=T))
	
	
	
	alltop <- unique(c(METtop,YAPtop))
	plot(YAPexprs[alltop], METexprs[alltop])
	outputList$realCor <- cor(YAPexprs[alltop], METexprs[alltop])
	
	# Do a permutation to determine the significance of the corralation value
	
#	randCor <- vector(mode='numeric', length=nperm)
#	for (i in 1:nperm) {
#		YAPexprsRand <- YAPexprs[sample(1:length(YAPexprs))]
#		METexprsRand <- METexprs[sample(1:length(METexprs))]
#		YAPtopRand <- which(YAPexprsRand > mean(YAPexprsRand, na.rm=T) + sd(YAPexprsRand, na.rm=T))
#		METtopRand <- which(METexprsRand > mean(METexprsRand, na.rm=T) + sd(METexprsRand, na.rm=T))
#		randTop <- unique(c(METtopRand, YAPtopRand))
#		randCor[i] <- cor(YAPexprsRand[randTop], METexprsRand[randTop])
#	}
	
#	outputList$randCor <- randCor	
	
	return(outputList)
	
}
