# Functions for use with KC smart analyses
# depends on chris_general_plotting.R

# -------------------------------------------------------------------#
chrisKCplotPeakRegions <- function (data, spm, peakChroms, ampcutoff=NULL, doPlot=T, titleName=NULL, pointCol='darkgreen') {

	# Load necessary library
	library(ggplot2)

	# Get the highest peaks from the Spm for given chromosomes
	peakLoc <- vector(mode='list', length(peakChroms))
	peakNames <- paste('chr', peakChroms, sep='')
	peakNames <- gsub('20', 'X', peakNames)
	names(peakLoc) <- peakNames
	peakBAC <- peakLoc
	for (i in 1:length(peakChroms)) {
		peakLoc[[i]] <- 50000 * which.max(spm@data[[as.character(peakChroms[i])]]$pos)
		peakBAC[[i]] <- data[data$maploc > peakLoc[[i]]-1000000 & data$maploc < peakLoc[[i]]+1000000 & data$chrom == peakChroms[i],3:ncol(data)]
	}

	BACmeans <- as.data.frame(sapply(peakBAC, colMeans))
	
	chrisFancyScatterPlot(BACmeans, plotTitle=titleName, pointCol=pointCol)
	
#	# Set values under a certain cutoff to 0, to specifically look at amplification and not small gains etc
#	if (!is.null(ampcutoff)) {
#		BACmeans[BACmeans < ampcutoff] <- 0
#	}
#	# Make scatterplots, uses ggplot2 library
	
#	if (doPlot) {
#		peakCombs <- combn(names(peakLoc), 2)
#		allPlotsPeaks <- vector(mode='list', length=ncol(peakCombs))

#		# Collect all possible combinations as plots
#		for (i in 1:length(allPlotsPeaks)) {
#			allPlotsPeaks[[i]] <- ggplot(aes_string(x = peakCombs[1,i], y = peakCombs[2,i]), data=BACmeans)
#		}

#		# Set up the grid for multiple plotting
#		numPlots <- (length(peakChroms)*(length(peakChroms)-1))/2
#		grid.newpage()
#		pushViewport(viewport(layout = grid.layout(ceiling(numPlots/3)+1, 3, heights=unit(c(.25, rep(1, times=ceiling(numPlots/3))), 'null'))))
#		vplayout <- function(x, y)
#		  viewport(layout.pos.row = x, layout.pos.col = y)

#		# loop de loop for plotting
#		finished <- F
#		row <- 1
#		column <- 1
#		count <- 1

#		grid.text(titleName, gp=gpar(fontsize=20, col="grey"), vp=vplayout(1, 2))

#		while (T) {
#			print(allPlotsPeaks[[count]]+geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+geom_point(col='darkred', cex=3), vp = vplayout(row+1, column))
#			count <- count + 1
#			if (count > length(allPlotsPeaks)) break
#			column <- column + 1
#			if (column > 3) {
#				row <- row + 1
#				column <- 1
#			}	
#		}
#	}
	return(BACmeans)
}

# -------------------------------------------------------------------#

chrisKCcluster <- function (data, spm, topn=list(gain=10, loss=10), meanWindow=10e6) {

# Get the highest values and the lowest values per chromosome into a dataframe
	topGains <- as.data.frame(sapply(spm@data, function (x) {c(x$pos[which.max(as.numeric(x$pos))],which.max(as.numeric(x$pos)))}))
	topLosses <- as.data.frame(sapply(spm@data, function (x) {c(x$neg[which.min(as.numeric(x$neg))],which.min(as.numeric(x$neg)))}))

	topGains <- topGains[,order(topGains[1,], decreasing=T)[1:topn$gain], ]
	topLosses <- topLosses[,order(topLosses[1,], decreasing=F)[1:topn$loss]]

	peakLocGains <- vector(mode='list', length=ncol(topGains))
	peakNamesGains <- colnames(topGains)
	peakBACgains <- peakLocGains
	peakChromGains <- as.numeric(gsub('X', '20', peakNamesGains))
	for (i in 1:length(peakLocGains)) {
		peakLocGains[[i]] <- 50000 * topGains[2,i]
		peakBACgains[[i]] <- data[data$maploc > peakLocGains[[i]]-meanWindow & data$maploc < peakLocGains[[i]]+meanWindow & data$chrom == peakChromGains[i],3:ncol(data)]
	}

	peakLocLosses <- vector(mode='list', length=ncol(topLosses))
	peakNamesLosses <- colnames(topLosses)
	peakBACLosses <- peakLocLosses
	peakChromLosses <- as.numeric(gsub('X', '20', peakNamesLosses))
	for (i in 1:length(peakLocLosses)) {
		peakLocLosses[[i]] <- 50000 * topLosses[2,i]
		peakBACLosses[[i]] <- data[data$maploc > peakLocLosses[[i]]-meanWindow & data$maploc < peakLocLosses[[i]]+meanWindow & data$chrom == peakChromLosses[i],3:ncol(data)]
	}

	BACmeanGains <- as.data.frame(sapply(peakBACgains, colMeans, na.rm=T))
	colnames(BACmeanGains) <- paste('GainChr', peakNamesGains, sep='')
	BACmeanLosses <- as.data.frame(sapply(peakBACLosses, colMeans, na.rm=T))
	colnames(BACmeanLosses) <- paste('LossChr', peakNamesLosses, sep='')

	BACmeansAll <- cbind(BACmeanGains, BACmeanLosses)
	
	BACmeansAll <- BACmeansAll[,colSums(is.na(BACmeansAll)) < ceiling(ncol(BACmeansAll)/10)]

	heatmap.2(as.matrix(BACmeansAll), trace = 'none', col = greenred(75), labRow=rownames(BACmeansAll), dendrogram='row', breaks=seq(-2,2,by=4/75), main='TITLE HERE', margins=c(12, 10))

}

# -------------------------------------------------------------------#

chrisKCclusterSpm <- function (spm, spmColl, sampleNames=NA, topn=list(pos=10, neg=10), plotTitle=NULL, colSide=NULL, clustMeth='euclidean') {
	
	library(gplots)
	
	allData <- spm@data
	peaksPos <- NULL
	peaksNeg <- NULL
	for (c in 1:length(allData)) {
		if (!(topn$pos == 0)) {
			peaksTemp <- KCfindPeaks(allData[[names(allData)[c]]]$pos, mode='pos')
			names(peaksTemp) <- paste(names(peaksTemp), 'Chr', names(allData)[c], sep='')
			peaksPos <- c(peaksTemp, peaksPos)
			peaksPos <- peaksPos[order(peaksPos, decreasing=T)[1:topn$pos]]
		}
		if (!(topn$neg == 0)) {
			peaksTemp <- KCfindPeaks(allData[[names(allData)[c]]]$neg, mode='neg')
			names(peaksTemp) <- paste(names(peaksTemp), 'Chr', names(allData)[c], sep='')
			peaksNeg <- c(peaksTemp, peaksNeg)
			peaksNeg <- peaksNeg[order(peaksNeg, decreasing=F)[1:topn$neg]]
		}
	}
	allPeaks <- c(peaksPos, peaksNeg)
	peaksChrom <- gsub('[0-9]*Chr', '', names(allPeaks))
	peaksMaploc <- as.numeric(gsub('Chr[0-9X]*', '', names(allPeaks)))
	ann <- spmColl@annotation
	indexPeaks <- vector(mode='integer', length=length(peaksChrom))

	for (i in 1:length(peaksChrom)) {
		indexPeaks[i] <- which(ann@chromosome == peaksChrom[i] & ann@maploc == peaksMaploc[i])
	}

	# For now: get rid of NA containting peaks
	indexPeaks <- indexPeaks[which(rowSums(is.na(spmColl@data[indexPeaks,])) == 0)]
	
	# Do the clustering using dist and hclust to get more control
	
	if (clustMeth == 'euclidean') {
		dataDist <- dist(t(spmColl@data[indexPeaks,]), method = "euclidean")
	}
	
	dataClust <- hclust(dataDist, method='complete')
	
	if (is.null(colSide)) { 
		output <- heatmap.2(spmColl@data[indexPeaks,], trace = 'none', col = greenred(75), labRow=names(allPeaks), 
				labCol=sampleNames, breaks=seq(-1,1,by=2/75), main=plotTitle, margins=c(12, 10),
				dendrogram='col',  Colv=as.dendrogram(dataClust))
	}
	else {
		output <- heatmap.2(spmColl@data[indexPeaks,], trace = 'none', col = greenred(75), labRow=names(allPeaks), 
				labCol=sampleNames, breaks=seq(-1,1,by=2/75), main=plotTitle, margins=c(12, 10),
				ColSideColors=colSide, dendrogram='col', Colv=as.dendrogram(dataClust))
	}
	return(dataClust)
}

# -------------------------------------------------------------------#
KCfindPeaks <- function(data, mode='pos', sigma=50000){
	# A modified version of the find peaks function in the KC smart 
	# package. Input here is either the $pos or $neg column in one of the
	# smp@data chromosome lists
	dir <- ifelse(mode=="pos", -2, 2)
	loc <- (seq(1, length(data)) * sigma) + 1
	names(data) <- loc
	data <- data[!is.na(data)]
	dataOriginal <- data[-1]

	data <- diff.default(data)
	data[data > 0] <- 1
	data[data < 0] <- -1
	nonFlats <- data != 0

	data2 <- diff.default(data[nonFlats])

	dataOriginal[nonFlats][data2 == dir]
}


