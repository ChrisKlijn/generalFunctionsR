# Useful functions for aCGH Data analysis

# -------------------------------------------------------------------
plotRawCghRegion <- function (start, end, data, plotch='o') {
	iter <- ceiling((ncol(data)-2)/24)
	ymax <- max(data[start:end, 3:ncol(data)])
	for (i in 1:iter) {
		x11(width=10, height=10)
		par(mfrow=c(6, 4))
		for (p in 1:24) {
			plot(data$maploc[start:end], data[start:end, (p+2)+(24*(i-1))], xlab='', ylab='', main=colnames(data)[(p+2)+(24*(i-1))], pch=plotch, ylim=c(0, ymax))
		}
	}
}

# -------------------------------------------------------------------
plotRawCghRegionProbes <- function (probes, data, plotch='o') {
	iter <- ceiling((ncol(data)-2)/24)
	ymax <- max(data[probes, 3:ncol(data)])
	for (i in 1:iter) {
		x11(width=10, height=10)
		par(mfrow=c(6, 4))
		for (p in 1:24) {
			plot(data[probes,2], data[probes, (p+2)+(24*(i-1))], xlab='', ylab='', main=colnames(data)[(p+2)+(24*(i-1))], ylim=c(0, ymax), pch=plotch)
		}
	}
}

# -------------------------------------------------------------------
plotRawCghRegionSamples <- function (start, end, chromosome, data, samples, plotch='o', gene=0, filename='testraw.pdf') {
	probes <- which(data$maploc > start & data$maploc < end & data$chrom == chromosome)
	ymax <- max(data[probes, 3:ncol(data)])
	ymin <- min(data[probes, 3:ncol(data)])
	pdf(width=10, height=(2 * length(samples)), file=filename)
	par(mfrow=c(length(samples), 1))
	for (i in 1:length(samples)) {
		plot(data$maploc[probes], data[probes, samples[i]], xlab='', ylab='', main=samples[i], pch=plotch, ylim=c(ymin, ymax))
		abline(h=0, col='black')
		if (gene[1] != 0) {
			points(gene, rep(ymax/2, times=length(gene)), col='red', pch=21, bg='red')
			for (t in 1:length(gene)) {
				text(gene[t], (ymax/2 + ymax/4), names(gene[t]), col='red')
			}
		}
	}
	dev.off()
}

# -------------------------------------------------------------------
getTopChromosomePeakSpm <- function (spm, chromosome, type='g') {
	spm.data <- spm@data
	if (is.numeric(chromosome)) { chromosome <- paste(chromosome) }
	spm.data.chromosome <- (spm.data[chromosome])[[1]]
	if (type == 'g') { 
		position <- which(spm.data.chromosome$pos == max(spm.data.chromosome$pos[!is.na(spm.data.chromosome$pos)]))
	}
	if (type == 'l') { 
		position <- which(spm.data.chromosome$neg == min(spm.data.chromosome$neg[!is.na(spm.data.chromosome$neg)]))
	}
	position * spm@sampleDensity
}

# -------------------------------------------------------------------
plotRawCghRegionGenes <- function (start, end, chromosome, dataSet, mart, sample, geneflag=1, runMean=F, filterWidth=5) {
	
	# Handle KCsmart data frames or Esets
	if (is.data.frame(dataSet)) {
		data.flag = TRUE
		probes <- which(dataSet$maploc > start & dataSet$maploc < end & dataSet$chrom == chromosome)
		ratios <- dataSet[,3:ncol(dataSet)]
		ymax <- max(ratios[probes,])
		ymin <- min(ratios[probes,])
	}
	else {
		data.flag = FALSE
		probes <- dataSet[featureData(dataSet)$chromosome == chromosome & (featureData(dataSet)$midpos > start& featureData(dataSet)$midpos < end),]
		ymax <- max(exprs(probes))
		ymin <- min(exprs(probes))
	}
# Get the correct probes and the genes in the region
	
	genes <- getBM(attributes=c('start_position', 'end_position', 'strand', 
		'mgi_symbol', 'ensembl_gene_id'), filters=c('chromosome_name', 'start', 'end'), 
		values=list(chromosome, start, end), mart=mart)
	genes$midpos <- (genes$start_position + genes$end_position)/2
	genes.plus <- genes[genes$strand == 1,]
	genes.min <- genes[genes$strand == -1,]

	# plotting
	mylayout <- layout(seq(1,length(sample)+1), heights=c(rep(1, times=length(sample)), 1))
	for (i in 1:length(sample)) {
		par(mar=c(1,4,4,2),bty="n", xaxt="n", yaxt="s")
		if (data.flag) {
			plot(dataSet$maploc[probes], ratios[probes,sample[i]], pch=20, ylab=sample[i],ylim=c(ymin, ymax))
			axis(2)
			if (runMean) {
				lines(dataSet$maploc[probes], filter(ratios[probes,sample[i]], rep(1, filterWidth)/filterWidth), col='red')
			}
		}
		else {
			plot(featureData(probes)$midpos, exprs(probes[,sample[i]]), pch=20, ylab=sample[i], ylim=c(ymin, ymax))
			axis(2)
			if (runMean) {
				lines(featureData(probes)$midpos, filter(exprs(probes[,sample[i]]), rep(1, filterWidth)/filterWidth), col='red')
			}
		}
		abline(h=0, col='black')
	}
	par(mar=c(5,4,2,2), bty="n", xaxt="s", yaxt="n")
	plot(0, 0, type='n', xlim=c(start, end), axes=F, xlab=NA, ylab=NA)
	rect(genes.plus$start_position, rep(0.1, times=nrow(genes.plus)), genes.plus$end_position, rep(0.5, times=nrow(genes.plus)), col='black')
	rect(genes.min$start_position, rep(-0.1, times=nrow(genes.min)), genes.min$end_position, rep(-0.5, times=nrow(genes.min)), col='black')
	if (geneflag==1) { 
	text(genes.plus$midpos, rep(c(0.65, 0.85), times=nrow(genes.plus)/2), labels=genes.plus$mgi_symbol, cex=1) 
	text(genes.min$midpos, rep(c(-0.65, -0.85), times=nrow(genes.plus)/2), labels=genes.min$mgi_symbol, cex=1)	
	}
	else {
	text(genes.plus$midpos, rep(c(0.6, 0.8, 1.0), times=nrow(genes.plus)/3), labels=genes.plus$mgi_symbol, cex=1) 
	text(genes.min$midpos, rep(c(-0.6, -0.8, -1.0), times=nrow(genes.plus)/3), labels=genes.min$mgi_symbol, cex=1)	
	}	
}

# -------------------------------------------------------------------
plotRawCghDotPlot <- function (KCdataSet, mirrorLocs, samples=1, doFilter=F, filterSize=10, 
chromosomes=NA, setcex=1, plotTitle=NA, setylim=NULL) {


# First remove missing values
	if (any(is.na(KCdataSet[, samples+2]))) {
		KCdataSet <- KCdataSet[-which(is.na(KCdataSet[,samples+2])),]
	}
## This function plot the raw ratios on a genomic axis

## Standard color for every chromosome
	chromCols <- rainbow(length(unique(KCdataSet$chrom)))

## If chromosomes are selected, use those. Otherwise use all chromosomes
## in the dataset
	if (any(is.na(chromosomes))) {
		uniqChroms <- unique(KCdataSet$chrom)
	}
	else {
		uniqChroms <- chromosomes
	}
	
	uniqChroms <- sort(uniqChroms)
	mirrorLocs <- mirrorLocs[uniqChroms]
	attr(mirrorLocs, 'chromNames') <- as.character(uniqChroms)
	
## Get the end coordinates of the chromosomes, cumulative
	chromEnds <- cumsum(unlist(lapply(mirrorLocs, max)))
## Construct a conlinear vector
	KCdataSet <- KCdataSet[KCdataSet$chrom %in% uniqChroms,]
	KCdataSet <- KCdataSet[order(KCdataSet$chrom, KCdataSet$maploc),]
	maplocLin <- KCdataSet$maploc
	

	for (i in 1:length(uniqChroms)) {
		maplocLin[KCdataSet$chrom == uniqChroms[i]] <- 
			maplocLin[KCdataSet$chrom == uniqChroms[i]] + c(0,chromEnds)[i]
	}

## Use a running mean filter if selected, with supplied (or standard)
## filter size
	if (doFilter) {
		filteredData <- filter(KCdataSet[,samples+2], rep(1, filterSize)/filterSize)
		dataRange <- range(filteredData[!is.na(filteredData)])
	}
	else {
		dataRange <- range(KCdataSet[,samples+2])
	}
	if (!is.null(setylim)) {
    dataRange <- setylim
  }
  
  if (is.na(plotTitle)) {
    plotTitle=colnames(KCdataSet)[samples+2]
  }
  	plot(c(0, max(maplocLin)), dataRange,type='n', xlab='Genomic Position (bp)', ylab='log2', main=plotTitle)


	
	for (i in 1:length(uniqChroms)) {
		if (doFilter) {
			points(maplocLin[KCdataSet$chrom == uniqChroms[i]], filter(KCdataSet[KCdataSet$chrom == uniqChroms[i],samples+2], 
				rep(1, filterSize)/filterSize), col=chromCols[uniqChroms[i]], pch='.', cex=setcex)
		}
		else {
			points(maplocLin[KCdataSet$chrom == uniqChroms[i]], KCdataSet[KCdataSet$chrom == uniqChroms[i],samples+2], 
				col=chromCols[uniqChroms[i]], pch='.', cex=setcex)
		}
	}
	abline(v=chromEnds, col=colors()[615])
	abline(h=0, col='black')
	textX <- chromEnds - (unlist(lapply(mirrorLocs, max)))/2
	text(textX, max(dataRange)-0.05*max(dataRange), labels=attr(mirrorLocs, 'chromNames'))

}

# -------------------------------------------------------------------
exportCNformat <- function (KCdataSet, fileName = 'KCdata.cn', species='mouse') {
  
  # Export a KC data frame as a CN data file
  
  if (!(species == 'mouse' | species == 'human')) {
    stop('species must be either mouse or human')
  }
  
  # Print the header line
  
  sampsText <- paste(colnames(KCdataSet)[3:ncol(KCdataSet)], collapse='\t')
  header <- paste('SNP', 'Chromosome', 'PhysicalPosition', sampsText, sep='\t')
  
  if (is.numeric(KCdataSet$chrom)) {
    KCdataSet$chrom <- as.character(KCdataSet$chrom)
    if (species == 'mouse') {
      KCdataSet$chrom <- gsub('20', 'X', KCdataSet$chrom)
      KCdataSet$chrom <- gsub('21', 'Y', KCdataSet$chrom)
    }
    if (species == 'human') {
      KCdataSet$chrom <- gsub('23', 'X', KCdataSet$chrom)
      KCdataSet$chrom <- gsub('24', 'Y', KCdataSet$chrom)
    }
  }

  # Write data to file - header first
  cat(file=fileName, header, append=F)
  cat(file=fileName, '\n', append=T)
  write.table(x=KCdataSet, file=fileName, col.names=F, quote=F, append=T, sep='\t')


}
