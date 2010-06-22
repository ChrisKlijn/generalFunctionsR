# General plotting functions.

chrisFancyScatterPlot <- function (vectorFrame, plotTitle=NULL, titleSz=20, titleCol='gray', pointCol='darkgreen', pointCex=2, scaleX=T, lineIntersect=0) {
	# This function plots scatterplots of the columns in an input dataframe
	# Only numeric data in the columns. All combinations of the columns
	# in the data frame will be plotted
	
	if (ncol(vectorFrame) > 5) {
		print('Too many columns, don\'t even think about it.\n')
		return(NULL)
	}
	
	if (scaleX) {
		setXLim = c(min(vectorFrame), max(vectorFrame))
		}
	
	
	# Load necessary library
	library(ggplot2)
	
	# Calculate the number of combinations
	numPlots <- ncol(vectorFrame)*(ncol(vectorFrame)-1)/2
	
	# Ugly way to limit the number of columns and the place of the title
	numCol <- numPlots
	titleColumn <- 1
	if (numCol >= 3) {
		numCol <- 3
		titleColumn <- 2
	}
	
	# Set up the grid-supported multiple plots
	# Row 1 will be the title
	# The rest of the rows will be the data. 3 plot every row.
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow=ceiling(numPlots/3)+1, ncol=numCol, 
		heights=unit(c(.25, rep(1, times=ceiling(numPlots/3))), 'null'))))
	vplayout <- function(x, y) {
		viewport(layout.pos.row = x, layout.pos.col = y)
	}
	
	# Set up the ggplot2 plots
	peakCombs <- combn(colnames(vectorFrame), 2)
	allPlotsPeaks <- vector(mode='list', length=ncol(peakCombs))

	# Collect all possible combinations as plots
	for (i in 1:length(allPlotsPeaks)) {
		allPlotsPeaks[[i]] <- ggplot(aes_string(x = peakCombs[1,i], y = peakCombs[2,i]), data=vectorFrame)
	}
	
	# plot the title
	grid.text(plotTitle, gp=gpar(fontsize=titleSz, col=titleCol), vp=vplayout(1, titleColumn))
	
	# loop for plotting
	finished <- F
	row <- 1
	column <- 1
	count <- 1

	while (T) {
		if (scaleX) {
			setXLim = c(min(vectorFrame), max(vectorFrame))
			print(allPlotsPeaks[[count]]+geom_vline(xintercept = lineIntersect)+geom_hline(yintercept = lineIntersect)+geom_point(col=pointCol, cex=pointCex)+xlim(setXLim)+ylim(setXLim), vp = vplayout(row+1, column))
		}
		else {
			print(allPlotsPeaks[[count]]+geom_vline(xintercept = lineIntersect)+geom_hline(yintercept = lineIntersect)+geom_point(col=pointCol, cex=pointCex), vp = vplayout(row+1, column))
		}
		count <- count + 1
		if (count > length(allPlotsPeaks)) break
		column <- column + 1
		if (column > 3) {
			row <- row + 1
			column <- 1
		}	
	}
}

chrisFancyJitterPlot <- function (vectorFrame, xCol, yCol, plotTitle=NULL, titleSz=20, titleCol='gray', pointCex=3, colOffset=50) {
	# This function plots boxplots of the columns in an input dataframe
	# Only numeric data in the columns. All combinations of the columns
	# in the data frame will be plotted
	
	# Load necessary library
	library(ggplot2)
	
	# How many plots?
	numPlots <- length(yCol)
	# Set up the grid-supported multiple plots
	# Row 1 will be the title
	# The rest of the rows will be the data. 3 plot every row.
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(ceiling(numPlots/3)+1, 3, heights=unit(c(.25, rep(1, times=ceiling(numPlots/3))), 'null'))))
	vplayout <- function(x, y) {
		viewport(layout.pos.row = x, layout.pos.col = y)
	}
	
	# Set up the ggplot2 plots
	allPlotsPeaks <- vector(mode='list', length=numPlots)

	# Collect all possible combinations as plots
	for (i in 1:length(allPlotsPeaks)) {
		allPlotsPeaks[[i]] <- ggplot(aes_string(x=xCol, y=yCol[i], colour=xCol), data=vectorFrame)
	}
	
	# plot the title
	grid.text(plotTitle, gp=gpar(fontsize=titleSz, col=titleCol), vp=vplayout(1, 2))
	
	# loop for plotting
	finished <- F
	row <- 1
	column <- 1
	count <- 1

	while (T) {
		print(allPlotsPeaks[[count]] + geom_jitter(size=pointCex) + 
			scale_colour_discrete(h=c(180, 270)+colOffset, l=50, c=50) + opts(legend.position = "none", title=paste(yCol[count], 'expression')) + 
			ylab('') + xlab(''), vp = vplayout(row+1, column))
		count <- count + 1
		if (count > length(allPlotsPeaks)) break
		column <- column + 1
		if (column > 3) {
			row <- row + 1
			column <- 1
		}	
	}
}

chrisFancyBoxPlot <- function (vectorFrame, xCol, yCol, plotTitle=NULL, titleSz=20, titleCol='gray', pointCex=3, colOffset=50) {
	# This function plots boxplots of the columns in an input dataframe
	# Only numeric data in the columns. All combinations of the columns
	# in the data frame will be plotted
	
	# Load necessary library
	library(ggplot2)
	
	# How many plots?
	numPlots <- length(yCol)
	# Set up the grid-supported multiple plots
	# Row 1 will be the title
	# The rest of the rows will be the data. 3 plot every row.
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(ceiling(numPlots/3)+1, 3, heights=unit(c(.25, rep(1, times=ceiling(numPlots/3))), 'null'))))
	vplayout <- function(x, y) {
		viewport(layout.pos.row = x, layout.pos.col = y)
	}
	
	# Set up the ggplot2 plots
	allPlotsPeaks <- vector(mode='list', length=numPlots)

	# Collect all possible combinations as plots
	for (i in 1:length(allPlotsPeaks)) {
		allPlotsPeaks[[i]] <- ggplot(aes_string(x=xCol, y=yCol[i], colour=xCol), data=vectorFrame)
	}
	
	# plot the title
	grid.text(plotTitle, gp=gpar(fontsize=titleSz, col=titleCol), vp=vplayout(1, 2))
	
	# loop for plotting
	finished <- F
	row <- 1
	column <- 1
	count <- 1

	while (T) {
		print(allPlotsPeaks[[count]] + geom_boxplot() + 
			scale_colour_discrete(h=c(180, 270)+colOffset, l=50, c=50) + opts(legend.position = "none", title=paste(yCol[count], 'expression')) + 
			ylab('') + xlab(''), vp = vplayout(row+1, column))
		count <- count + 1
		if (count > length(allPlotsPeaks)) break
		column <- column + 1
		if (column > 3) {
			row <- row + 1
			column <- 1
		}	
	}
}
