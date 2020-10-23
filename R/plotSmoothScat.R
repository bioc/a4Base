# Genome wide plots with density coloring
# 
# Author: Willem Talloen, Suzy Van Sanden
###############################################################################



#' Plots the correlation in gene expression between two samples
#' 
#' Plots the correlation in gene expression between two samples. Each dot represents
#' a gene, and the dots have a density-dependent coloring.
#' Genes with exceptional behavior can be highlighted by showing their gene symbol. 
#' @param object ExpressionSet object for the experiment
#' @param x String containing the name of the first sample. This should be a 
#'  the name of a column in the \code{exprs} data of the \code{expressionSet} object.
#' @param y String containing the name of the second sample. See \code{x}
#' @param trsholdX Vector of two values specifying the X-axis thresholds within which
#' genes should be highlighted by their gene symbol.
#' @param trsholdY Vector of two values specifying the Y-axis thresholds within which
#' genes should be highlighted by their gene symbol.
#' @param probe2gene Boolean indicating whether the probeset should be translated to a gene symbol
#' (used for the default title of the plot)
#' @param ... Possibility to add extra plot options. See \code{\link{par}}
#' @return No returned value, a plot is drawn to the current device.
#' @examples 
#' if (require(ALL)){
#'   data(ALL, package = "ALL")
#'   ALL <- addGeneInfo(ALL)
#'  plotComb2Samples(ALL,"84004", "01003",
#'     trsholdX = c(10,12), trsholdY = c(4,6),
#' 	   xlab = "a B-cell", ylab = "a T-cell")
#' }
#' @seealso \code{\link{plotCombMultSamples}}
#' @author W. Talloen
#' @importFrom Biobase exprs featureData
#' @importFrom graphics plot axis box points text
#' @importFrom grDevices densCols
#' @export
plotComb2Samples <- function(object, x, y,
		trsholdX = NULL, trsholdY = NULL,
		probe2gene = TRUE, ...){
	
	x <- exprs(object)[, as.character(x)]
	y <- exprs(object)[, as.character(y)]
	
	gSymbol <- featureData(object)$`SYMBOL`

	### determine points to label
	if (is.null(trsholdX) & is.null(trsholdY)){
		pointsToLabel <- NULL
	} else {
		XpointsToLabel <- if (!is.null(trsholdX)){
		  min(trsholdX) < x & x < max(trsholdX)
		} else {
			FALSE
		}
        YpointsToLabel <- if (!is.null(trsholdY)){
		  min(trsholdY) < y & y < max(trsholdY)
		} else {
			FALSE
		}
		pointsToLabel <- XpointsToLabel & YpointsToLabel
	}
	
	plot(x, y, axes = FALSE, type="n", ...)
	axis(1, lwd = 1.5, las = 1); axis(2, lwd = 1.5, las = 1)
	box(bty='l',lwd = 1.5)
	
	if (!is.null(trsholdX) | !is.null(trsholdY)) {
		dotColors <- densCols(x[-pointsToLabel], y[-pointsToLabel])
		points(x[-pointsToLabel], y[-pointsToLabel],pch = 20, cex = 1, col = dotColors)
		text(x[pointsToLabel], y[pointsToLabel], labels = gSymbol[pointsToLabel],
				cex = 0.65, col = "black")
	} else {
		dotColors <- densCols(x, y)
		points(x, y,pch = 20, cex = 1, col = dotColors)
	}
}

#### Scatterplot matrix with density-dependent coloring
#' @importFrom graphics points abline
#' @importFrom grDevices densCols
panel.plotSmoothScat <- function(x, y, ...) {
	points(x, y, type="n", main="", xlab="", ylab="", ...)
	dotColors <- densCols(x, y)
	points(x, y, pch = 20, cex = 1, col = dotColors)
	abline(a=0, b=1, col="red")
}

#' @importFrom graphics par text strwidth
#' @importFrom stats cor
#' @importFrom stats cor.test symnum
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- abs(cor(x, y))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
	
	test <- cor.test(x,y)
	# borrowed from printCoefmat
	Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
			cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
			symbols = c("***", "**", "*", ".", " "))
	
	text(0.5, 0.5, txt, cex = cex * r)
	text(.8, .8, Signif, cex=cex, col=2)
}

#' Plots the correlation in gene expression between more than 2 samples
#' @param exprsMatrix ExpressionSet object to plot. For larger datasets,
#' this will typically be a subset of the data.
#' @param ... Further arguments, e.g. to add extra plot options. See \code{\link{pairs}}
#' @return no returned value, a plots is drawn in the current device
#' @author Willem Talloen
#' @seealso \code{\link{plotComb2Samples}}
#' @examples 
#' if (require(ALL)){
#'  data(ALL, package = "ALL")
#'  ALL <- addGeneInfo(ALL)
#'  plotCombMultSamples(exprs(ALL)[,c("84004", "11002", "01003")])
#' }
#' @importFrom graphics pairs
#' @export
plotCombMultSamples <- function(exprsMatrix, ...){
	pairs(
		exprsMatrix, 
		lower.panel = panel.cor, 
		upper.panel = panel.plotSmoothScat, 
		...
	)
}
