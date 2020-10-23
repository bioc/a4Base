#' Draw a Volcano Plot
#' 
#' Generic function to draw a volcano plot. A volcano plot is a graph that allows to 
#'   simultaneously assess the P values (statistical significance) and log ratios
#'   (biological difference) of differential expression for the given genes.
#' @param x either an object of class 'tTest', of class 'limma' or a numeric vector of 
#'   log ratios, i.e. the log of the fold change values; the names of the logRatio vector 
#'    will be used to display the names of the most interesting gene
#' @param y should not be given if an object of class 'tTest' or 'limma' is passed as 
#'     argument 'x'; if 'x' is a numeric vector of log ratios, 'y' should be given and 
#'     should be a numeric vector of P-values indicating the statistical significance
#' @param pointLabels Labels for points on the volcano plot that are interesting
#'    taking into account both the x and y dimensions; typically this is a
#'    vector of gene symbols; most methods can access the gene symbols directly from 
#'    the object passed as 'x' argument; the argument allows for custom labels if
#'    needed
#' @param ... further arguments to specific methods
#' @examples
#' if (require(ALL)){
#'   data(ALL, package = "ALL")
#'   ALL <- addGeneInfo(ALL)
#'   ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'   tTestRes <- tTest(object = ALL,	groups = "BTtype", probe2gene = TRUE)
#'   volcanoPlot(tTestRes)  
#' }
#' @references Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
#'   Microarrays, Chapman \& Hall/CRC, pp. 148 - 153.
#' @return The volcano plot is drawn to the current device.
#' @seealso See \code{\link{volcanoplotter}}
#' @author Tobias Verbeke, based on code by Willem Talloen
#' @keywords dplot
#' @exportMethod volcanoPlot
setGeneric("volcanoPlot", function(x, y, pointLabels, ...){
      standardGeneric("volcanoPlot")
})

#' Draw a Volcano Plot
#' 
#' This function draws a volcano plot, a graph that allows to simultaneously
#'  assess the statistical and biological significance of differential expression
#'  for the given genes.
#' @details   
#' The set of genes for which labels are displayed is the \emph{union} of the set of
#'  genes that have lowest P-values (\code{topPValues}) and the set of genes
#'  that display the highest absolute values for the log ratios (\code{topLogRatios}).
#' @param x either an object of class 'tTest', of class 'limma' or a numeric vector of 
#'   log ratios, i.e. the log of the fold change values; the names of the logRatio vector 
#'    will be used to display the names of the most interesting gene
#' @param y should not be given if an object of class 'tTest' or 'limma' is passed as 
#'     argument 'x'; if 'x' is a numeric vector of log ratios, 'y' should be given and 
#'     should be a numeric vector of P-values indicating the statistical significance
#' @param pointLabels Labels for points on the volcano plot that are interesting
#'    taking into account both the x and y dimensions; typically this is a
#'    vector of gene symbols; most methods can access the gene symbols directly from 
#'    the object passed as 'x' argument; the argument allows for custom labels if
#'    needed
#' @param topPValues top n points that will be included in the points to label based
#' on their low P Values
#' @param topLogRatios top n points that will be included in the points to label based
#'  on their high absolute values of the log ratio
#' @param smoothScatter use color saturation to indicate dots that are in densely
#' populated regions of the graph; defaults to \code{TRUE} 
#' @param xlab label for the x axis (string)
#' @param ylab label for the y axis (string)
#' @param main main title for the graph (string)
#' @param sub subtitle for the graph (string)
#' @param newpage should the graph be drawn to a new grid page? Defaults to
#' \code{TRUE}. This argument is useful for including several volcano plots 
#'  in one layout.
#' @param additionalPointsToLabel Entrez IDs of genes of interest, that will be highlighted on the plot; the color of highlighting is determined
#' by the 'additionalLabelColor' argument.
#' @param additionalLabelColor Color used to highlight the 'additionalPointsToLabel'; defaults to "red"
#' @return The volcano plot is drawn to the current device.
#' @author Tobias Verbeke, based on code by Willem Talloen
#' @keywords dplot
#' @docType methods
#' @name volcanoPlot-methods
#' @aliases volcanoPlot,tTest,missing
#' @aliases volcanoPlot,tTest,character
#' @aliases volcanoPlot,limma,missing,missing
#' @aliases volcanoPlot,limma,missing,character
#' @aliases volcanoPlot,numeric,numeric,character
#' @aliases volcanoPlot,numeric,numeric,missing
NULL

#' volcanoPlot for an object resulting from \code{tTest}
#' @rdname volcanoPlot-methods
#' @export
setMethod("volcanoPlot",
	signature(x = "tTest", 
		y = "missing",
		pointLabels = "missing"),
	function(x, y, pointLabels, topPValues = 10, 
		topLogRatios = 10,
		smoothScatter = TRUE, xlab = NULL, ylab = NULL,
		main = NULL, sub = NULL, newpage = TRUE, 
		additionalPointsToLabel = NULL, additionalLabelColor = "red"){
	
		logRatio <- x[,"logRatio"]
		pValue <- x[,"p"]
		pointLabels <- x[,"gSymbol"]
			
		volcanoplotter(logRatio = logRatio, pValue = pValue, 
			pointLabels = pointLabels, topPValues = topPValues,
			topLogRatios = topLogRatios, smoothScatter = smoothScatter, 
			xlab = xlab, ylab = ylab, main = main, sub = sub,
			newpage = newpage, additionalPointsToLabel = additionalPointsToLabel, 
			additionalLabelColor = additionalLabelColor)
			
})

#' volcanoPlot for an object resulting from \code{tTest}
#' @rdname volcanoPlot-methods
#' @export
setMethod("volcanoPlot",
    signature(x = "tTest", 
        y = "missing",
        pointLabels = "character"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10,
        smoothScatter = TRUE, xlab = NULL, ylab = NULL,
        main = NULL, sub = NULL, newpage = TRUE, 
        additionalPointsToLabel = NULL, additionalLabelColor = "red"){
      
      logRatio <- x[,"logRatio"]
      pValue <- x[,"p"]
      
      if (length(pValue) != length(pointLabels))
        stop("'pointLabels' should have the same length as the number of rows of 'x'")
      
      volcanoplotter(logRatio = logRatio, pValue = pValue, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, smoothScatter = smoothScatter, 
          xlab = xlab, ylab = ylab, main = main, sub = sub,
          newpage = newpage, additionalPointsToLabel = additionalPointsToLabel, 
          additionalLabelColor = additionalLabelColor)
      
})


#' volcanoPlot for an object resulting from \code{limma2Groups}
#' @rdname volcanoPlot-methods
#' @export
setMethod("volcanoPlot",
    signature(x = "limma", 
        y = "missing",
        pointLabels = "missing"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10,
        smoothScatter = TRUE, xlab = NULL, ylab = NULL,
        main = NULL, sub = NULL, newpage = TRUE, 
        additionalPointsToLabel = NULL, additionalLabelColor = "red"){

      logRatio <- as.matrix(x@MArrayLM$coef)[, 2]
      pValue <- as.matrix(x@MArrayLM$p.value)[, 2]
      pointLabels <- x@geneSymbols
      
      volcanoplotter(logRatio = logRatio, pValue = pValue, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, logTransformP = TRUE, # p values
          smoothScatter = smoothScatter, xlab = xlab, ylab = ylab, 
          main = main, sub = sub,
          newpage = newpage, additionalPointsToLabel = additionalPointsToLabel, 
          additionalLabelColor = additionalLabelColor)

})

#' volcanoPlot for an object resulting from \code{limma2Groups}
#' @rdname volcanoPlot-methods
#' @export
setMethod("volcanoPlot",
    signature(x = "limma", 
        y = "missing",
        pointLabels = "character"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10, smoothScatter = TRUE, 
        xlab = NULL, ylab = NULL, main = NULL, sub = NULL, newpage = TRUE, 
        additionalPointsToLabel = NULL, additionalLabelColor = "red"){
      
      logRatio <- as.matrix(x@MArrayLM$coef)[, 2]
      pValue <- as.matrix(x@MArrayLM$p.value)[, 2]
      
      volcanoplotter(logRatio = logRatio, pValue = pValue, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, logTransformP = TRUE, # p values  
          smoothScatter = smoothScatter, 
          xlab = xlab, ylab = ylab, main = main, sub = sub,
          newpage = newpage, additionalPointsToLabel = additionalPointsToLabel, 
          additionalLabelColor = additionalLabelColor)      
})


#' volcanoPlot for arbitrary numeric vectors 
#' containing log ratio values and p values respectively
#' @rdname volcanoPlot-methods
#' @export
setMethod("volcanoPlot",
    signature(x = "numeric", 
        y = "numeric",
        pointLabels = "character"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10,
        smoothScatter = TRUE, xlab = NULL, ylab = NULL,
        main = NULL, sub = NULL, newpage = TRUE, 
        additionalPointsToLabel = NULL, additionalLabelColor = "red"){
      
      if ((length(x) != length(y)) | (length(x) != length(pointLabels)))
        stop("'x', 'y' and 'pointLabels' should have equal length")
      
      volcanoplotter(logRatio = x, pValue = y, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, smoothScatter = smoothScatter, 
          xlab = xlab, ylab = ylab, main = main, sub = sub,
          newpage = newpage, additionalPointsToLabel = additionalPointsToLabel, 
          additionalLabelColor = additionalLabelColor)
})

#' volcanoPlot for arbitrary numeric vectors containing 
#' log ratio values and p values respectively
#' @rdname volcanoPlot-methods
#' @export
setMethod("volcanoPlot",
    signature(x = "numeric", 
        y = "numeric",
        pointLabels = "missing"),
    function(x, y, pointLabels, topPValues = 10, 
        topLogRatios = 10,
        smoothScatter = TRUE, xlab = NULL, ylab = NULL,
        main = NULL, sub = NULL, newpage = TRUE, 
        additionalPointsToLabel = NULL, additionalLabelColor = "red"){
      
      if (length(x) != length(y))
        stop("'x' and 'y' should have equal length")
      
      if (is.null(names(x))){
        if (is.null(names(y))){
          stop(paste("nor 'x' nor 'y' have names that can be used to use as default 'pointLabels'\n",
                     "please make sure either 'x' or 'y' has names or, alternatively, \n",
                     "explicitly provide a 'pointLabels' argument", sep = "")) 
        } else {
          pointLabels <- names(y)
        } 
      } else {
        pointLabels <- names(x)
      }
      
      
      volcanoplotter(logRatio = x, pValue = y, 
          pointLabels = pointLabels, topPValues = topPValues,
          topLogRatios = topLogRatios, smoothScatter = smoothScatter, 
          xlab = xlab, ylab = ylab, main = main, sub = sub,
          newpage = newpage, additionalPointsToLabel = additionalPointsToLabel, 
          additionalLabelColor = additionalLabelColor)
    })

#' Workhorse function for the different volcanoPlot methods
#' 
#' Workhorse function for the different volcanoPlot methods. A volcano plot 
#'   is a graph that allows to simultaneously assess the P values (statistical 
#'   significance) and log ratios (biological difference) of differential 
#'   expression for the given genes.
#' @param logRatio numeric vector of log ratios
#' @param pValue numeric vector of P values
#' @param pointLabels Labels for points on the volcano plot that are interesting
#'     taking into account both the x and y dimensions; typically this is a
#'     vector of gene symbols; most methods can access the gene symbols directly from 
#'     the object passed as 'x' argument; the argument allows for custom labels if
#'     needed
#' @param topPValues top n points that will be included in the points to label based
#'     on their low P Values
#' @param topLogRatios top n points that will be included in the points to label based
#'     on their high absolute values of the log ratio
#' @param logTransformP if \code{TRUE} (default) -log10(pValue) is used for the plot 
#'     instead of the raw P values
#' @param smoothScatter use color saturation to indicate dots that are in densely
#'     populated regions of the graph; defaults to \code{TRUE}
#' @param xlab label for the x axis (string)
#' @param ylab label for the y axis (string)
#' @param main main title for the graph (string)
#' @param sub subtitle for the graph (string)
#' @param newpage should the graph be drawn to a new grid page? Defaults to
#'     \code{TRUE}. This argument is useful for including several volcano plots 
#'     in one layout.
#' @param additionalPointsToLabel Entrez IDs of genes of interest, that will be highlighted on the plot; the color of highlighting is determined
#'     by the 'additionalLabelColor' argument.
#' @param additionalLabelColor Color used to highlight the 'additionalPointsToLabel'; defaults to "red"
#' @return a volcanoplot is drawn to the current device
#' @author Tobias Verbeke
#' @keywords dplot
#' @importFrom grid grid.newpage plotViewport pushViewport 
#'  textGrob gpar unit grobWidth convertHeight dataViewport 
#'  grid.pretty xaxisGrob grid.yaxis 
#'  editGrob gEditList gEdit grid.draw grid.points grid.text
#'  current.viewport
#' @importFrom grDevices densCols
#' @export
volcanoplotter <- function(logRatio, pValue, pointLabels,
    topPValues = 10, topLogRatios = 10, logTransformP = TRUE,
    smoothScatter = TRUE, xlab = NULL, ylab = NULL, main = NULL, 
    sub = NULL, newpage = TRUE, additionalPointsToLabel = NULL, additionalLabelColor = "red"){
  ### checks                      
  if (!is.numeric(logRatio))
    stop("'logRatio' should be numeric")
  if (!is.numeric(pValue))
    stop("'pValue' should be numeric")
  # test if pValue is between 0 and 1
    if (any(pValue < 0 | pValue > 1))
      stop("'pValue' should be >= 0 and <= 1") # prevent from being already on log scale
  
  pVals <- if (logTransformP) -log10(pValue) else pValue 
  
  ### compute which points to label on the graph
  topLR <- order(abs(logRatio), decreasing = TRUE)[seq(length.out = topLogRatios)] 
  topP <- order(pVals, decreasing = TRUE)[seq(length.out = topPValues)]#logTransformP)[seq(length.out = topPValues)]
  pointsToLabel <- union(topP, topLR)
  
  pointsToLabel <- union(pointsToLabel, which(names(logRatio) %in% additionalPointsToLabel))
  if(is.null(additionalPointsToLabel)){
    colPointsToLabel <- rep("black", length(pointsToLabel))
  } else{
    if(is.null(names(logRatio))){
      stop("labeling additional points required a named vector for the logRatios")
    }
    else{
      colPointsToLabel <- ifelse(names(logRatio)[pointsToLabel] %in% additionalPointsToLabel, additionalLabelColor, "black")
    }
  }
  
  ### set up graph
  if (newpage)
    grid.newpage()
  pvp <- plotViewport(c(5, 6, 5, 3))
  pushViewport(pvp)
  
  tg <- textGrob(label = pointLabels[pointsToLabel], 
      x = unit(logRatio[pointsToLabel], "native"), 
      y = unit(pVals[pointsToLabel], "native"), 
      gp = gpar(cex = 0.65, col = colPointsToLabel))

  # compute maximum grobwidth
  
  maxLabelWidth <- max(grobWidth(tg))
  nMaxLabelWidth <- convertHeight(maxLabelWidth, "native", 
      valueOnly = TRUE)
  
  dvp <- dataViewport(xscale = range(logRatio) + c(-nMaxLabelWidth/2.2, nMaxLabelWidth/2.2),
      yscale = range(pVals, na.rm = TRUE))
  
  pushViewport(dvp)
  
  atPositionsY <-  (current.viewport()$yscale)
  # atPositionsX <- grid.pretty(current.viewport()$xscale)
  xa <- xaxisGrob(name = "xa")# , at = atPositionsX)# , vp = pvp)
  
  grid.yaxis(name = "ya", at = atPositionsY, label = signif(10^(-atPositionsY), 1))
  
  # move the x axis down a bit
  moveUnit <- unit(-0.5, "char")
  xa <- editGrob(xa, edits = gEditList(
          gEdit("major", y = moveUnit),   # defaults to 0npc
          gEdit("ticks", y0 = moveUnit),  # defaults to 0npc
          gEdit("ticks", y1 = unit(-0.5, "lines") + moveUnit),   # defaults to -0.5lines
          gEdit("labels", y = unit(-1.5, "lines") + moveUnit)))  # defaults to -1.5lines 
  grid.draw(xa)
  
  # box(bty = "l", lwd = 1.5)
  dotColors <- if (smoothScatter){ 
        densCols(x = logRatio[-pointsToLabel], y = pVals[-pointsToLabel])
      } else {
        "#9ECAE1" # brewer.pal(9, "Blues")[4]
      } 
  grid.points(x = unit(logRatio[-pointsToLabel], "native"),
      y = unit(pVals[-pointsToLabel], "native"),
      pch = 20, gp = gpar(col = dotColors))
  
  grid.draw(tg)
  
  if (!is.null(main)){
    grid.text(label = main, y = unit(1, "npc") + unit(2, "lines"),
        gp = gpar(fontface = "bold"))
  }
  if (!is.null(sub)){
    grid.text(label = sub, y = unit(-4.25, "lines") + 0.5 * moveUnit)
  }
  if (!is.null(xlab)){
    grid.text(label = xlab, y = unit(-3, "lines") + 0.5 * moveUnit)
  }
  if (!is.null(ylab)){
    grid.text(label = ylab, x = unit(-4.5, "lines"), rot = 90)
  }          
}
