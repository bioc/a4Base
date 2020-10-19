#' Plot the Distribution of P Values
#' 
#' This function displays the distribution of the p values using
#' a histogram; the horizontal line represents a uniform distribution 
#' based on the p value distribution between 0.5 and 1. This represents
#' the hypothetical p value distribution arising just by chance.
#' This uniform distribution is used to estimate the proportion of differentially
#' expressed genes.
#' @param object either a numeric vector of p-values,
#'  or an object of class \code{tTest}, \code{limma} or \code{MArrayLM}
#' @param ... further arguments passed to the method
#' @return The histogram is displayed on the current device. 
#' @references
#' Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
#'  Microarrays, Chapman \& Hall/CRC, p. 253.
#' @author Willem Talloen and Tobias Verbeke
#' @examples 
#' if (require(ALL)){
#'  data(ALL, package = "ALL")
#'  ALL <- addGeneInfo(ALL)
#'  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'  tTestResult <- tTest(ALL, "BTtype")
#'  histPvalue(tTestResult[,"p"], addLegend = TRUE)
#'  propDEgenesRes <- propDEgenes(tTestResult[,"p"])  
#' }
#' @keywords dplot
#' @export 
#' @rdname histPvalue
setGeneric("histPvalue", function(object, ...){
      standardGeneric("histPvalue")
})    

#' @export 
#' @inheritParams histPvalue
#' @importFrom limma topTable
#' @rdname histPvalue
setMethod("histPvalue", "limma",
    function(object, ...){
    
    nGenes <- length(object@geneSymbols)
    
    # currently default of coef = 2 is OK (as only limmaTwoLevels generates an object of class 'limma')
    pValue <- topTable(object, coef = 2, n = nGenes)$P.Value
    
    histpvalueplotter(pValue = pValue, ...)      
})

#' @param index of the coefficient for which the p values should be plotted;
#'  only applies to the MArrayLM method
#' @inheritParams histPvalue
#' @export
#' @importFrom limma topTable
#' @rdname histPvalue
setMethod("histPvalue", "MArrayLM",
    function(object, coef, ...){
      
      if (missing(coef))
        stop("Please specify a 'coef' argument to select a coefficient for which the P values should be displayed.")
      
      pValue <- topTable(object, coef = coef, n = nrow(object))$P.Value
      
      histpvalueplotter(pValue = pValue, ...)      
    })

#' @inheritParams histPvalue
#' @export
#' @rdname histPvalue
setMethod("histPvalue", "numeric",
    function(object, ...){
    histpvalueplotter(pValue = object, ...)      
})


#' Workhorse function for the histPvalue function
#' 
#' Workhorse function for the histPvalue function. 
#' This function displays the distribution of the p values using
#'   a histogram; the horizontal line represents a uniform distribution 
#'   based on the p value distribution between 0.5 and 1. This represents
#'   the hypothetical p value distribution arising just by chance.
#'   This uniform distribution is used to estimate the proportion of differentially
#'   expressed genes.
#' @param pValue numeric vector of p values
#' @param addLegend logical; should a legend be added (TRUE) or not (FALSE; default
#' @param xlab label for the x axis; defaults to NULL (no label)
#' @param ylab label for the y axis; defaults to NULL (no label
#' @param main main title for the plot; if NULL (default) no main title is displayed
#' @param ... further arguments for the \code{hist} call; currently none are used
#' @return no returned value, a plot is drawn to the current device.
#' @seealso \code{\link{histPvalue}}, \code{\link{propdegenescalculation}}
#' @author Willem Talloen and Tobias Verbeke
#' @examples 
#'  if (require(ALL)){
#'    data(ALL, package = "ALL")
#'    ALL <- addGeneInfo(ALL)
#'    ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'    tTestResult <- tTest(ALL, "BTtype")
#'    histPvalue(tTestResult[,"p"], addLegend = TRUE, xlab = "Adjusted P Value")
#'    histPvalue(tTestResult[,"p"], addLegend = TRUE, main = "Histogram of Adjusted P Values")
#'    propDEgenesRes <- propDEgenes(tTestResult[,"p"])
#'  }
#' @importFrom graphics hist abline legend
#' @export
histpvalueplotter <- function(pValue, addLegend = FALSE, xlab = NULL, 
	ylab = NULL, main = NULL, ...){
  
  mainTitle <- if (is.null(main)) "" else  main
  
  histOutput <- hist(pValue, 50, col = "skyblue", main = mainTitle, 
	xlab = xlab, ylab = ylab, ...)
  lengthHist <- length(histOutput$counts)
  meanNonDE <- mean(histOutput$counts[(lengthHist/2):lengthHist])
  abline(h = meanNonDE, col = 'goldenrod', lwd = 2)
  
  if (addLegend){
    legend("topright", bty = "n", 
        legend = paste(propDEgenes(pValue), "% DE genes", sep = ""),
        text.col = "goldenrod", cex = 1.2)
  }
}
