
#' Generic function to compute the proportion of differentially expressed
#'   genes that are present
#' @param object object of class \code{propDEgene}
#' @param ... further arguments for the method (currently none implemented)
#' @return numeric of length one giving the proportion of differentially expressed genes
#' @author Willem Talloen and Tobias Verbeke
#' @keywords htest
#' @exportMethod propDEgenes
setGeneric("propDEgenes", function(object, ...){
      standardGeneric("propDEgenes")
    })

#' Generic function to compute the proportion of differentially expressed
#'   genes that are present
#' @section Methods:
#' \describe{
#' limma
#' \item{object = "limma"}{propDEgenes method for a limma object}
#' numeric
#' \item{object = "numeric"}{propDEgenes method for a numeric vector, i.e. a vector
#'  of P Values}
#'}
#' @param object object of class \code{propDEgene}
#' @param ... further arguments for the method (currently none implemented)
#' @return numeric of length one giving the proportion of differentially expressed genes
#' @author Willem Talloen and Tobias Verbeke
#' @keywords htest
#' @docType methods
#' @name propDEgenes-methods
#' @aliases propDEgenes,limma
#' @aliases propDEgenes,numeric
NULL

#' @export
#' @rdname propDEgenes-methods
setMethod("propDEgenes", "limma",
    function(object, ...){
      
      nGenes <- length(object@geneSymbols)
      pValue <- a4Base::topTable(object, n = nGenes)$P.Value
      
      propdegenescalculation(pValue = pValue)
})

#' @export
#' @rdname propDEgenes-methods
setMethod("propDEgenes", "numeric",
    function(object, ...){
      
      propdegenescalculation(pValue = object)
      
    })

#' Estimation of proportion of differentially expressed genes
#' 
#' Estimation of proportion of differentially expressed genes.
#'  This estimation is based on a histogram of the p-values. More specifically,
#'   based on the horizontal line representing a uniform distribution 
#'   based on the p value distribution between 0.5 and 1. This represents
#'   the hypothetical p value distribution arising just by chance.
#'  All genes with small p-values above this line reflect
#'   the expected number of differentially expressed genes not by chance.
#' @param pValue a vector of p-values
#' @return proportion of differential genes
#' @examples
#' if (require(ALL)){
#' 	data(ALL, package = "ALL")
#' 	ALL <- addGeneInfo(ALL)
#' 	ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#' 	
#' 	tTestResult <- tTest(ALL, "BTtype")
#' 	histPvalue(tTestResult[,"p"], addLegend = TRUE)
#' 	propDEgenesRes <- propDEgenes(tTestResult[,"p"])
#' }
#' @seealso \code{\link{histPvalue}}
#' @author Willem Talloen and Tobias Verbeke
#' @export
propdegenescalculation <- function(pValue){
  NbDEgenes <- length(pValue) - (sum(pValue > 0.5)*2)
  return( round((100/ length(pValue) )* NbDEgenes, 1) )
}

