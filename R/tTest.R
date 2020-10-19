# 
# ordinary t-test
# 
# Author: wtalloen
###############################################################################

#' Use t Test to Compare Two Groups
#' 
#' Use a (modified) t test to compare two groups
#' @param object ExpressionSet object
#' @param groups string indicating the name of the variable of the phenoData
#' containing the group information
#' @param probe2gene logical; if \code{TRUE} Affymetrix probeset IDs are translated
#' into gene symbols; if \code{FALSE} no such translation is conducted
#' @return Object of class \code{"tTest"}, a data frame with the following columns
#' \item{gSymbol}{Gene Symbol}
#' \item{p}{p-value of the difference between the groups}
#' \item{logRatio}{Log ratio of the expression between the groups}
#' \item{pBH}{p-value of the difference between the groups,
#' with Benjamini-Hochberg multiplicity correction}
#' \item{tStat}{Student t-statistic of the different between groups}
#' @details For multiple testing the \code{mt.rawp2adjp} function of package
#' \code{multtest} is used.
#' @examples 
#' if (require(ALL)){
#'   data(ALL, package = "ALL")
#'   ALL <- addGeneInfo(ALL)
#'   ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'   tTestRes <- tTest(object = ALL,	groups = "BTtype", probe2gene = TRUE)
#'   volcanoPlot(tTestRes)  
#' }
#' @author Willem Talloen, Tobias Verbeke
#' @seealso \code{\link[genefilter]{rowttests}}
#' @keywords htest
#' @importFrom genefilter rowttests
#' @importFrom multtest mt.rawp2adjp
#' @importFrom Biobase featureData
#' @export
tTest <- function(object, groups, probe2gene = TRUE){

	# t-test for differential expression
	ttests <- rowttests(object, groups)
	pTtest <- data.frame(rownames(ttests), ttests[,'p.value'])
	statTtest <- data.frame(rownames(ttests), ttests[,'statistic'])
	
	# adjustment for multiple testing
	pAdjusted <- mt.rawp2adjp(ttests[, "p.value"], proc = c("BH"))
	pTtestBH <- pAdjusted$adjp[order(pAdjusted$index), "BH"]
	
	# log-ratio of differential expression
	labels <- as.numeric(factor(pData(object)[,groups]))-1
	logRatio <- rowMeans(exprs(object)[, labels == 1]) - rowMeans(exprs(object)[,labels == 0])
	
	if (probe2gene){
		gSymbol <- featureData(object)$`SYMBOL`
		if (is.null(gSymbol))
			stop("There is no variable named'SYMBOL' in the pData of the object.\n
				You may want to set the argument 'probe2gene' to FALSE (the default is TRUE)")
		
		pvalues <- data.frame(gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("gSymbol", "p", "logRatio", "pBH", "tStat")
	} else {
		pvalues <- data.frame(# gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("p", "logRatio", "pBH", "tStat")
	}
	
	class(pvalues) <- c("tTest", class(pvalues))
	return(pvalues)
}

fTest <- function(object, groups, probe2gene = TRUE, varEqual = FALSE){
	# t-test for differential expression
	ttests <- rowFtests(object, groups, var.equal = varEqual)
	pTtest <- data.frame(rownames(ttests), ttests[,'p.value'])
	statTtest <- data.frame(rownames(ttests), ttests[,'statistic'])
	
	# adjustment for multiple testing
	pAdjusted <- mt.rawp2adjp(ttests[, "p.value"], proc = c("BH"))
	pTtestBH <- pAdjusted$adjp[order(pAdjusted$index), "BH"]
	
	# log-ratio of differential expression
	labels <- as.numeric(factor(pData(object)[,groups]))-1
	logRatio <- rowMeans(exprs(object)[, labels == 1]) - rowMeans(exprs(object)[,labels == 0])
	
	if (probe2gene){
		gSymbol <- featureData(object)$`SYMBOL`
		pvalues <- data.frame(gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("gSymbol", "p", "logRatio", "pBH", "fStat")
	} else {
		pvalues <- data.frame(# gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("p", "logRatio", "pBH", "fStat")
	}
	
	class(pvalues) <- c("fTest", class(pvalues))
	return(pvalues)
}

#' @importFrom a4Core topTable
setOldClass("tTest")


#' Methods for topTable
#' 
#' Methods for topTable. topTable extracts the top n most important features
#' for a given classification or regression procedure 
#' \section{Methods}{
#' \describe{
#' glmnet
#' \item{fit = "glmnet", n = "numeric"}{glmnet objects are produced by \code{lassoClass} or \code{lassoReg}}
#' limma
#' \item{fit = "limma", n = "numeric"}{limma objects are produced by \code{limma2Groups}}
#' MarrayLM
#' \item{fit = "limma", n = "numeric"}{MarrayLM objects are produced by \code{lmFit} of the \code{limma package}}
#' pamClass
#' \item{fit = "pamClass", n = "numeric"}{pamClass objects are produced by \code{pamClass}}
#' rfClass
#' \item{fit = "rfClass", n = "numeric"}{rfClass objects are produced by \code{rfClass}}
#' tTest
#' \item{fit = "tTest", n = "numeric"}{tTest objects are produced by \code{tTest}}
#' fTest
#' \item{fit = "fTest", n = "numeric"}{fTest objects are produced by \code{fTest}}
#' }
#' }
#' @param fit object resulting from a classification or regression procedure
#' @param n number of features that one wants to extract from a table that
#' ranks all features according to their importance in the classification
#' or regression model; defaults to 10 for limma objects
#' @seealso \code{\link[a4Core]{topTable,glmnet-method}},
#' \code{\link[a4Core]{topTable,lognet-method}},
#' \code{\link[a4Core]{topTable,elnet-method}},
#' \code{\link[a4Classif]{pamClass,lognet-method}},
#' \code{\link[a4Classif]{rfClass,elnet-method}},
#' @keywords methods manip
#' @docType methods
#' @name topTable-methods
NULL

#' @docType methods
#' @rdname topTable-methods
#' @export
setMethod("topTable", "tTest",
    function(fit, n){
      head(fit, n = n)
})

setOldClass("fTest")


#' @rdname topTable-methods
#' @docType methods
#' @export
setMethod("topTable", "fTest",
		function(fit, n){
			head(fit, n = n)
})


tTest2 <- function(object, groups, probe2gene = TRUE){
	groups <- pData(object)[, groups]
	testData <- exprs(object)
	# t-test for differential expression
	ttestfun <- function(y) t.test(y ~ groups, var.equal = TRUE)
	ttests <- apply(testData, 1, ttestfun)
	rowttests(object, groups)
	pTtest <- data.frame(rownames(ttests), ttests[,'p.value'])
	statTtest <- data.frame(rownames(ttests), ttests[,'statistic'])
	
	# adjustment for multiple testing
	pAdjusted <- mt.rawp2adjp(ttests[, "p.value"], proc = c("BH"))
	pTtestBH <- pAdjusted$adjp[order(pAdjusted$index), "BH"]
	
	# log-ratio of differential expression
	labels <- as.numeric(factor(pData(object)[,groups]))-1
	logRatio <- rowMeans(exprs(object)[, labels == 1]) - rowMeans(exprs(object)[,labels == 0])
	
	if (probe2gene){
		gSymbol <- featureData(object)$`SYMBOL`
		pvalues <- data.frame(gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("gSymbol", "p", "logRatio", "pBH", "tStat")
	} else {
		pvalues <- data.frame(# gSymbol,
				ttests[,"p.value"],
				logRatio,
				pTtestBH,
				ttests[,"statistic"])[pAdjusted$index, ]
		colnames(pvalues) <- c("p", "logRatio", "pBH", "tStat")
	}
	
	class(pvalues) <- c("tTest", class(pvalues))
	return(pvalues)
}
