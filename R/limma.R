#' Wrapper for the limma function for the comparison of two groups 
#' (two factor levels)
#' 
#' @param object object of class ExpressionSet
#' @param group string indicating the variable defining the two groups
#'               to be compared
#' @param probe2gene logical; if \code{TRUE} Affymetrix probeset IDs are translated
#' into gene symbols; if \code{FALSE} no such translation is done
#' @return 
#' S4 object of class 'limma' with the following two components:
#' \item{MArrayLM}{S4 object of class MArrayLM as returned by the limma
#' 			function of the limma package}
#' \item{geneSymbols}{character vector of gene symbols; this slot is only
#' 	populated if \code{probe2gene=TRUE} (and if the ExpressionSet object
#' 	is appropriately annotated by \code{addGeneInfo} for gene symbols to
#' 	be extracted)}
#' @author Tobias Verbeke and Willem Talloen
#' @note A 'topTable' method is defined for 'limma' objects.
#' @keywords models
#' @importFrom Biobase pData featureData
#' @importFrom stats model.matrix
#' @importFrom limma lmFit eBayes
#' @importFrom methods new
#' @export
limmaTwoLevels <- function(object, group, probe2gene = TRUE){
  f <- factor(pData(object)[, group])[,drop = TRUE]
  if (nlevels(f) != 2)
    stop("Use 'limmaTwoGroups' only with a 'group' variable having two group levels")
  
  design <- model.matrix(~ f)
  fit <- lmFit(object, design)
  fit <- eBayes(fit)
  res <- if (probe2gene){
        new("limma", MArrayLM = fit, geneSymbols = featureData(object)$`SYMBOL`)
      } else {
        fit
      }
  # use gene symbols from ExpressionSet and not the ones
  # that are in (the third column of) limmaObj@MArrayLM$genes
  
  return(res) 
}

#' Wrapper for the limma function for the comparison of two groups 
#' (two factor levels)
#' 
#' @param covariable string indicating the variable defining the continuous covariate
#' @inheritParams limmaTwoLevels
#' @importFrom Biobase pData featureData
#' @importFrom limma lmFit eBayes
#' @importFrom methods new
limmaReg <- function(object, covariable, probe2gene = TRUE){
	
	covar <- pData(object)[, covariable]
	tfidx <- !is.na(covar)
	eset2 <- object[,tfidx]
	covar <- pData(eset2)[, covariable]

	design <- cbind(mean = 1, slope = covar)
	fit <- lmFit(exprs(eset2), design)
	fit <- eBayes(fit)

	res <- if (probe2gene){
				new("limma", MArrayLM = fit, geneSymbols = featureData(object)$`SYMBOL`)
			} else {
				fit
			}
	# use gene symbols from ExpressionSet and not the ones
	# that are in (the third column of) limmaObj@MArrayLM$genes
	
	return(res) 
}

# check with Willem; no default coef (with message) / default n = 10
# coef = 2 because we are not interested whether the intercept is significant 
# but whether group 2 is significantly different from group 1

#' @param eb subset of \code{fit} containing
#' Empirical Bayesian estimates, so 
#' columns: 't', 'p-value' and 'lods' by default.
#' For expert use only.
#' @inheritParams limma::topTableF
#' @importFrom limma eBayes topTableF
#' @rdname topTable-methods
#' @export
setMethod("topTable", "limma",
    function(fit, n = 10, coef = 2, genelist = fit$genes, 
        eb = fit[c("t", "p.value", "lods")], adjust.method = "BH",
        sort.by = "B", resort.by = NULL, p.value = 1, lfc = 0){
      
      fit <- fit@MArrayLM
      
      ### from limma:::topTable
      if (length(coef) > 1) {
        coef <- unique(coef)
        if (length(fit$coef[1, coef]) < ncol(fit)) 
          fit <- eBayes(fit[, coef])
        if (sort.by == "B") 
          sort.by <- "F"
        return(topTableF(fit, number = n, genelist = genelist, 
                adjust.method = adjust.method, sort.by = sort.by, 
                p.value = p.value))
      }
      fit <- unclass(fit)
	  limma:::.topTableT(fit = fit[c("coefficients", "stdev.unscaled")], 
          coef = coef, number = n, genelist = fit$genes, A = fit$Amean, 
          eb = eb, adjust.method = "BH", 
          sort.by = sort.by, resort.by = resort.by, p.value = p.value, 
          lfc = lfc)
})

### redefine as well for MArrayLM objects

#' @export
#' @importFrom limma eBayes topTableF
#' @rdname topTable-methods
setMethod("topTable", "MArrayLM",
    function(fit, n, coef = 2, genelist = fit$genes, 
        eb = fit[c("t", "p.value", "lods")], adjust.method = "BH",
        sort.by = "B", resort.by = NULL, p.value = 1, lfc = 0){
      
      # coef = 2 because we are not interested whether the intercept is significant 
      # but whether group 2 is significantly different from group 1

      ### from limma:::topTable
      
      if (length(coef) > 1) {
        coef <- unique(coef)
        if (length(fit$coef[1, coef]) < ncol(fit)) 
          fit <- eBayes(fit[, coef])
        if (sort.by == "B") 
          sort.by <- "F"
        return(topTableF(fit, number = n, genelist = genelist, 
                adjust.method = adjust.method, sort.by = sort.by, 
                p.value = p.value))
      }
      fit <- unclass(fit)
	  limma:::.topTableT(
			fit = fit[c("coefficients", "stdev.unscaled")], 
          coef = coef, number = n, genelist = fit$genes, A = fit$Amean, 
          eb = eb, adjust.method = "BH", 
          sort.by = sort.by, resort.by = resort.by, p.value = p.value, 
          lfc = lfc
  	)
    })


