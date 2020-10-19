#' Multiple regression using the Lasso algorithm as implemented in the glmnet package
#' 
#' Multiple regression using the Lasso algorithm as implemented in the glmnet package.
#'  This is a theoretically nice approach to see which combination of genes predict best 
#'  a continuous response. Empirical evidence that this actually works with high-dimensional
#'  data is however scarce.
#' @param object object containing the expression measurements; currently the
#'  only method supported is one for ExpressionSet objects
#' @param covariate character string indicating the column containing 
#'  the continuous covariate.
#' @return object of class \code{glmnet}
#' @references 
#' Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
#'     Microarrays, Chapman \& Hall/CRC, pp. 211.
#' @author Willem Talloen
#' @seealso \code{\link[a4Classif]{lassoClass}}
#' @examples 
#' if (require(ALL)){
#'  data(ALL, package = "ALL")
#'  ALL <- addGeneInfo(ALL)
#'  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'  resultLasso <- lassoReg(object = ALL[1:100,], covariate = "age")
#'  plot(resultLasso, label = TRUE,
#' 	   main = "Lasso coefficients in relation to degree of penalization.")
#'   featResultLasso <- topTable(resultLasso, n = 15)
#' }
#' @importFrom Biobase pData exprs featureData
#' @importFrom glmnet glmnet
#' @export
lassoReg <- function(object, covariate){
  covariateVector <- pData(object)[, covariate]
  if (!is.numeric(covariateVector))
	  stop("The argument 'covariate' needs to refer to CONTINUOUS variable
		from the 'phenoData(object)'")

  object <- object[,!is.na(covariateVector)]
  covariateVector <- pData(object)[,covariate]

  fit <- glmnet(t(exprs(object)), covariateVector, family="gaussian", alpha = 1)
  fit$featureData <- pData(featureData(object))
  return(fit) 
}




