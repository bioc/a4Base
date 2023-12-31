% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lassoReg.R
\name{lassoReg}
\alias{lassoReg}
\title{Multiple regression using the Lasso algorithm as implemented in the glmnet package}
\usage{
lassoReg(object, covariate)
}
\arguments{
\item{object}{object containing the expression measurements; currently the
only method supported is one for ExpressionSet objects}

\item{covariate}{character string indicating the column containing 
the continuous covariate.}
}
\value{
object of class \code{glmnet}
}
\description{
Multiple regression using the Lasso algorithm as implemented in the glmnet package.
 This is a theoretically nice approach to see which combination of genes predict best 
 a continuous response. Empirical evidence that this actually works with high-dimensional
 data is however scarce.
}
\examples{
if (require(ALL)){
 data(ALL, package = "ALL")
 ALL <- addGeneInfo(ALL)
 ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
 resultLasso <- lassoReg(object = ALL[1:100,], covariate = "age")
 plot(resultLasso, label = TRUE,
	   main = "Lasso coefficients in relation to degree of penalization.")
  featResultLasso <- topTable(resultLasso, n = 15)
}
}
\references{
Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
    Microarrays, Chapman \& Hall/CRC, pp. 211.
}
\seealso{
\code{\link[a4Classif]{lassoClass}}
}
\author{
Willem Talloen
}
