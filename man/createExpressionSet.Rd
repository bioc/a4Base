% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exprSetManipulations.R
\name{createExpressionSet}
\alias{createExpressionSet}
\title{combine gene expression and phenotype data onto a ExpressionSet object}
\usage{
createExpressionSet(
  exprs = matrix(nrow = 0, ncol = 0),
  phenoData = AnnotatedDataFrame(),
  varMetadata = NULL,
  dimLabels = c("rowNames", "colNames"),
  featureData = NULL,
  experimentData = MIAME(),
  annotation = character(0),
  changeColumnsNames = TRUE,
  ...
)
}
\arguments{
\item{exprs}{gene expression matrix}

\item{phenoData}{phenotype data associated with exprs columns, as a matrix or data.frame}

\item{varMetadata}{optional metadata on phenotype data}

\item{dimLabels}{see \code{\linkS4class{ExpressionSet}}}

\item{featureData}{see \code{\linkS4class{ExpressionSet}}}

\item{experimentData}{see \code{\linkS4class{ExpressionSet}}}

\item{annotation}{see \code{\linkS4class{ExpressionSet}}}

\item{changeColumnsNames}{Change exprs columns names -- see details}

\item{...}{\code{\dots}}
}
\value{
An object of class ExpressionSet
}
\description{
Basically a wrapper for \code{new('ExpressionSet',...)}, this function gathers gene
expression and phenotype data, after having checked their compatibility.
}
\details{
If \code{changeColumnsNames} is \code{TRUE}, then the procedure is the following: first one checks if phenoData contains a column named 'colNames'. If so, content will be used to rename exprs colums. On the other case, one uses combinations of phenoData columns to create new names. In any case, old columns names
 are stored within a column named 'oldcolnames' in the pData.
}
\examples{
# simulate expression data of 10 features (genes) measured in 4 samples
x <- matrix(rnorm(40), ncol = 4)
colnames(x) <- paste("sample", 1:4, sep = "_")
rownames(x) <- paste("feature", 1:10, sep = "_")
# simulate a phenodata with two variables
ToBePheno <- data.frame(Gender = rep(c('Male','Female'), 2), 
		Treatment = rep(c('Trt','Control'), each=2))
rownames(ToBePheno) <- paste("sample", 1:4, sep = "_")
eset <- createExpressionSet(exprs = x, phenoData = ToBePheno)
}
\seealso{
\code{\linkS4class{ExpressionSet}}
}
\author{
Eric Lecoutre
}
\keyword{data}
