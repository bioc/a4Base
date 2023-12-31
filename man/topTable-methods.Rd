% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/limma.R, R/tTest.R
\docType{methods}
\name{topTable,limma-method}
\alias{topTable,limma-method}
\alias{topTable,MArrayLM-method}
\alias{topTable-methods}
\alias{topTable,tTest}
\alias{topTable,fTest}
\alias{topTable,limma}
\alias{topTable,MArrayLM}
\alias{topTable,tTest-method}
\alias{topTable,fTest-method}
\title{Methods for topTable}
\usage{
\S4method{topTable}{limma}(
  fit,
  n = 10,
  coef = 2,
  genelist = fit$genes,
  eb = fit[c("t", "p.value", "lods")],
  adjust.method = "BH",
  sort.by = "B",
  resort.by = NULL,
  p.value = 1,
  lfc = 0
)

\S4method{topTable}{MArrayLM}(
  fit,
  n,
  coef = 2,
  genelist = fit$genes,
  eb = fit[c("t", "p.value", "lods")],
  adjust.method = "BH",
  sort.by = "B",
  resort.by = NULL,
  p.value = 1,
  lfc = 0
)

\S4method{topTable}{tTest}(fit, n)

\S4method{topTable}{fTest}(fit, n)
}
\arguments{
\item{fit}{object resulting from a classification or regression procedure}

\item{n}{number of features that one wants to extract from a table that
ranks all features according to their importance in the classification
or regression model; defaults to 10 for limma objects}

\item{coef}{column number or column name specifying which coefficient or contrast of the linear model is of interest. For \code{topTable}, can also be a vector of column subscripts, in which case the gene ranking is by F-statistic for that set of contrasts.}

\item{genelist}{data frame or character vector containing gene information.
  For \code{topTable} only, this defaults to \code{fit$genes}.}

\item{eb}{subset of \code{fit} containing
Empirical Bayesian estimates, so 
columns: 't', 'p-value' and 'lods' by default.
For expert use only.}

\item{adjust.method}{method used to adjust the p-values for multiple testing.  Options, in increasing conservatism, include \code{"none"}, \code{"BH"}, \code{"BY"} and \code{"holm"}.
  See \code{\link{p.adjust}} for the complete list of options. A \code{NULL} value will result in the default adjustment method, which is \code{"BH"}.}

\item{sort.by}{
    character string specifying which statistic to rank the genes by.
    Possible values for \code{topTable} are \code{"logFC"}, \code{"AveExpr"}, \code{"t"}, \code{"P"}, \code{"p"}, \code{"B"} or \code{"none"}.
    (Permitted synonyms are \code{"M"} for \code{"logFC"}, \code{"A"} or \code{"Amean"} for \code{"AveExpr"}, \code{"T"} for \code{"t"} and \code{"p"} for \code{"P"}.)
    Possible values for \code{topTableF} are \code{"F"} or \code{"none"}.
    \code{topTreat} accepts the same values as \code{topTable} except for \code{"B"}.
  }

\item{resort.by}{
    character string specifying statistic to sort the selected genes by in the output data.frame.  Possibilities are the same as for \code{sort.by}.
  }

\item{p.value}{cutoff value for adjusted p-values. Only genes with lower p-values are listed.}

\item{lfc}{
    minimum absolute log2-fold-change required.
    \code{topTable} and \code{topTableF} include only genes with (at least one) absolute log-fold-change greater than \code{lfc}.
    \code{topTreat} does not remove genes but ranks genes by evidence that their log-fold-change exceeds \code{lfc}.
  }
}
\description{
Methods for topTable. topTable extracts the top n most important features
for a given classification or regression procedure
}
\section{Methods}{

\describe{
glmnet
\item{fit = "glmnet", n = "numeric"}{glmnet objects are produced by \code{lassoClass} or \code{lassoReg}}
limma
\item{fit = "limma", n = "numeric"}{limma objects are produced by \code{limma2Groups}}
MarrayLM
\item{fit = "limma", n = "numeric"}{MarrayLM objects are produced by \code{lmFit} of the \code{limma package}}
pamClass
\item{fit = "pamClass", n = "numeric"}{pamClass objects are produced by \code{pamClass}}
rfClass
\item{fit = "rfClass", n = "numeric"}{rfClass objects are produced by \code{rfClass}}
tTest
\item{fit = "tTest", n = "numeric"}{tTest objects are produced by \code{tTest}}
fTest
\item{fit = "fTest", n = "numeric"}{fTest objects are produced by \code{fTest}}
}
}

\seealso{
\itemize{
\item{\code{\link[a4Core]{topTable-methods}} for: }{ glmnet, lognet and elnet}
\item{\code{\link[a4Classif]{topTable,pamClass-method}}} 
\item{\code{\link[a4Classif]{topTable,rfClass-method}}}
}
}
\keyword{manip}
\keyword{methods}
