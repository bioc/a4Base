% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addQuantilesColors.R
\name{addQuantilesColors}
\alias{addQuantilesColors}
\title{Compute quantiles for plotGeneDE function}
\usage{
addQuantilesColors(e, ngroups = 3)
}
\arguments{
\item{e}{ExpressionSet object to use for computation}

\item{ngroups}{Number of groups to be created}
}
\value{
The ExpressionSet object e is returned, 
with a new column called colorsQuantilesVector in its slot featureData
}
\description{
Compute quantiles on mean expression level for plotGeneDE function. 
Colors of bars in the plot could 
then be allocated using buckets defined by those quantiles.
}
\details{
Number of computed quantiles is equal to (ngroups - 1).
}
\examples{
if (require(ALL)){
  data(ALL, package = "ALL")
  ALLQ <- addQuantilesColors(ALL)
  fData(ALLQ)
}
}
\seealso{
\code{\link{plotLogRatio}}
}
\author{
Eric Lecoutre
}
