\name{a4palette}
\alias{a4palette}
\title{Utility function that defines a color palette for use in a4}
\description{
  Utility function that defines a color palette for use in a4.
}
\usage{
a4palette(n, alpha = 1, Janssen = FALSE)
}
\arguments{
  \item{n}{Number of color levels the palette should provide}
  \item{alpha}{alpha transparency level of the colors}
  \item{Janssen}{logical.  If \code{TRUE}, Janssen Pharmaceutical colors
   are used (with a maximum of 6 possible colors).}
}
\details{
  For n = 1, \code{"blue"} is returned; for n = 2 
  \code{c("red", "blue")} is returned; for n = 3 
  \code{c("red", "green", "blue"} is returned; for n = 4 
  \code{c("red", "green", "blue", "purple")} is returned and for n > 2, 
  the output of \code{rainbow(n)} is returned.
}
\value{
  a character vector of colors
}
\author{Steven Osselaer, Tobias Verbeke}
\seealso{\code{\link[grDevices]{rainbow}}}
\examples{
  op <- par(mfrow = c(2, 3))
  for (nGroups in 1:6)
    pie(rep(1, nGroups), a4palette(nGroups))
  par(op)
}
\keyword{dplot}

