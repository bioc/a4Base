#
# Filtering on Intensity and Variance
# 
# Author: wtalloen
###############################################################################



#' Filtering on Intensity and Variance
#' 
#' Function to filter on intensity and variance as
#' typically used in gene expression studies
#' @param object ExpressionSet object
#' @param IntCutOff cut-off value used for the intensity filt
#' @param IntPropSamples proportion of samples used by the intensity
#'  filter; by default \code{IntPropSamples} is set to 0.25.
#' @param VarCutOff cut-off value used for the variance filter
#' @details 
#'  The intensity filter implies that (by default) the intensity
#'  levels must be greater than log2(100) in at least 25 percent
#'  of the samples.
#'  The variance filter requires that the features have an interquartile
#'  range (IQR) greater than 0.5. Note that the IQR is quite insensitive 
#'  to outliers such that genes with outlying expression values in a few 
#'  samples are excluded as long as their overall variation is small.
#' @return   
#' Object of class ExpressionSet containing only the features that
#' pass the variance and intensity filter.
#' @references 
#' Gentleman, R. et al. (2005). Bioinformatics and Computational Biology Solutions 
#'  using R and BioConductor, New York: Springer.
#' Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
#'  Microarrays, Chapman \& Hall/CRC, p. 128.
#' @examples 
#' if (require(ALL)){
#'   data(ALL, package = "ALL")
#'   fALL <- filterVarInt(ALL)
#'   fALL
#' }
#' @author Willem Talloen
#' @seealso \code{\link[genefilter]{pOverA}}, \code{\link[genefilter]{filterfun}}
#' @keywords manip
#' @export
filterVarInt <- function(object,
    IntCutOff = log2(100), # exclude low-expressed genes: genes with less than 6.6 on log2 scale 
    IntPropSamples = 0.25, # in more than 3/4 of the samples
    VarCutOff = 0.5){			 # exclude small-variance genes: genes with a InterQuartileRange smaller than 0.5
  f1 <- pOverA(IntPropSamples, IntCutOff)     
  f2 <- function(x) IQR(x) > VarCutOff   
  ff <- filterfun(f1, f2)
  selected <- genefilter(object, ff)
  esSel <- object[selected, ]
  return(esSel)
}
