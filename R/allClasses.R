# TODO: Add comment
# 
# Author: tobias
###############################################################################


setClass("limma",
    representation = list(MArrayLM = "MArrayLM",
        geneSymbols = "character")
)

### S4 class for output of computeLogRatio function
#' Class "ExpressionSetWithComputation"
#' 
#' This class adds statistical information to the exprs of the ExpressionSet
#' as well as descriptive information to the pData of the ExpressionSet
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{new("ExpressionSetWithComputation", assayData, phenoData, featureData, experimentData, annotation, exprs, ...)}.
#' @section Extends:
#' \itemize{
#' \item{Class \code{\linkS4class{ExpressionSet}}, directly.}
#' \item{Class \code{\linkS4class{eSet}}, by class "ExpressionSet", distance 2.}
#' \item{Class \code{\linkS4class{VersionedBiobase}}, by class "ExpressionSet", distance 3.}
#' \item{Class \code{\linkS4class{Versioned}}, by class "ExpressionSet", distance 4.}
#' }
#' @section Methods:
#' No methods defined with class "ExpressionSetWithComputation" in the signature.
#' @slot assayData Object of class \code{"AssayData"}
#' @slot phenoData Object of class \code{"AnnotatedDataFrame"}
#' @slot featureData Object of class \code{"AnnotatedDataFrame"}
#' @slot experimentData Object of class \code{"MIAME"}
#' @slot annotation Object of class \code{"character"}
#' @slot .__classVersion__ Object of class \code{"Versions"}
#' @author Tobias Verbeke
#' @seealso \code{\linkS4class{ExpressionSet}}, \code{\link{computeLogRatio}}
#' @keywords classes
#' @export
setClass("ExpressionSetWithComputation",
		contains = "ExpressionSet")


