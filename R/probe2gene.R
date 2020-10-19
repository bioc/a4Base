#' Translate Affymetrix probeset IDs into gene symbols
#' 
#' Auxiliary function for (currently) spectralMap
#' allowing the conversion of Affy probeset IDs to gene symbols
#' @param probesetIds Affymetrix probeset IDs
#' @param chipPkg string indicating the annotation package for the chip
#' @return Vector containing the respective gene symbols
#' @author Tobias Verbeke
#' @seealso \code{\link{spectralMap}}, \code{\link[a4Classif]{lassoClass}}, ...
#' @examples
#' if (require(ALL)){
#' 	data(ALL, package = "ALL")
#' 	chip <- annotation(ALL)
#' 	chipAnnotationPkg <- paste(chip, "db", sep = ".")
#' 	res <- probe2gene(featureNames(ALL), chipAnnotationPkg)
#' 	head(res)
#' }
#' @keywords manip
#' @importFrom annaffy aafSymbol getText
probe2gene <- function(probesetIds, chipPkg){
  # DB-based annotation
  featureSymbols <- aafSymbol(probesetIds, chipPkg)
  return(getText(featureSymbols))
}
