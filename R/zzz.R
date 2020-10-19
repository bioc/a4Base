#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname){
  options(error = NULL)
  packageStartupMessage(paste("\na4Base version ", packageDescription("a4Base")$Version, 
          "\n", sep = ""))
}
