#' computes replicates across a vector
#' 
#' Given a vector, returns the replicates in order
#' @param x character or numeric vector
#' @return numeric vector 
#' @author Henrique Dallazuanna
#' @references R-help mailing list
#' @seealso \code{\link{rle}}
#' @examples   
#' x <- c('a','b','a','a','b','a','c','c','c')
#'   data.frame(val=x,rep=replicates(x))
#' @keywords manip
#' @export
replicates <- function(x)
{
  x <- as.character(x)
  x[is.na(x)] <- ''
  replicate <- vector("numeric", length = length(x))
  replicate[order(x)] <- unlist(sapply(rle(sort(x))$lengths, seq_len))
  return(replicate)
}

