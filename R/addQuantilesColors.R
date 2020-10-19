
#' Compute quantiles for plotGeneDE function
#' 
#' Compute quantiles on mean expression level for plotGeneDE function. 
#' Colors of bars in the plot could 
#' then be allocated using buckets defined by those quantiles.
#' @details Number of computed quantiles is equal to (ngroups - 1).
#' @param e ExpressionSet object to use for computation
#' @param ngroups Number of groups to be created
#' @return The ExpressionSet object e is returned, 
#' with a new column called colorsQuantilesVector in its slot featureData
#' @seealso \code{\link{plotLogRatio}}
#' @author Eric Lecoutre
#' @examples
#' if (require(ALL)){
#'   data(ALL, package = "ALL")
#'   ALLQ <- addQuantilesColors(ALL)
#'   fData(ALLQ)
#' }
#' @importFrom Biobase exprs `fData<-`
#' @importFrom stats quantile
#' @export
addQuantilesColors <- function(e,ngroups=3)
{
#D ngroups=5
### will add quantiles levels in ExpressionSet featureData slot
  stopifnot(is(e,'ExpressionSet'))
  means <- apply(exprs(e),1,mean,na.rm=TRUE)
  quantiles <- quantile(means,probs=seq(0,1,length=ngroups+1))
  colorsVector <- cut(means,quantiles,include.lowest=TRUE,labels=FALSE)
  if ('colorsQuantilesVector' %in% colnames(fData(e))) 
    fData(e)[,'colorsQuantilesVector'] <- colorsVector else 
    fData(e) <- data.frame(cbind(as.matrix(fData(e)),colorsQuantilesVector=colorsVector))
    invisible(e)    
  }
