#' Generic function to draw a spectral map, according to JnJ Standards
#' @param object object of class ExpressionSet
#' @param groups string indicating the name of the column in the phenoData that
#'  defines the groups
#' @param ... further arguments to be passed to the methods
#' @author Tobias Verbeke
#' @return Object of class \code{plot.mpm}, i.e. the S3 output object of the \code{plot.mpm}
#' function of the \code{mpm} package
#' @note Coloring of groups on the spectralMap uses the a4 palette as produced
#'  by \code{a4palette}
#' @seealso \code{\link[mpm]{plot.mpm}}
#' @references 
#' Wouters, L., Goehlmann, H., Bijnens, L., Kass, S.U., Molenberghs, G.,
#'   Lewi, P.J. (2003). Graphical exploration of gene expression data: a
#'   comparative study of three multivariate methods. \emph{Biometrics}
#'   \bold{59}, 1131-1140.
#'   Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
#'     Microarrays, Chapman \& Hall/CRC, pp. 148 - 153.
#' @examples 
#' if (require(ALL)){
#' 	data(ALL, package = "ALL")
#' 	ALL <- addGeneInfo(ALL)
#' 	spectralMap(object = ALL, groups = "BT", legendPos = 'bottomright')
#' 	spectralMap(object = ALL, groups = "BT",
#' 	 plot.mpm.args = list(label.tol = 10, rot = c(-1, 1), sub = "", lab.size = 0.65,
#' 					dim = c(1,2), sampleNames = FALSE, zoom = c(1,5), col.size = 2, 
#' 					do.smoothScatter = TRUE))
#' 	spectralMap(object = ALL, groups = "BT",
#' 			plot.mpm.args = list(label.tol = 10, rot = c(-1, 1), sub = "", lab.size = 0.65,
#' 					dim = c(1,2), sampleNames = as.character(pData(ALL)$BT),
#' 					zoom = c(1,5), col.size = 2, do.smoothScatter = TRUE))
#' }
#' @keywords hplot
#' @exportMethod spectralMap
setGeneric("spectralMap", function(object, groups, ...){
      standardGeneric("spectralMap")
    })    

#' Methods for Function spectralMap according to JnJ Standards
#' 
#' Methods for spectralMap
#' @param object object of class ExpressionSet
#' @param groups string indicating the name of the column in the phenoData that
#'  defines the groups
#' @param makeLognormal boolean indicating whether one wants to exponentiate the
#' 	data to make them lognormally shaped (\code{TRUE}; the default) or not 
#' (\code{FALSE})
#' @param mpm.args list of arguments that can be passed to the \code{mpm} function
#' @param plot.mpm.args list of arguments that can be passed to the 
#' 	\code{plot.mpm} function that actually draws the plot
#' @param probe2gene boolean indicating whether one wants to display the gene symbols
#' 	for the labeled points (\code{TRUE}) or not (\code{FALSE}; the default)
#' @param addLegend Boolean indicating whether a legend for the colors of the dots should be added.
#' @param legendPos Specify where the legend should be placed. Typically either \code{topright},
#' @aliases spectralMap,ExpressionSet,character
#' @param ... further arguments to be passed to the methods,
#' currently not used.
#' @return the plot is returned invisibly 
#' @author Tobias Verbeke
#' @docType methods
#' @keywords methods hplot
#' @importFrom Biobase `pData<-` pData exprs annotation sampleNames
#' @importFrom mpm mpm plot.mpm
#' @importFrom graphics par legend
#' @importFrom stats na.omit
#' @name spectralMap-methods
#' @aliases volcanoPlot,ExpressionSet,character
NULL

#' @rdname spectralMap-methods
#' @export
setMethod("spectralMap",
    signature(object = "ExpressionSet", 
        groups = "character"),
    function(object, groups, makeLognormal = TRUE,
        mpm.args = list(row.weight = "mean", # mpmObject 
            col.weight = "constant", 
            logtrans = TRUE),
        plot.mpm.args = list(
            zoom = c(1,2),  # only these arguments are included that differ from plot.mpm defaults     
            label.tol = 10,  # please refer to ?plot.mpm for more information
            rot = c(-1, 1), 
            sub = "",
            lab.size = 0.85,
            col.group = pData(object)[, groups],
            # colors = c("orange1", "red", rainbow(length(unique(col.group)), start=2/6, end=4/6)),
            colors = c("wheat", # gene color (if no smoothScatter is used)
                "darkgrey", # color for genes considered to be outlying 
                a4palette(nlevels(pData(object)[, groups]))), # colors for the groups 
            col.size = 2,
            do.smoothScatter = TRUE),
        probe2gene = TRUE, addLegend = TRUE, legendPos = "topleft", ...){
      
      if (!is.list(mpm.args))
        stop("'mpm.args' should be a list of arguments to pass to the 'mpm' function")
      
      if (!is.list(plot.mpm.args))
        stop("'plot.mpm.args' should be a list of arguments to pass to the 'plot.mpm' function")
      
      if (length(groups) > 1){
        stop("'groups' should be a string (character vector of length one)")
      }

      if (any(is.na(pData(object)[, groups]))){
        stop("'groups' variable contains missing values")
      }
      
      if (!is.factor(pData(object)[, groups])){
        warning("'groups' should refer to a factor variable \n
                The variable has been transformed into factor variable")
        pData(object)[, groups] <- factor(pData(object)[, groups])
      }
      
      expressionData <- exprs(object)
      chip <- annotation(object)
      chipAnnotationPkg <- paste(chip, "db", sep = ".")
      
      mpmInput <- if (makeLognormal){ 
            data.frame(rownames(expressionData), 2^expressionData) 
          } else {
            data.frame(rownames(expressionData), expressionData)
          }
      mpmInput <- na.omit(mpmInput)
      mpm.args$data <- mpmInput
      
      # compute the projection
      plot.mpm.args$x <- do.call(mpm, mpm.args)    
      
      # adjust sample names (to escape the constraints of data frame column names)
      #   otherwise 'X' will have been prepended by the mpm function and displayed as such
      plot.mpm.args$x$col.names <- sampleNames(object)
      
      plot.mpm.args$zoom <- if (is.null(plot.mpm.args$zoom)) c(1,2) 
          else plot.mpm.args$zoom
      plot.mpm.args$label.tol <- if (is.null(plot.mpm.args$label.tol)) 10
          else plot.mpm.args$label.tol
      plot.mpm.args$rot <- if (is.null(plot.mpm.args$rot)) c(-1, 1)
          else plot.mpm.args$rot
      plot.mpm.args$sub <- if (is.null(plot.mpm.args$sub)) ""
          else plot.mpm.args$sub
      plot.mpm.args$lab.size <- if (is.null(plot.mpm.args$lab.size)) 0.85
          else plot.mpm.args$lab.size
      plot.mpm.args$col.group <- if (is.null(plot.mpm.args$col.group)) pData(object)[, groups]
          else plot.mpm.args$col.group
      plot.mpm.args$colors <- if (is.null(plot.mpm.args$colors)) 
            c("wheat", "black", a4palette(nlevels(pData(object)[, groups])))
          else
            plot.mpm.args$colors
      plot.mpm.args$col.size <- if (is.null(plot.mpm.args$col.size)) 2
          else plot.mpm.args$col.size
      plot.mpm.args$scale <- if (is.null(plot.mpm.args$scale)) "uvc"
          else plot.mpm.args$scale
      
      if (probe2gene){
        plot.mpm.args$labels <- pData(featureData(object))[plot.mpm.args$x$row.names,"SYMBOL"]
        if (is.null(plot.mpm.args$labels))
          stop("There is no variable named 'SYMBOL' in the pData of the object.\n
                  You may want to set the argument 'probe2gene' to FALSE (the default is TRUE)")
        
      }
      
      mpmPlot <- do.call(plot.mpm, plot.mpm.args)
      
      # add legend
		if (addLegend){
			colorsLegend <- plot.mpm.args$colors[-c(1, 2)]
			if(length(colorsLegend) > 0){
				par(font = 2)
				legend(legendPos, bty = "n", 
					legend = levels(pData(object)[, groups]),
					text.col = colorsLegend,
					cex = 1)
				par(font = 1)
			}
		}
      invisible(mpmPlot)
    })
