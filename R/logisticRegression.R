#' Logistic regression for predicting the probability to belong to a certain class
#'  in binary classification problems
#' 
#' Logistic regression for predicting the probability to belong to a certain class
#'  in binary classification problems.
#' @param object ExpressionSet object for the experiment
#' @param groups String containing the name of the grouping variable. This should be a 
#'  the name of a column in the \code{pData} of the \code{expressionSet} object.
#' @param probesetId The probeset ID. These should be stored in the \code{featureNames}
#'  of the \code{expressionSet} object.
#' @param geneSymbol The gene symbol. These should be stored in the column \code{`Gene Symbol`}
#'  in the \code{featureData} of the \code{expressionSet} object.
#' @param main Main title on top of the gra
#' @param probe2gene Boolean indicating whether the probeset should be translated to a gene symbol
#'  (used for the default title of the plot)
#' @param ... Possibility to add extra plot options. See \code{\link[base]{plot}}
#' @details 
#' It will always estimate probability scores to belong to the second level
#' of the factor variable. If a probability score to other level is preferred,
#' then you need to change the order of the levels of the factor.
#' @return 
#' A data.frame object with three columns and rownames
#'   \item{rownames }{The 'sampleNames' of the expressionSet object}
#'   \item{x }{The expression values for the specified gene for all samples}
#'   \item{y }{The labels of the samples}
#'   \item{fit }{The fitted probability score to belong to one of the two classes.}
#' @seealso \code{\link[a4Classif]{ROCcurve}},\code{\link{probabilitiesPlot}}
#' @author Willem Talloen
#' @examples 
#' \dontrun{
#' if (require(ALL)){
#'   data(ALL, package = "ALL")
#'   ALL <- addGeneInfo(ALL)
#'   ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'   logRegRes <- logReg(geneSymbol = "HLA-DPB1", object = ALL, groups = "BTtype")
#'   # scoresplot
#'   probabilitiesPlot(proportions = logRegRes$fit, classVar = logRegRes$y,
#'       sampleNames = rownames(logRegRes), main = 'Probability of being a T-cell type ALL')
#'   # barplot
#'   probabilitiesPlot(proportions = logRegRes$fit, classVar = logRegRes$y, barPlot=TRUE,
#'       sampleNames = rownames(logRegRes), main = 'Probability of being a T-cell type ALL')
#'  }
#' }
#' @importFrom Biobase exprs featureData featureNames pData
#' @importFrom stats glm fitted
#' @importFrom graphics par lines axis
#' @export
logReg <- function(object, groups, probesetId = NULL, 
		geneSymbol = NULL, main = NULL, probe2gene = TRUE, ...){
	
	dotCol <- 'goldenrod'
	lineCol <- 'blue'
	
	if (!is.null(probesetId) & !is.null(geneSymbol))
		stop("Please provide either a 'probeset' or a 'gene'")
	
	### create gene expression vector
	if (is.null(geneSymbol)){ # probeset given
		probesetId <- as.character(probesetId)     # use names not position !!
		exprGene <- exprs(object)[probesetId, ]
	} else { # gene given
		probesetPos <- which(geneSymbol == featureData(object)$`SYMBOL`)
		if (!length(probesetPos))
			stop("gene 'gene' does not occur in ExpressionSet 'object'")
		
		probesetId <- featureNames(object)[probesetPos]
		if (length(probesetId) > 1)
			warning(paste("Gene", geneSymbol, "corresponds to", length(probesetId), 
							"probesets; only the first probeset (",probesetId[1],") has been displayed on the plot."))
		exprGene <- exprs(object)[probesetId[1], ] # use names not position !!
	}
	
	labels <- as.factor(as.character(pData(object)[,groups]))
	
	# prepare title
	if (probe2gene){
		gSymbol <- featureData(object)[probesetId[1],]$`SYMBOL`
	}
	
	mainTitle <- if (is.null(main)){
				if (probe2gene)
					paste(gSymbol, " (", probesetId[1], ")", sep = "")
				else probesetId[1]
			} else { 
				main 
			}
	
	logResGene <- glm(labels ~ exprGene, family = 'binomial')
	logRes <- data.frame(x = exprGene, y = labels, fit = fitted(logResGene))
	logRes$y <- labels[]
	logRes <- logRes[order(logRes$x),]
	
	par(mar = c(5, 4, 4, 7) + 0.1)
	plot(x = logRes$x, y = as.numeric(logRes$y)-1, axes=FALSE,
			pch = 21, bg = dotCol, 
			xlab=expression(log[2] ~ intensity), 
			ylab= paste('Probability of being',levels(labels)[2]),
			main = mainTitle, ...)
	lines(logRes$x, logRes$fit, col = lineCol, lwd=2)
	axis(1, las=1); axis(2, las=1); box(bty='l') 
	axis(4, at= c(0,1), labels= levels(labels), 
			col.axis='black', tick = FALSE, las=1)
	par(mar = c(5, 4, 4, 2) + 0.1)
	
	return(logRes)
	#-----------
}

#' Function to plot the probabilities to belong to a certain class
#'  in binary classification problems.
#' 
#'  Function to plot the probabilities to belong to a certain class
#'   in binary classification problems. These probabilities are often calculated 
#'   using a logistic regression model. The class membership of the samples is
#'   displayed using a colored strip (with legend below the plot).
#' @param proportions A vector containing the calculated probabilities to belong to a certain class
#'  in binary classification problems. These probabilities are often calculated 
#'  using a logistic regression model.
#' @param classVar A vector containing the class where the sample belongs t
#' @param sampleNames A vector with the names of the samp
#' @param plot logical.  If \code{FALSE}, nothing is plotted
#' @param barPlot Should a barplot be drawn (\code{TRUE}) or a scatterplot like
#'   MCREstimate-type scores plot (the default, \code{FALSE}
#' @param layout boolean indicating whether \code{mcrPlot} should prespecify
#' a layout for a single plot (default, \code{TRUE}) or whether the user
#' takes care of the layout (\code{FALSE})
#' @param main Main title for the scores plot; if not supplied, 'Scores Plot'
#' is used as a defaul
#' @param sub Subtitle for the scores plot; if not supplied, the classification
#' technique and the chosen number of features are displayed
#' @param ... Additional graphical parameters to pass to the plot functio
#' @return no returned value, a plot is drawn in the current device.
#' @examples 
#' \dontrun{
#'   if (require(ALL)){
#'   data(ALL, package = "ALL")
#'     ALL <- addGeneInfo(ALL)
#'     ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'     logRegRes <- logReg(geneSymbol = "HLA-DPB1", object = ALL, groups = "BTtype")
#'     # scoresplot
#'     probabilitiesPlot(proportions = logRegRes$fit, classVar = logRegRes$y,
#'       sampleNames = rownames(logRegRes), main = 'Probability of being a T-cell type ALL')
#'     # barplot
#'     probabilitiesPlot(proportions = logRegRes$fit, classVar = logRegRes$y, barPlot=TRUE,
#'       sampleNames = rownames(logRegRes), main = 'Probability of being a T-cell type ALL')
#'   }
#' }
#' @seealso \code{\link{logReg}}
#' @author Willem Talloen and Tobias Verbeke
#' @importFrom graphics par layout plot axis abline points title rect barplot
#' @importFrom grDevices rgb
#' @export
probabilitiesPlot <- function(proportions,
		classVar,
		sampleNames,
		plot = TRUE,
		barPlot = FALSE, # barplot or MCREstimate-type scores plot
		layout = TRUE, 
		main = NULL,
		sub = NULL,
		...){ # additional arguments such as main, sub etc.
	
	def.par <- par(no.readonly = TRUE) # save default, for resetting...
	
	plotData <- proportions # vector of proportions
	names(plotData) <- sampleNames
	classVar <- as.factor(classVar) # after using names
	classVarLevels <- levels(classVar)
	
	##  data values are grouped according to their observed labels
	##  (e.g. responders / non responders)
	plotData <- plotData[order(classVar)] 
	classVar <- classVar[order(classVar)] # order the labels as well
	
	if (plot){
		if (!barPlot) {
			### layout
			if (layout) layout(matrix(1:2, ncol = 1), heights = c(6, 1))
			
			### upper plot
			plot(x = seq(along = plotData), # indices (names) of the samples 
					y = plotData,   # vector of proportions of misclassification for each sample
					ylim = c(0, 1),
					type = "n", las =  0, ann = FALSE, axes = FALSE, ...)
			
			axis(1, las = 3, at = seq(along = plotData), 
					labels = abbreviate(names(plotData)), cex.axis = 0.5)
      
			axis(2, las = 2)
			
			# draw grid first...
			abline(h = 0.5, col = "grey")
			abline(v = seq(length(plotData)), lty = "dashed", col = "grey")
			
			# ... then add dots  
			points(x = seq(along = plotData), y = plotData,
					pch = ifelse(plotData >= 0.5, 19, 17), cex = 1,
					col = ifelse(plotData >= 0.5, "darkblue", "orange"))
			
			title(main = if(is.null(main)){ 
								"Proportions Plot"
							} else {
								main
							})
			
			### add observed class membership
			nClasses <- length(classVarLevels)
			classColors <- switch(nClasses+1,   # based on RColorBrewer code            
          rgb(c(247,173,49),
              c(252,221,163),
              c(185,142,84),maxColorValue=255),
          rgb(c(255,194,120,35),
              c(255,230,198,132),
              c(204,153,121,67),maxColorValue=255),
          rgb(c(255,194,120,49,0),
              c(255,230,198,163,104),
              c(204,153,121,84,55),maxColorValue=255),
          rgb(c(255,217,173,120,49,0),
              c(255,240,221,198,163,104),
              c(204,163,142,121,84,55),maxColorValue=255),
          rgb(c(255,217,173,120,65,35,0),
              c(255,240,221,198,171,132,90),
              c(204,163,142,121,93,67,50),maxColorValue=255),
          rgb(c(255,247,217,173,120,65,35,0),
              c(255,252,240,221,198,171,132,90),
              c(229,185,163,142,121,93,67,50),maxColorValue=255),
          rgb(c(255,247,217,173,120,65,35,0,0),
              c(255,252,240,221,198,171,132,104,69),
              c(229,185,163,142,121,93,67,55,41),maxColorValue=255)
      )
      classColors <- classColors[2:(nClasses+1)]
          # brewer.pal(nClasses+3,"YlGn")[2:(nClasses+1)] # from MCREstimate
			
			sampleLocations <- seq(along = classVar)
			rect(xleft = sampleLocations-0.5, ybottom = -0.5, 
					xright = sampleLocations + 0.5, ytop = -0.015, 
					col = classColors[as.numeric(classVar)], 
					border = classColors[as.numeric(classVar)])
			abline(v = sampleLocations, lty = 2, col = "grey")
			
			### lower plot (with legends)
#			op <- par(mar = c(0,4,0,2))
#			plot(c(0, 1), type = "n", ann = FALSE, axes = FALSE)
#			legend("left", legend = classVarLevels, fill = classColors,
#					bty = "n")
#			legend("right", legend = c("0.5 <= score <=   1", 
#							"   0 <= score <  0.5"), 
#					pch = c(19, 17), pt.cex = 1.5, col = c("darkblue", "orange"), 
#					bty = "n")
#			par(op)
#						
		} else {
			barplot(height = plotData, col = ifelse(plotData >= 0.5, "green", "red"),  
					las = 3, ...)
			abline(h = 0.5, col = "grey")    
		}
	}
	
	par(def.par)#- reset to default
	invisible(plotData) # named vector of proportion correctly classified
}

