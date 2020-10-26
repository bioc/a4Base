#' Pick One or More OA Colors
#' @param color a character vector of color names; possible values are "red", "orange", "yellow", "green", "cyan", "blue", "pink", 
#'       "limegreen", "purple", "black", "white", "grey" or "gray" 
#' @param alpha transparency level for the color(s)
#' @return character vector of colors 
#' @author Tobias Verbeke
#' @importFrom grDevices hcl
oaColors <- function(color = NULL, alpha = 1.0){
	
	colPalette1 <- hcl(seq(0, 360, by = 10 ), l = 50, c = 250, alpha = alpha)
	colPalette2 <- hcl(seq(250, 270, by = 1), l = 50, c = 360, alpha = alpha)
	
	colorValues <- c(colPalette1[2], colPalette2[6], colPalette2[15], colPalette1[12], 
			colPalette1[21], colPalette1[27], colPalette1[29], 
			rainbow(60, alpha = alpha)[17], rainbow(360, alpha = alpha)[278], 
			hcl(180, c = 0, l = 0, alpha = alpha), hcl(269, l = 100, c = 90, alpha = alpha), 
			hcl(180, l = 60, c = 0, alpha = alpha), hcl(180, l = 60, c = 0, alpha = alpha)
	)
	names(colorValues) <- c("red", "orange", "yellow", "green", "cyan", "blue", "pink", 
			"limegreen", "purple", "black", "white", "grey", "gray")
	
	out <- colorValues[color]
	
	if(any(is.na(out)))
		stop(
				cat("Error: ", color[which(is.na(out))], " is not in the list of oaColors. Current colors include: ", "\n         ", 
						paste(paste(names(colorValues), collapse = ", "), sep = ""), "\n", sep = "")
		)
	
	
	return(out)
	
}



#' Generate a Palette of OA Colors
#' @param numColors number of colors to be contained in the palette
#' @param alpha transparency level of the colors
#' @return vector of colors
#' @author Jason Waddell
oaPalette <- function(numColors = NULL, alpha = 1.0){
	
	# check numColors = NULL, test range
	
	fullPaletteNames <- c("red", "blue", "green", "orange", "pink", "cyan", "yellow", "limegreen", "purple")
	
	samplePaletteNames <- fullPaletteNames[1:numColors]
	samplePalette <- oaColors(samplePaletteNames, alpha = alpha)
	return(samplePalette)
	
}

#' Utility function that defines a color palette for use in a4
#' @param n Number of color levels the palette should provid
#' @param alpha alpha transparency level of the colors
#' @param Janssen logical.  If \code{TRUE}, Janssen Pharmaceutical colors
#' are used (with a maximum of 6 possible colors).
#' @details
#'  For n = 1, \code{"blue"} is returned; for n = 2 
#'  \code{c("red", "blue")} is returned; for n = 3 
#'  \code{c("red", "green", "blue"} is returned; for n = 4 
#'  \code{c("red", "green", "blue", "purple")} is returned and for n > 2, 
#'  the output of \code{rainbow(n)} is returned.
#' @return a character vector of colors
#' 
#' @author Steven Osselaer, Tobias Verbeke
#' @seealso \code{rainbow} palette in \code{\link[grDevices]{palettes}}
#' @examples 
#'  op <- par(mfrow = c(2, 3))
#'  for (nGroups in 1:6)
#'   pie(rep(1, nGroups), a4palette(nGroups))
#'  par(op)
#' @keywords dplot
#' @importFrom grDevices rgb rainbow
#' @export
a4palette <- function(n, alpha = 1.0, Janssen = FALSE){
	if (!is.numeric(n) | n < 1)
		stop("'n' should be a positive integer")
	
	if(Janssen){
		# Janssen ppt template colors
		col1 <- rgb(9,53,122, alpha = alpha, maxColorValue = 255)
		col2 <- rgb(42,142,191, alpha = alpha, maxColorValue = 255)
		col3 <- rgb(148,198,223, alpha = alpha, maxColorValue = 255)
		col4 <- rgb(35,125,38, alpha = alpha, maxColorValue = 255)
		col5 <- rgb(127,195,28, alpha = alpha, maxColorValue = 255)
		col6 <- rgb(191,225,141, alpha = alpha, maxColorValue = 255)
		
		res <- if (n==1) col2
				else if (n==2) c(col2, col4)
				else if (n==3) c(col1, col2, col4)
				else if (n==4) c(col1, col2, col4, col5)
				else rainbow(n)
		return(res)
		
	} else {if (n < 10) oaPalette(n, alpha)	else if (n==2) rainbow(n)}
}

#' Create a Profile Plot for a given Gene
#' 
#'  Create a profile plot for a given gene. A profile plot displays the expression values (y-axis)
#'    by samples (x-axis), sorted by group. This is a useful working graph as samples can be
#'    directly identified. For presentation purposes, a \code{boxPlot} can also be considered. with jittered for readability of the plot.
#' @param probesetId The probeset ID. These should be stored in the \code{featureNames}
#'  of the \code{expressionSet} object.
#' @param geneSymbol The gene symbol. These should be stored in the column \code{`Gene Symbol`}
#'   in the \code{featureData} of the \code{expressionSet} object
#' @param object ExpressionSet object for the experiment
#' @param groups String containing the name of the grouping variable. This should be a 
#'  name of a column in the \code{pData} of the \code{expressionSet} object.
#' @param main Main title on top of the graph
#' @param colvec Vector of colors to be used for the groups. If not specified, the default colors of
#'   \code{a4palette} are used.
#' @param colgroups String containing the name of the variable to color the superimposed dots.
#'   This should be a the name of a column in the \code{pData} of the \code{expressionSet} object
#' @param probe2gene Boolean indicating whether the probeset should be translated to a gene symbol
#'    (used for the default title of the plot)
#' @param sampleIDs A boolean or a string to determine the labels on the x-axis. Setting it to FALSE
#'  results in no labels (interesting when the labels are unreadable due to large sample sizes).
#'  Setting it to a string will put the values of that particular \code{pData} column as labels.
#'  The string should be a name of a column in the \code{pData} of the \code{expressionSet} object."
#' @param addLegend Boolean indicating whether a legend for the colors of the dots should be added.
#' @param legendPos Specify where the legend should be placed. Typically either \code{topright},
#'   \code{bottomright}, \code{topleft} (the default) or \code{bottomleft}
#' @param cex character expansion used for the plot symbols; defaults to 1.5
#' @param ... Further arguments, e.g. to add extra plot options. See \code{\link{par}}
#' @return   If a \code{geneSymbol} is given that has more than one probeSet,
#' the plots for only the first probeSet is displayed.
#'   A character vector of corresponding probeset IDs is returned invisibly,
#'   so that one can check the profiles of the other related probeset IDs with
#'   an extra \code{plot1gene} statement 
#'  If a \code{probesetId} is given, one single profile plot for the probeset is 
#'   displayed.
#' @seealso \code{\link{plotCombination2genes}}, \code{\link{boxPlot}}
#' @author S. Osselaer, W. Talloen, T. Verbeke
#' @examples 
#' if (require(ALL)){
#'   data(ALL, package = "ALL")
#'   ALL <- addGeneInfo(ALL)
#'   # one variable (specified by groups)
#'   plot1gene(geneSymbol = 'HLA-DPB1', object = ALL, groups = "BT",
#' 	    addLegend = TRUE, legendPos = 'topright')
#'   # two variables (specified by groups and colGroups)
#'  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'   plot1gene(probeset = '1636_g_at', object = ALL, groups = 'BT',
#'       colgroups = 'mol.biol', legendPos='topright', sampleIDs = 'BT')    
#' }
#' @keywords dplot
#' @importFrom Biobase pData featureData featureNames exprs
#' @importFrom graphics points axis mtext lines legend
#' @export
plot1gene <- function (probesetId = NULL, 
    geneSymbol = NULL, 
    object, groups, main = NULL, colvec = NULL,
    colgroups = NULL, 
    probe2gene = TRUE, sampleIDs = TRUE, 
    addLegend = TRUE, legendPos = "topleft", cex = 1.5, ...) {
  
  if (!is.null(probesetId) & !is.null(geneSymbol))
    stop("Please provide either a 'probesetId' or a 'geneSymbol'")
  
  if((length(sampleIDs) > 1) | !(is.logical(sampleIDs) | is.character(sampleIDs)))
    stop("'sampleIDs' should either be a logical or a character of length one")
  
  dotList <- list(...)
  
  groups <- pData(object)[, groups]
  if (!is.factor(groups)){
    groups <- as.factor(groups)
    warning("The variable referenced by 'groups' was coerced to a factor for the purpose of this visualization")
  }
  
  groups <- groups[, drop=TRUE] # remove unused groups
  nByGroup <- as.numeric(table(groups))
  
  if (is.null(geneSymbol)){ # probeset given
    probesetId <- as.character(probesetId)     # use names not position !!
    plotData <- exprs(object)[probesetId, ]
  } else { # gene given
    probesetPos <- which(geneSymbol == featureData(object)$`SYMBOL`)
    if (!length(probesetPos))
      stop("gene 'gene' does not occur in ExpressionSet 'object'")
    
    probesetId <- featureNames(object)[probesetPos]
    if (length(probesetId) > 1)
      warning(paste("Gene", geneSymbol, "corresponds to", length(probesetId), 
              "probesets; only the first probeset (",probesetId[1],") has been displayed on the plot."))
    plotData <- exprs(object)[probesetId[1], ] # use names not position !!
  }
  
  # order data by factor level
  orderGroups <- order(groups)
  groups <- groups[orderGroups]
  numericGroups <- as.numeric(groups) # after ordering
  
  plotData <- plotData[orderGroups]
  nc <- length(plotData)
  
  # define colors
  if (!is.null(colgroups)) {
    colgroups <- pData(object)[, colgroups]
    if (!is.factor(colgroups))
      stop("The variable referenced by 'colgroups' should be a factor")
    colgroups <- colgroups[, drop=TRUE] # remove unused groups
    colgroups <- colgroups[orderGroups]
    numericColgroups <- as.numeric(as.factor(colgroups)) # after ordering
    colGroupsNotDefinedAsArgument <- FALSE
  } else {
    colgroups <- groups
    numericColgroups <- numericGroups
    colGroupsNotDefinedAsArgument <- TRUE
  }
  
  if (is.null(colvec)){
    colvec <- a4palette(nlevels(colgroups))	  
  } else {
    if(length(colvec) != nlevels(colgroups))
      stop("'colvec' should contain as many elements as there are levels in 'groups' or 'colgroups'")
  }
  
  # prepare title
  if (probe2gene){
    gSymbol <- featureData(object)[probesetId[1],]$`SYMBOL`
  }
  
  # cexMain <- if (is.null(dotList$cex.main)) par("cex.main") else dotList$cex.main
  # cexLab <- if (is.null(dotList$cex.lab)) par("cex.lab") else dotList$cex.lab
  
  mainTitle <- if (is.null(main)){
        if (probe2gene)
          paste(gSymbol, " (", probesetId[1], ")", sep = "")
        else probesetId[1]
      } else { 
        main 
      }
  
  argList <- list(x = 1:(nc + 1), y = c(mean(plotData), plotData), type = "n", axes = FALSE, xlab = "", ylab = expression(log[2] ~ intensity),
      main = mainTitle)
  argList <- c(argList, dotList)
  do.call(plot, args = argList)
  
  #  plot(1:(nc + 1), c(mean(plotData), plotData), type = "n",
  #      axes = FALSE, xlab = "", ylab = expression(log[2] ~ intensity),
  #      main = mainTitle, cex.main = cexMain, cex.lab = cexLab)
  
  points(2:(nc + 1), plotData, bg = colvec[numericColgroups], pch = 21, 
      cex = cex)
  
  axis(2, las = 2, cex.axis = 0.7, lwd = 1.5)
  
  if (is.logical(sampleIDs)){
    if (sampleIDs){
      axis(1, las = 3, at = c(2:(nc + 1)), labels = names(plotData),
          cex.axis = 0.7, lwd = 1.5, las = 2)  
    } else {
      emptyLabels <- rep("", length(plotData))
      axis(1, las = 3, at = c(2:(nc + 1)), labels = emptyLabels,
          cex.axis = 0.7, lwd = 1.5, las = 2)
    }
  } else {
    sampleIDs <- pData(object)[, sampleIDs]
    axis(1, las = 3, at = c(2:(nc + 1)), labels = sampleIDs[orderGroups],
        cex.axis = 0.7, lwd = 1.5, las = 2)
  }
  mtext("labels Means", 1, at = 1, font = 2, cex = 0.7, las = 3)
  
  # add lines
  j <- 0
  if (colGroupsNotDefinedAsArgument) {
    for (i in 1:max(numericGroups)) {
      med <- mean(as.vector(plotData[(j + 1):(j + nByGroup[i])], mode = "numeric"))
      
      lines(x = c(0, 1.5), y = c(med, med), col = colvec[i], lty = 1, lwd = 2)
      lines(x = c(j + 2, j + 1 + nByGroup[i]), y = c(med, med), col = colvec[i], lty = 1, lwd = 2)
      
      j <- j + nByGroup[i] 
    }
  } else {			# if colored by an extra factor, the median lines should be in grey 
    for (i in 1:max(numericGroups)) {
      med <- mean(as.vector(plotData[(j + 1):(j + nByGroup[i])], mode = "numeric"))
      
      lines(x = c(0, 1.5), y = c(med, med), col = 'grey', lty = 1, lwd = 2)
      lines(x = c(j + 2, j + 1 + nByGroup[i]), y = c(med, med), col = 'grey', lty = 1, lwd = 2)
      
      j <- j + nByGroup[i] 
    } 
  } 
  if (addLegend){
    legend(legendPos, bty='n', 
        legend = levels(colgroups),
        text.col = colvec, cex=1)
  }
  invisible(probesetId)
}

# plot boxplot per gene with raw data superimposed
#' Create a boxplot for a given gene.
#' 
#'  Create a boxplot for a given gene. The boxplot displays the expression values (y-axis)
#'    by groupss (x-axis). The raw data are superimposed as dots, jittered for readability of the plot.
#'    Optionally, the dots can be colored by another variable.
#' @param probesetId The probeset ID. These should be stored in the \code{featureNames}
#'   of the \code{expressionSet} object.
#' @param geneSymbol The gene symbol. These should be stored in the column \code{`Gene Symbol`}
#'    in the \code{featureData} of the \code{expressionSet} object.
#' @param object ExpressionSet object for the experiment
#' @param groups String containing the name of the grouping variable. This should be a 
#'   the name of a column in the \code{pData} of the \code{expressionSet} object.
#' @param main Main title on top of the graph
#' @param colvec Vector of colors to be used for the groups. If not specified, the default colors of
#'    \code{a4palette} are used.
#' @param colgroups String containing the name of the variable to color the superimposed dots.
#'    This should be a the name of a column in the \code{pData} of the \code{expressionSet} object
#' @param probe2gene Boolean indicating whether the probeset should be translated to a gene symbol
#'     (used for the default title of the plot)
#' @param addLegend Boolean indicating whether a legend for the colors of the dots should be added.
#' @param legendPos Specify where the legend should be placed. Typically either \code{topright},
#'    \code{bottomright}, \code{topleft} (the default) or \code{bottomleft}
#' @param ... Possibility to add extra plot options. See \code{\link{par}}
#' @return A plot is drawn to the current device and
#' \code{probesetId} are returned invisibly.
#' @seealso \code{\link{plot1gene}}
#' @examples 
#' # simulated data set
#' esSim <- simulateData()
#' boxPlot(probesetId = 'Gene.1', object = esSim, groups = 'type', addLegend = FALSE)
#' # ALL data set
#' if (require(ALL)){
#'   data(ALL, package = "ALL")
#'   ALL <- addGeneInfo(ALL)
#'   ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#'   boxPlot(geneSymbol = 'HLA-DPB1', object = ALL, boxwex = 0.3,
#' 		  groups = 'BTtype', colgroups = 'BT', legendPos='topright')
#' }
#' @author Willem Talloen
#' @importFrom Biobase pData exprs featureData
#' @importFrom graphics boxplot points legend
#' @export
boxPlot <- function(probesetId = NULL, 
    geneSymbol = NULL, 
    object, groups, main = NULL, colvec = NULL,
    colgroups = NULL, probe2gene = TRUE, 
    addLegend = TRUE, legendPos = "topleft", ...) {
  
  if (!is.null(probesetId) & !is.null(geneSymbol))
    stop("Please provide either a 'probeset' or a 'gene'")
  
  groups <- pData(object)[, groups]
  if (!is.factor(groups))
    stop("The variable referenced by 'groups' should be a factor")
  groups <- groups[, drop=TRUE] # remove unused groups
  numericGroups <- as.numeric(groups) # after ordering
  
  if (is.null(geneSymbol)){ # probeset given
    probesetId <- as.character(probesetId)     # use names not position !!
    plotData <- exprs(object)[probesetId, ]
  } else { # gene given
    probesetPos <- which(geneSymbol == featureData(object)$`SYMBOL`)
    if (!length(probesetPos))
      stop("gene 'gene' does not occur in ExpressionSet 'object'")
    
    probesetId <- featureNames(object)[probesetPos] 
    if (length(probesetId) > 1)
      warning(paste("Gene", geneSymbol, "corresponds to", length(probesetId), 
              "probesets; only the first probeset (",probesetId[1],") has been displayed on the plot."))
    plotData <- exprs(object)[probesetId[1], ] # use names not position !!
  }
  
  # define colors
  if (!is.null(colgroups)) {
    colgroups <- pData(object)[, colgroups]
    if (!is.factor(colgroups))
      stop("The variable referenced by 'colgroups' should be a factor")
    colgroups <- colgroups[, drop=TRUE] # remove unused groups
    numericColgroups <- as.numeric(as.factor(colgroups)) # after ordering
    colGroupsNotDefinedAsArgument <- FALSE
  } else {
    colgroups <- groups
    numericColgroups <- numericGroups
    colGroupsNotDefinedAsArgument <- TRUE
  }
  
  if (is.null(colvec)){
    colvec <- a4palette(nlevels(colgroups))	  
  } else {
    if(length(colvec) != nlevels(colgroups))
      stop("'colvec' should contain as many elements as there are levels in 'groups' or 'colgroups'")
  }
  
  # prepare title
  if (probe2gene){
    gSymbol <- featureData(object)[probesetId[1],]$`SYMBOL`
  }
  
  mainTitle <- if (is.null(main)){
        if (probe2gene)
          paste(gSymbol, " (", probesetId[1], ")", sep = "")
        else probesetId
      } else { 
        main 
      }
  
  # make plot
  boxplot(plotData~groups, outline=FALSE, col='grey', type='n',
      ylab=expression(log[2]~concentration), las=1, main = mainTitle, ...)
  points(jitter(as.numeric(groups)),plotData,
      bg = colvec[numericColgroups], pch = 21, cex=1.5)
  
  # add legend
  if (addLegend){
    legend(legendPos, bty='n', 
        legend = levels(colgroups),
        text.col = colvec, cex=1)
  }
  invisible(probesetId)
}

# plot combination of two genes

#' Plot a Combination of Two Genes
#' 
#' Plot a Combination of Two Genes
#' @param probesetId1 First probeset id, plotted in the x-axis
#' @param probesetId2 Second probeset id, plotted in the y-axis
#' @param geneSymbol1 First gene symbol, plotted in the x-axis
#' @param geneSymbol2 Second gene symbol, plotted in the y-axi
#' @param object ExpressionSet object for the experiment
#' @param groups string containing the name of the grouping variable
#' @param addLegend Logical value to indicate whether a legend needs to be draw
#' @param legendPos Position on the graph where to put the legend
#' @param probe2gene should the probeset be translated to a gene symbol
#'  (used for the default title of the plot
#' @param colvec a character vector of colors. If not specified it will be 
#'  automatically generated by \code{a4palette}
#' @param ... This allows to specify typical arguments in the \code{plot} functio
#' @return If a gene id is given, the plots for only the first probeset is displayed 
#'   and a character vector of corresponding probeset IDs is returned invisibly.  
#'   It is a list containing
#'   \item{probeset1 }{Probeset ids measuring 'gene1'}
#'   \item{probeset1 }{Probeset ids measuring 'gene1'}
#'   If a probeset id is given, one single profile plot for the probeset is 
#'  displayed.
#' @seealso \code{\link{plot1gene}}
#' @author W. Talloen, T. Verbeke
#' @examples 
#' if (require(ALL)){
#'  data(ALL, package = "ALL")
#'  ALL <- addGeneInfo(ALL)
#'  aa <- plotCombination2genes(geneSymbol1 = 'HLA-DPB1', geneSymbol2 = 'CD3D',
#'			object = ALL, groups = "BT",
#'			addLegend = TRUE, legendPos = 'topright')
#'  aa
#'}
#' @importFrom Biobase exprs featureNames pData
#' @importFrom graphics plot points legend
#' @export
plotCombination2genes <- function(probesetId1 = NULL, probesetId2 = NULL,
    geneSymbol1 = NULL, geneSymbol2 = NULL,
    object, groups, addLegend = TRUE, legendPos = "topleft",
    probe2gene = TRUE, colvec = NULL, ...) {
  
  if (!is.null(probesetId1) & !is.null(probesetId2) & !is.null(geneSymbol1) & !is.null(geneSymbol2))
    stop("Please provide either two probesets or two genes")
  
  # change gene into probeset
  if (is.null(geneSymbol1)){ # probeset given
    probesetId1 <- as.character(probesetId1)     # use names not position !!
  } else { # gene given
    probesetPos1 <- which(geneSymbol1 == featureData(object)$`SYMBOL`)
    if (!length(probesetPos1))
      stop("gene 'gene1' does not occur in ExpressionSet 'object'")
    
    probesetId1 <- featureNames(object)[probesetPos1]
    if (length(probesetId1) > 1)
      warning(paste("Gene1", geneSymbol1, "corresponds to", length(probesetId1), 
              "probesets; only the first probeset (",probesetId1[1],") has been displayed on the plot."))
  }
  
  if (is.null(geneSymbol2)){ # probeset given
    probesetId2 <- as.character(probesetId2)     # use names not position !!
  } else { # gene given
    probesetPos2 <- which(geneSymbol2 == featureData(object)$`SYMBOL`)
    if (!length(probesetPos2))
      stop("gene 'gene2' does not occur in ExpressionSet 'object'")
    
    probesetId2 <- featureNames(object)[probesetPos2]
    if (length(probesetId2) > 1)
      warning(paste("Gene2", geneSymbol2, "corresponds to", length(probesetId2), 
              "probesets; only the first probeset (",probesetId2[1],") has been displayed on the plot."))
  }
  groups <- factor(pData(object)[, groups])[, drop=TRUE]
  exprGene1 <- exprs(object)[probesetId1[1], ]
  exprGene2 <- exprs(object)[probesetId2[1], ] 
  
  if (probe2gene){
    gSymbol1 <- featureData(object)[probesetId1[1],]$`SYMBOL`
    gSymbol2 <- featureData(object)[probesetId2[1],]$`SYMBOL`
  }
  
  if (is.null(colvec))
    colvec <- a4palette(nlevels(groups))
  
  # make plot
  plot(exprGene1, exprGene2, type = "n", 
      xlab = if (probe2gene) gSymbol1 else probesetId1[1], 
      ylab = if (probe2gene) gSymbol2 else probesetId2[1], 
      las = 1, ...)
  points(exprGene1, exprGene2, pch=21, cex=1.5,
      bg = colvec[as.numeric(groups)], ...)
  
  # add legend
  if (addLegend){
    legend(legendPos, bty='n', 
        legend = levels(groups),
        text.col = colvec, cex=1)
  }
  invisible(list(probeset1=probesetId1, probeset2=probesetId2))
}

# parallel coordinate plots
#' Plot expression profiles of multiple genes or probesets
#' Plot expression profiles of multiple genes or probesets. Each line depicts a gene,
#' and the color legend can be used to identify the gene.
#' @param object ExpressionSet object for the experiment
#' @param probesetIds The probeset ID. These should be stored in the \code{featureNames}
#'   of the \code{expressionSet} object.
#' @param sampleIDs A boolean or a string to determine the labels on the x-axis. Setting it to FALSE
#'  results in no labels (interesting when the labels are unreadable due to large sample sizes).
#'  Setting it to a string will put the values of that particular \code{pData} column as labels.
#'  The string should be a name of a column in the \code{pData} of the \code{expressionSet} object
#' @param addLegend Boolean indicating whether a legend for the colors of the dots should be added.
#' @param legendPos Specify where the legend should be placed. Typically either \code{topright},
#'   \code{bottomright}, \code{topleft} (the default) or \code{bottomleft}
#' @param colvec Vector of colors to be used for the groups. If not specified, the default colors of
#'   \code{a4palette} are used
#' @param orderGroups String containing the name of the grouping variable to order the samples 
#'  in the x-axis accordingly. This should be a name of a column in the \code{pData} of
#'   the \code{expressionSet} object.
#' @param ... Possibility to add extra plot options. See \code{\link{par}}
#' @examples 
#' if (require(ALL)){
#' data(ALL, package = "ALL")
#' ALL <- addGeneInfo(ALL)
#' ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
#' myGeneSymbol <- c("LCK")	# a gene 
#' probesetPos <- which(myGeneSymbol == featureData(ALL)$SYMBOL)
#' myProbesetIds <- featureNames(ALL)[probesetPos]
#' profilesPlot(object = ALL, probesetIds = myProbesetIds, 
#'       orderGroups = "BT", sampleIDs = "BT")
#' }
#' @return No returned value, a plot is drawn in the current device.
#' @seealso \code{\link{plot1gene}}, \code{\link{boxPlot}}
#' @author W. Talloen
#' @importFrom graphics matplot
#' @export
profilesPlot <- function (object, probesetIds, sampleIDs = TRUE, 
    addLegend = TRUE, legendPos = "topleft", colvec = NULL,
    orderGroups = NULL, ...) {
  
  if (length(probesetIds) < 2)
    stop("Please provide at least two 'probesetIds'")
  
  if((length(sampleIDs) > 1) | !(is.logical(sampleIDs) | is.character(sampleIDs)))
    stop("'sampleIDs' should either be a logical or a character of length one")
  
  plotData <- t(exprs(object)[probesetIds,])
  
  # order data by factor level
  if (!is.null(orderGroups)){
    orderGroups <- as.factor(as.character(pData(object)[, orderGroups]))
    orderGroups2 <- order(orderGroups)
    plotData <- plotData[orderGroups2,]
  }
  
  if (is.null(colvec)){
    colvec <- a4palette(ncol(plotData))	  
  } else {
    if(length(colvec) != ncol(plotData))
      stop("'colvec' should contain as many elements as there are levels in 'groups' or 'colgroups'")
  }
  
  # plot
  matplot(plotData, type = "l",
      xlab = "", ylab = expression(log[2]~concentration),
      axes = FALSE, lwd = 1, lty = 1, col = colvec, ...)
  axis(2, las = 2, cex.axis = 0.7, lwd = 1.5)
  
  # x-axis
  if (is.logical(sampleIDs)){
    if (sampleIDs){
      axis(1, labels = rownames(plotData), las = 3, at = 1:nrow(plotData),
          cex.axis = 0.7, lwd = 1.5)
    } else {
      emptyLabels <- rep("", length(plotData))
      axis(1, labels = emptyLabels, las = 3, at = 1:nrow(plotData),
          cex.axis = 0.7, lwd = 1.5)
    }
  } else {
    sampleIDs <- pData(object)[, sampleIDs]
    axis(1, labels = sampleIDs[orderGroups2], las = 3, at = 1:nrow(plotData),
        cex.axis = 0.7, lwd = 1.5)
  }
  
  if (addLegend){
    legend(legendPos, bty = "n", 
        legend = colnames(plotData),
        text.col = colvec, cex = 1)
  }
}
