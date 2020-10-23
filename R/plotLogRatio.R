#' Plot a summary gene expression graph
#' 
#' Plot ratios of expression values observed in a treatment versus those of a reference. 
#' First the ratios and variances are computated on the gene expression data.
#' @section Ordering:
#' orderBy: A list with two components, rows and cols, each one possibly being NULL (no ordering 
#' on the specific dimension). Ordering on cols can be done according to (a) pData column(s) 
#' (for example: \code{c('cellline','compound','dose'}. Ordering on rows can be done using of the 
#' following values:
#'  \itemize{
#'   \item{NULL}{no reordering on rows}
#'   \item{numeric vector}{use the vector values to sort rows}
#'   \item{alpha}{use genes names alphabetice order}
#'   \item{effect}{try to assess global gene expression level by taking sum(abs(values)) on specified exprs columns)}
#'  \item{hclust}{use the ordering returned by \code{hclust} invoked on specified exprs columns}
#'  }
#' @section Colors:
#'  The management of colors is very flexible but is a little bit tricky, as a variety of parameters 
#'  are available to the user. Basically, combinations of arguments allow to set colors for columns headers (text), 
#'  columns as a whole (different colors for the different columns) or for each of the inividual horizontal bars.
#'  By default, everything is red. There are four main different arguments that can be used and that are 
#'  applied in a consecutive order. Each one may override a previous argument value. Below is a list of 
#'  arguments and their consecutive actions:
#' \itemize{
#'  \item{\code{colorsColumns}}{ The first way to assign colors is to provide a vector of colors that will 
#'    be used for each column (headers and its horizontal bars). This vector is recycled so that providing one unique 
#'    value will color all columns, whereas providing a vector of length 2 will alternate columnns colors.}
#'  \item{\code{colorsColumnsBy}}{To be used when the experiment involves groupings for pData, for example dose, 
#'   cellline or treatment. In order to see the effects of such variables, one can color columns using 
#'   combinations of those. The argument is a vector of pData columns such as \code{c('cellline','dose')}. 
#'    Unique combinations will be computed and a color will be assigned for each group of columns. 
#'    The vector that is provided with the argument \code{colorsColumnsByPalette} is used to assign colors. 
#'    If the argument \code{colorColumnsBy} is not \code{NULL} then it overrides the previous argument \code{colorsColumns}.}
#'  \item{\code{colorsUseMeanQuantiles}}{ A logical value. The default plotGeneDE displays for each gene the expression value difference 
#'    between treatment and reference, but does not reveal any information about the expression levels in these conditions. 
#'    Parameter \code{colorsUseMeanQuantiles} allows to color the horizontal bars according to expression level that 
#'    is derived from quantiles computed on averages of the complete ExpressionSet object. 
#'    As it involves the expression data of all probesets, computations must be done 
#'    before subsetting the ExpressionSet object and the plotGeneDEting. The function \code{\link{addQuantilesColors}} 
#'    computes quantiles and corresponding mean expression level intervals. If \code{colorsUseMeanQuantiles} 'TRUE', 
#'    previous coloring parameters are overriden. The parameter \code{colorsMeanQuantilesPalette} is used to assign 
#'    colors for average-quantiles-groups. Note that columns headers are still given by previous arguments.}
#'  \item{\code{colorsBarsMatrix}}{The most flexible way to assign colors as the matrix will be used to color each bar 
#'    of the plot individually. A check is done to ensure that the number of rows and columns are not less than the number of 
#'    probesets and columns. If not \code{NULL}, this parameter overrides the previous ones.}
#' }
#' @param e ExpressionSet object to use
#' @param filename Name of the filename to use. No need to specify extension which 
#'    will be added according to device.
#' @param device One of 'pdf', 'X11', 'png', 'svg'. For svg device, one X11 device 
#'    is also opened.
#' @param orderBy See details
#' @param colorsColumns A vector of colors to be used for plotting columns; default value 
#'   is NULL which ends up with red -- see Colors section
#' @param colorsColumnsBy A vector of pData columns which combinations 
#' specify different colors to be used -- see Colors section
#' @param colorsColumnsByPalette If colorsColumns is NULL, 
#' vector of colors to be used for coloring columns potentially 
#'  splitted by colorsColumnsBy
#' @param colorsUseMeanQuantiles Boolean to indicate if the quantile groups computed 
#'    on averages over all treatments should be used for coloring -- see Colors section
#' @param colorsMeanQuantilesPalette if colorsUseMeanQuantiles is TRUE, these colors 
#'    will be used for the different groups -- see Colors section
#' @param colorsBarsMatrix Matrix of colors to be used for each individual bar; 
#'    colors are provided for genes in data order and thus are possibly reordered 
#'    according to orderBy -- see Colors sectio
#' @param colorsGenesNames Vector of colors to be used for gene names; will be recycled 
#'    if necessary; colors are provided for genes in data order and thus are possibly 
#'    reordered according to orderBy
#' @param main Main title
#' @param shortvarnames ector or pData column to be used to 
#' display in graph columns. If NULL, those names 
#'    will be used from the coded names added 
#'	to pData during computations (list of columns values pasted with a dot).
#'	Warning: shortvarnames must be defined in the order columns are 
#'	present in the ExpressionSet object so that they will be 
#'	reordered if one asks to order columns.
#' @param longvarnames pData column to be used in SVG tooltip title.
#'  If NULL, shortvarnames will be used.
#' Same warning than shortvarnames about ordering
#' @param gene.length Maximum number of characters that will be printed of the gene names
#' @param gene.fontsize Font size for the gene names , default = 
#' @param main.fontsize Font size for the main, default = 9
#' @param columnhead.fontsize Font size for the column headers, default = 8
#' @param mx Expansion factor for the 
#' width of the bars that represent the expression ratios
#' @param exp.width Expansion factor for global graph width, and the space between the plotted colum
#' @param exp.height Expansion factor for global graph height, and the space between the plotted row
#' @param log2l.show A logical value. If 'TRUE', the line for 
#' log2 values on each column (when max(data) > 2) is draw
#' @param log4l.show A logical value. If 'TRUE', the line for log4 
#' values on each column (when max(data) > 4) is drawn
#' @param quantiles.show A logical value. If 'TRUE', a line 
#' is drawn for quantiles computed separately on each column
#' @param quantiles.compute A logical value. If 'TRUE', the vector quantiles will be computed and displayed
#'    provided that \code{quantile.show} is \code{TRUE}
#' @param error.show A logical value. If 'TRUE', errors bars are displayed on the graph 
#'    (only for those columns for which they are available
#' @param view.psid A logical value. If 'TRUE', the genes psid is displayed on the gene name
#' @param errorLabel A character vector describing the error bars, printed at the bottom of the figu
#' @param closeX11 If \code{device} is SVG, do we close the required X11 device at the end?
#' @param openFile A logical value. If 'TRUE', the produced output file is opened
#' @param tooltipvalues If device is SVG, one can choose to display each bar separately, with data values as tooltips. 
#'    Note however that each bar will be considered as a distinct object instead of a column, which will takes much 
#'    more time to create the graph and produces a much bigger SVG file
#' @param probe2gene Boolean indicating whether the probeset should be translated to a gene symbol
#'    (used for the default title of the plot
#' @param ... \code{\dots}
#' @inheritParams computeLogRatio
#' @return The ExpressionSet object with the computated variables is returned.
#' @seealso \code{\link{computeLogRatio}},\code{\link{addQuantilesColors}}
#' @author Hinrich Goehlmann and Eric Lecoutre
#' @example inst/examples/plotLogRatio-example.R
#' @importFrom Biobase featureNames featureData pData fData exprs openPDF package.version
#' @importFrom stats dist hclust quantile
#' @importFrom grDevices rainbow pdf png dev.off x11 colors
#' @importFrom grid viewport pushViewport popViewport grid.rect grid.lines grid.segments grid.text gpar gPath
#' @importFrom utils browseURL
#' @importFrom a4Preproc addGeneInfo
#' @export
plotLogRatio <- function(
		e,
		reference,
		within=NULL,
		across=NULL,
		nReplicatesVar=3, ## see parameters of compute function
		filename = "Rplots", # extension svg or pdf will be added according to device -- please provide full path if you want to automatically open the file
		device="svg", # svg, pdf, X11 for SVG, requires X11 and so interactive session
		
		# if col=NULL then color according to quantiles of ##TODOchange to columnsColors#test to see if quantiles over whole objects are already computed and available in expressionSetObject
		orderBy=list(rows='hclust',cols=NULL), # list with 2 parameters or NULL for nothing done at all
		# for rows: alpha, effect, hclust
		
		colorsColumns=NULL, # by default will be red  
		
		colorsColumnsBy=NULL,  # use colorsColumnsByPalette
		colorsColumnsByPalette=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),# from brewer.pal(8,'accent') -- package RColorBrewer, 
		
		colorsUseMeanQuantiles = FALSE, #### If TRUE colors of bars will depends on following argument
		# test wether column quantilesColors is available in expressionSet object
		colorsMeanQuantilesPalette = c('orange','red','darkred'), # we will use the vector quantilesColors that have to be stored in expressionSet object (first use an helper function to compute quantiles on whole object)
		
		#<todo>
		colorsBarsMatrix = NULL, ## could be provided by user  -- overides colors colorsColumns and colorsColumnsBy
		colorsGenesNames = c('black'),  # colors to be used to give a color to genes names -- are recycled
		# warning: colors must be provided for genes in same order than in raw data as we don't know how they will be sorted...
		
		main = paste("log2 ratio's"),
		shortvarnames = NULL, # specify a variable from pData that will be used to display names 
		longvarnames = NULL, # if NULL, then .withingroups will be used
		gene.length = 50, 
		gene.fontsize = 6, 
		main.fontsize = 9,
		columnhead.fontsize = 8,
		mx = 1.5,
		exp.width = 1.8,
		exp.height = 0.2, 
		log2l.show = TRUE,
		log4l.show = FALSE,
		quantiles.show = FALSE,
		quantiles.compute = c(0.9),
		error.show = TRUE,
		view.psid = FALSE, 
		errorLabel = "Error bars show the pooled standard deviation", 
		closeX11=FALSE,
		openFile=FALSE,
		tooltipvalues=FALSE,
		probe2gene = TRUE,
		...){
	
	stopifnot(is(e, "ExpressionSet") | is(e, "ExpressionSetWithComputation"))
	
	if(!is.null(orderBy$cols)){
		if(!orderBy$cols %in% colnames(pData(e))){
			stop( "The column name by which the data has to be sorted cannot be found in the phenodata!\n")
		}
	}
	
	if (!inherits(e, "ExpressionSetWithComputation")){
		e <- try(computeLogRatio(e,reference=reference,within=within,across=across,nReplicatesVar=nReplicatesVar,...))
		if (inherits(e, "try-error")) 
      stop("Error when doing required computation.")
		
	}  
	device <- tolower(device)
	if (!(device %in% c('svg','pdf','x11','png','cairopng', 'javagd'))) 
    stop("Unknown device - use one of: svg, pdf, x11, png, javagd")
	if (device == "svg"){
		stopifnot(requireNamespace("gridSVG"))
		stopifnot(interactive())
	}
	if (device == "cairopng"){
    	stopifnot(requireNamespace("Cairo"))
  	}
	if(device == "javagd"){
		stopifnot(requireNamespace("JavaGD"))
	}
	
	# we now work with the complete ExpressionSet object enriched with statistics (name: e)
	# there are two parts: data with averages, variances and so on and pData to handle metadata
	# we try to work as long as possible with the entire object with ExpressionSet class
	# as it automatically subset correctly rows and cols
	
	igenes <- featureNames(e)
	if (probe2gene){
		if(is.null(featureData(e)$ENTREZID) | is.null(featureData(e)$SYMBOL) | is.null(featureData(e)$GENENAME)){
			e <- addGeneInfo(e)
		}
		igenes.ll <- featureData(e)$ENTREZID
		igenes.symbol <- featureData(e)$SYMBOL
		igenes.name <- featureData(e)$GENENAME
		igenes.name <- paste(igenes.symbol, igenes.name, sep = " - ")
	} else {
		igenes.ll <- rep(NA,length(igenes))
		igenes.symbol <- rep(NA,length(igenes))
		igenes.name <- igenes
	}
	# if annotation is used, ensure we select only found genes
	e <- e[igenes, ]
	
	# to be able to sort short/long varnames according to columns reordering (if those arguments are provided)
	if (!is.null(shortvarnames)) names(shortvarnames)<-rownames(pData(e)[e$statistic=='diffref',])
	if (!is.null(longvarnames)) names(longvarnames)<-rownames(pData(e)[e$statistic=='diffref',])
	
	# reorder columns/rows 
	if (!is.null(orderBy)){
		if (is.list(orderBy)){
			if (!is.null(orderBy$cols)) 
				colorder <- do.call('order',
						as.list(as.data.frame(pData(e)[,orderBy$cols])))    else colorder <- 1:ncol(exprs(e))
			if (!is.null(orderBy$rows)){
				if (is.numeric(orderBy$rows))  {
					roworder <- orderBy$rows
				} else {
					if (orderBy$rows=='alpha'){  
						roworder <- order(igenes.name)
					} else if (orderBy$rows == 'effect'){
						if(sum(e$statistic == 'diffref') == 1){roworder <- order(exprs(e)[, e$statistic == 'diffref'], decreasing = TRUE)
						} else {roworder <- order(apply(exprs(e)[, e$statistic == 'diffref'], 1, FUN = function(vec){
												sum(vec)
											}), decreasing = TRUE)}
					} else if (orderBy$rows=='hclust'){
						hc=hclust(dist(exprs(e[,e$statistic=='diffref'])))
						roworder=hc$order 
					}        #D print('OK')
				}} else roworder <- 1:nrow(exprs(e))
		} else roworder <- 1:nrow(exprs(e))
		
	} else{ 
		roworder <- 1:nrow(exprs(e))
		colorder <- 1:ncol(exprs(e))
		
	}
# reorder ExpressionSet object
	e <- e[roworder,colorder]
# also order information on genes: names, colors,...
	igenes.name <- igenes.name[roworder]
	if(exists("igenes.ll")){
		igenes.ll <- igenes.ll[roworder]}
	if(exists("igenes.symbol")){
		igenes.symbol <- igenes.symbol[roworder]}
	
	
# genes names colors (recycle colors) -- sort them on roworders ! 
	colorsGenesNames <- rep(colorsGenesNames,nrow(exprs(e)))[1:nrow(exprs(e))][roworder]
	
	e.diff <- e[,e$statistic=='diffref']
	
# to be used in graph
	e.quantiles <- matrix(apply(exprs(e.diff), 2,
					function(tr) quantile(tr,probs=quantiles.compute)),nrow=length(quantiles.compute),byrow=FALSE)
	
	nc <- ncol(e.diff) # number of columns
	nr <- nrow(e.diff)
	mx <- max(abs(exprs(e.diff))) / mx
	
	
	if (is.null(shortvarnames)) shortvarnames <- pData(e.diff)$.withingroups # use adequate order as ExpressionSet object is sorted
	else shortvarnames <- shortvarnames[rownames(pData(e.diff))] #reorder 
	if (is.null( longvarnames))  longvarnames <- shortvarnames 
	else longvarnames <- longvarnames[ rownames(pData(e.diff))]      #reorder
	
	
	
	### management of colors 
	
	
# columns headers
	## order colorsColumns
	if(!is.null(colorsColumns) & all(names(colorsColumns) %in% shortvarnames)){
		colorsColumns <- colorsColumns[shortvarnames]
	}
	if (is.null(colorsColumns)|(!is.null(colorsColumnsBy))){
		if (!is.null(colorsColumnsBy)){
			uniquecolors <- as.numeric(as.factor(
							do.call('paste',c(pData(e.diff)[,colorsColumnsBy,drop=FALSE],sep='.')))) 
			if (is.null(colorsColumnsByPalette)) colorsColumns <- rainbow(max(uniquecolors))[uniquecolors]
			else colorsColumns <- rep(colorsColumnsByPalette,max(uniquecolors))[uniquecolors]
		}
		else colorsColumns <- rep('red',nc)
	}  else colorsColumns <- rep(colorsColumns,nc)[1:nc]  
	
# colorsColumns will be used for headers (by default, all red)
	
	
# bars colors 
	if (colorsUseMeanQuantiles | !is.null(colorsBarsMatrix)){
		# use colorsBarsMatrix if not NULL, else create it using the vector of colors previously built
		if (!is.null(colorsBarsMatrix)){
			# check for enough information in matrix
			stopifnot((nrow(colorsBarsMatrix)>=nr)|(ncol(colorsBarsMatrix)>=nc))
			
		}
		else {
			# use colorsMeanQuantiles recycled for each column
			# have to get genes colors from featuresData slot of ExpressionSet object
			stopifnot( "colorsQuantilesVector" %in% colnames(fData(e)))
			colorsBarsGenes <- fData(e)[,'colorsQuantilesVector']
			stopifnot(length(colorsMeanQuantilesPalette)>=max(colorsBarsGenes))
			colorsBarsMatrix <- matrix(rep(matrix(colorsMeanQuantilesPalette[colorsBarsGenes]),nc),ncol=nc)
			# we don't have to sort anymore colorsBarsMatrix as fData are already sorted
			roworder <- 1:nrow(exprs(e))
		}
	} else
	{
		colorsBarsMatrix <-  matrix(rep(colorsColumns,nr),ncol=nc,byrow=TRUE)
		# no argument specific on bars, duplicates columns headers for each bar    
		# recycle  colorsColumns vector
	}
#reorder according to genes order
	colorsBarsMatrix <- colorsBarsMatrix[roworder,,drop=FALSE]  
	
	if (rev(strsplit(filename,split="")[[1]])[4] != '.'){
		if (device=="cairopng") filename <- paste(filename,'png',sep='.')
		filename <- paste(filename,device,sep='.')
	}
	if (dirname(filename)=='.'){
		filename <- file.path(getwd(),filename) # just to ensure viewer will find it
	}
	
# open window of specified size using either X11 of PDF device
	if (device %in% c('svg', 'x11')) 
		x11(width = 3 + (nc * exp.width), height = nr * exp.height,...) else {
		if (device=='javagd') JavaGD::JavaGD(width = (3 + (nc * exp.width))*100, height = (nr * exp.height)*100) ## FLAG
		if (device=='pdf') pdf(filename, width = 3 + (nc * exp.width), height = nr * exp.height,...)
		if (device=='png') png(filename, width = round(100*(3 + (nc * exp.width))), height = round(100*(nr * exp.height)),...)
		if (device=='cairopng')   Cairo::CairoPNG(filename, width = round(100*(3 + (nc * exp.width))), height = round(100*(nr * exp.height)),...)
	}
	##E width=min(xx,maxwidth)
	
# create main viewport
	vp1 <- viewport(x = 0, y = 0, width = 1, height = 1, just = c("left", "bottom"),
			layout = grid.layout(nrow = nr + 14, ncol = 1))
	pushViewport(vp1)
	
# draw top section of the graph (vp2)
	vp2 <- viewport(layout.pos.row = 1:6)
	pushViewport(vp2)
	
# draw gray box on top
	grid.rect(x = 0, y = 0, height = 1, width = 1, just = c("left", "bottom"),
			gp = gpar(fill = "lightgray"))
	
# add title
	grid.text(x = 0.5, y = 0.9, just = c("center", "top"), 
			gp = gpar(fontsize = main.fontsize), label = main)
	
# I reserve the first three columns of the graph for displaying the names
# of the genes, hence (i + 3)
# since I need one empty column on the right side I have a total of (nc + 4) columns
	
	
	
	
	
	for (i in 1:nc) {
		pData.i <- lapply(pData(e.diff)[i,],as.character)
		
		grid.text(x = (i + 3) / (nc + 4), y = 0.2, 
				label = shortvarnames[i],               ### color to be changed (first color)
				gp = gpar(col = colorsColumns[i], fontsize = columnhead.fontsize), just = c("center", "center"),
				name=paste("treatment_",i,sep=""))
		
		if (device=='svg'){
			pDatacol <- colnames(pData(e))[-grep('^\\.',colnames(pData(e)))]
			pDatacol <- pDatacol[pDatacol != 'statistic']
			tmp <- pData(e.diff)[,pDatacol]
			for (cn in pDatacol){
				tmp[,cn] <- paste(cn,tmp[,cn],sep='=')
			}
			details <- apply(tmp,1,paste,collapse=';')
			gridSVG::grid.garnish(paste("treatment_",i,sep=""),
					onmouseover=paste("highLightTreatment(evt,'",longvarnames[i],"', '",details[i],"');",sep=""),
					onmouseout="toolTip.setAttributeNS(null, 'display', 'none');" )
		}
	}
	popViewport()
	
# draw main section of the graph (vp3)
	vp3 <- viewport(layout.pos.row = 9:(nr + 9))
	pushViewport(vp3)
	for (i in 1:abs(nr / 2)) {
		grid.rect(x = 0, y = (1 / nr) * (i * 2),
				height = (1 / nr), width = 1,
				just = c("left", "centre"),
				gp = gpar(fill = colors()[63], col = 'white'))
	}
# define a vector with zeros to place multiple objects efficiently later on
	j <- rep (0, nr)
	
# main loop for placing all elements in each column
	for (i in 1:nc) {
		###E here there is a way to use a mxi rescale parameter individully for each
		### column instead of mx -- mxi to be computed using max(col)
		# make a gray box which marks 1 log2, so that scientists later can easily
		# see, how much (how many log2s) a gene was up- or down-regulated
		grid.rect(x = (i + 3) / (nc + 4), y = 0 + (1 / (nr * 2)),
				height = 1, width = (1 / (nc + 4)) / mx,
				just = c("centre", "bottom"),
				gp = gpar(fill = "lightgray", col = NULL))
		
		# draw the zero line - boxes right of this line indicate an up regulation
		# of a gene while boxes to the left indicate a down-regulation
		grid.lines(x = c((i + 3) / (nc + 4), (i + 3) / (nc + 4)),
				y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
				gp = gpar(col = colorsColumns[i], lwd = 1.5))
		
		if ((log2l.show == TRUE) & (max(abs(exprs(e.diff)[,i]))>2)) {
			# draw the first 2 log2 line 
			grid.lines(x = c((i + 3) / (nc + 4) - (2 / (nc + 4)) / (mx * 2),
							(i + 3) / (nc + 4) - (2 / (nc + 4)) / (mx * 2)),
					y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
					gp = gpar(col = "lightgrey", lwd = 1))
			# draw the second 2 log2 line 
			grid.lines(x = c((i + 3) / (nc + 4) + (2 / (nc + 4)) / (mx * 2),
							(i + 3) / (nc + 4) + (2 / (nc + 4)) / (mx * 2)),
					y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
					gp = gpar(col = "lightgrey", lwd = 1))  
		}
		
		if ((log4l.show == TRUE)& (max(exprs(e.diff)[,i])>4)) {
			# draw the first 4 log2 line 
			grid.lines(x = c((i + 3) / (nc + 4) - (4 / (nc + 4)) / (mx * 2),
							(i + 3) / (nc + 4) - (4 / (nc + 4)) / (mx * 2)),
					y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
					gp = gpar(col = "lightgrey", lwd = 1))
			# draw the second 4 log2 line 
			grid.lines(x = c((i + 3) / (nc + 4) + (4 / (nc + 4)) / (mx * 2),
							(i + 3) / (nc + 4) + (4 / (nc + 4)) / (mx * 2)),
					y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
					gp = gpar(col = "lightgrey", lwd = 1))  
		}
		
		
		
		if (quantiles.show) {
			# draw lines for computed quantiles 
			tmp <- sapply(e.quantiles[,i],function(qq) {
						grid.lines(x = c((i + 3) / (nc + 4) + (qq / (nc + 4)) / (mx * 2),
										(i + 3) / (nc + 4) + (qq / (nc + 4)) / (mx * 2)),
								y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
								gp = gpar(col = "lightgrey", lwd = 1))  
						
					})
			
		}
		
		# draw box according to the size of the ratio in one go for all rows
		# here comes the vector j in...
		# and it took me long to figure this one out properly...
		#E I change it to a loop to be able to add tooltips
		for (jj in 1:nrow(e.diff)){
			
#D cat('\n jj:',jj,' i:',i)
#D cat(' color:', colorsBarsMatrix[jj,i])
#D cat(' name in matrix',rownames(e.diff[jj]), ' full name:',igenes.name[jj] )
			
			grid.rect(x = (i + 3) / (nc + 4) +  
							( ( exprs(e.diff)[jj,i]  ) * ((1 / (nc + 4)) / (mx * 2)) ) / 2,
					#       y = (1 / nr)*jj, 
					y = (nr-jj+1)/nr, 
					height = (1 / (nr * 2)),
					width = abs((  exprs(e.diff)[jj,i] ) * ((1 / (nc + 4)) / (mx * 2))),
					just = c("centre", "centre"), ### change to use matrix
					
					gp = gpar(fill = colorsBarsMatrix[jj,i], col = NULL),name=paste('val',jj,i,sep='_') )
			
			if (device=='svg' & tooltipvalues){
				
				gridSVG::grid.garnish(paste('val',jj,i,sep='_'),
						onmouseover=paste("highLightTreatment(evt,'",round(exprs(e.diff)[jj,i],1),"', '",details[i],"');",sep=""),
						onmouseout="toolTip.setAttributeNS(null, 'display', 'none');" )
			}
			
		}  
		# draw error bar
		# only when available...
		if (error.show) {
			# retrieve column for pooled sd IF available
			#   withinData <- split(pData(e),pData(e)$.withingroups)[[i]]
			#   withinData <- withinData[withinData$statistic=='signedpooledSD',]
			errorname <- paste(
					substr(rownames(pData(e.diff))[i],1,nchar(rownames(pData(e.diff)))[i]-8)
					,'spooledSD',sep='.')
			
			
			#  if (nrow(withinData)==1){
			if (errorname %in% rownames(pData(e))){
#          psd.i <- exprs(e)[,rownames(withinData)]
				psd.i <- exprs(e)[,errorname]
				
				# horizontal line of the error bar
				# data.sd contains the data for the adjusted standard error
				# from the Dunnett's test
				grid.segments(x0 = (j + i + 3) / (nc + 4) + (exprs(e.diff)[,i]) 
								* ((1 / (nc + 4)) / (mx * 2)),
						x1 = (j + i + 3) / (nc + 4) + (exprs(e.diff)[,i]+ psd.i) 
								* ((1 / (nc + 4)) / (mx * 2)),
						#  y0 = (1 / nr) * seq(nr),
						#  y1 = (1 / nr) * seq(nr),
						
						y0 = (1 / nr) * rev(seq(nr)),
						y1 = (1 / nr) * rev(seq(nr)),
						
						
						gp = gpar(lwd = 0.5, col = "black") )
				# vertical line of the error bar
				grid.segments(x0 = (j + i + 3) / (nc + 4) +
								(exprs(e.diff)[,i] + psd.i) * ((1 / (nc + 4)) / (mx * 2)),
						x1 = (j + i + 3) / (nc + 4) + (exprs(e.diff)[,i] +psd.i) * ((1 / (nc + 4)) / (mx * 2)),
#            y0 = (1 / nr) * seq(nr) - (1 / (nr * 4)),
#            y1 = (1 / nr) * seq(nr) + (1 / (nr * 4)),
						y0 = (1 / nr) * rev(seq(nr)) - (1 / (nr * 4)),
						y1 = (1 / nr) * rev(seq(nr)) + (1 / (nr * 4)),
						
						gp = gpar(lwd = 0.5, col = "black") )
			}    
		}
	}
	
  # draw the names of the genes and hyperlink them to the current annotation
  # igenes.name contains a vector with the names of all the "interesting" genes
  # since they are often too long, I cut them off
  # igenes.ll contains the number for the entrez gene annotation
  # both data are all provided by the bioconductor project
	g <- substring(igenes.name, 1, gene.length)
  #  print(g)
	for (i in 1:nr) {
		grid.text(x = 0.005, 
				#y = (i/ nr), 
				y=(nr-i+1)/nr, 
				label = as.character(g[i]), name = paste("gene", i, sep = ""),
				just = c("left", "centre"), gp = gpar(fontsize = gene.fontsize, col = colorsGenesNames[i]))
		
		if (!probe2gene & device=='svg'){
			gridSVG::grid.hyperlink(gPath(paste("gene", i, sep = "")),  paste("http://www.ncbi.nih.gov/entrez/query.fcgi?db=gene&amp;cmd=Retrieve&amp;dopt=summary&amp;list_uids=",as.integer(igenes.ll[i]), sep = "") )
		}
	}
	popViewport()
	
	
  # draw bottom section of the graph (vp4)
	vp4 <- viewport(layout.pos.row = (nr + 11):(nr + 14))
	pushViewport(vp4)
	
  # draw gray box at the bottom
	grid.rect(x = 0, y = 0, height = 1, width = 1, just = c("left", "bottom"),
			gp = gpar(fill = "lightgrey"))
	
  # the following commands simply draw the legend
	###<CHANGE>
	if (error.show) grid.text(x = 0.05, y = 0.5, 
				just = c("left", "bottom"), gp = gpar(fontsize = 9),label = errorLabel)
	if (colorsUseMeanQuantiles){
		grid.text(x = 0.8, y = 0.5, 
				just = c("right", "bottom"), gp = gpar(fontsize = 9),label = "Increasing quantiles ")
		k=length(colorsMeanQuantilesPalette)
		for (ik in 1:k){
			grid.rect(x = 0.8+(ik-1)*((0.9-0.8)/k),
					y = 0.5, 
					height = 0.2,
					width = (0.9-0.8)/k,
					just = c("left", "bottom"),gp = gpar(fill =colorsMeanQuantilesPalette[ik] , col = 'black'))
			
		}
	}
	
	grid.text(x = 0.05, y = 0.05, just = c("left", "bottom"),
			gp = gpar(fontsize = 4, col = "black"),
			label = paste(date(),";",R.version.string,"; Biobase version ", 
					package.version("Biobase")))
	
	popViewport()
	
  # finish up the graph
  # just a box around the complete graph
	grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar(fill = NA, col = "black"))
	popViewport()
	
	if (device == 'svg') {        #E includes script for tooltips
		gridSVG::grid.script(file=file.path(path.package('a4'), 'etc', 'tooltip.script'))  
		gridSVG::gridToSVG(name = filename)
		if (openFile) browseURL(url = filename)
	}
	if ((device %in% c('pdf','png','cairopng')) | (device %in% c('svg','x11') & closeX11)){
		if(device %in% 'svg'){
			gridSVG::dev.off()
		}else	grDevices::dev.off()
	}
	if (device == "pdf" & openFile)  openPDF(filename)
	
	invisible(e)       
	
}