### utility function to generate an ExpressionSet object
### from expression data and phenodata
### 

#' combine gene expression and phenotype data onto a ExpressionSet object
#' 
#' Basically a wrapper for \code{new('ExpressionSet',...)}, this function gathers gene
#' expression and phenotype data, after having checked their compatibility.
#' @param exprs gene expression matrix
#' @param phenoData phenotype data associated with exprs columns, as a matrix or data.frame
#' @param varMetadata optional metadata on phenotype data
#' @param dimLabels see \code{\link[Biobase]{ExpressionSet}}
#' @param featureData see \code{\link[Biobase]{ExpressionSet}}
#' @param experimentData see \code{\link[Biobase]{ExpressionSet}}
#' @param annotation see \code{\link[Biobase]{ExpressionSet}}
#' @param changeColumnsNames Change exprs columns names -- see details
#' @param ... \code{\dots}
#' @details
#'  If \code{changeColumnsNames} is \code{TRUE}, then the procedure is the following: first one checks if phenoData contains a column named 'colNames'. If so, content will be used to rename exprs colums. On the other case, one uses combinations of phenoData columns to create new names. In any case, old columns names
#'  are stored within a column named 'oldcolnames' in the pData.
#' @examples 
#' # simulate expression data of 10 features (genes) measured in 4 samples
#' x <- matrix(rnorm(40), ncol = 4)
#' colnames(x) <- paste("sample", 1:4, sep = "_")
#' rownames(x) <- paste("feature", 1:10, sep = "_")
#' # simulate a phenodata with two variables
#' ToBePheno <- data.frame(Gender = rep(c('Male','Female'), 2), 
#' 		Treatment = rep(c('Trt','Control'), each=2))
#' rownames(ToBePheno) <- paste("sample", 1:4, sep = "_")
#' eset <- createExpressionSet(exprs = x, phenoData = ToBePheno)
#' @return An object of class ExpressionSet
#' @author Eric Lecoutre
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @keywords data
#' @importFrom Biobase AnnotatedDataFrame ExpressionSet MIAME
#' @export
createExpressionSet <-
		function(
				exprs = matrix(nrow = 0, ncol = 0),
				phenoData = AnnotatedDataFrame(),
				varMetadata= NULL,
				dimLabels = c("rowNames","colNames"),
				featureData = NULL, 
				experimentData = MIAME(), 
				annotation = character(0), 
				changeColumnsNames = TRUE,...){  
	

	if (nrow(phenoData) != ncol(exprs)) 
    stop('phenoData must have the same number of rows than exprs number of columns')  

	if (all(rownames(phenoData) != colnames(exprs))){
		stop("rownames of phenoData are not identical to colnames of exprs.\n")
	}
	
	if (!is(phenoData, "AnnotatedDataFrame")){
		# we must prepare AnnotatedDataFrame
		# check varMetadata consistency with phenoData matrix
		if (!is.null(varMetadata)){
			if ((nrow(varMetadata)!=ncol(phenoData)) | (!'labelDescription' %in% colnames(varMetadata))){
				warning(paste("varMetadata not compliant with phenoData", "maybe there is not a column called 'labelDescription'","check ?AnnotatedDataFrame","--- we will not use it",sep='\n'))
				phenoData <- AnnotatedDataFrame(data=phenoData,dimLabels=dimLabels)
			}
			else {
				phenoData <-AnnotatedDataFrame(
						data=phenoData,
						varMetadata=varMetadata,
						dimLabels=dimLabels)
			}
		} else phenoData <- AnnotatedDataFrame(
					data=phenoData,
					dimLabels=dimLabels)
	}
	
	if (changeColumnsNames){
		oldcolnames <- colnames(exprs)
		if ("colNames" %in% colnames(phenoData)){
			if (any(duplicated(phenoData$colNames))) {
				warnings("Cant' use 'colNames' as new colnames as it has duplicates -- we use V1-Vn new names")
				colnames(exprs) <- paste('V',1:ncol(exprs),sep='')
			}
			else 
				colnames(exprs) <- phenoData$colNames
		}
		if ((!"colNames" %in% colnames(phenoData))){
			newcolnames <- do.call('paste',c(as.list(pData(phenoData)),sep='.'))
			if (any(duplicated(newcolnames))) {
				newcolnames <- paste(newcolnames,replicates(newcolnames),sep='.')            
				colnames(exprs) <- newcolnames
			}
		}
		
		phenoData$.oldcolnames <- oldcolnames
	}
	rownames(pData(phenoData)) <- colnames(exprs)  
	if (is(exprs,'data.frame')) exprs <- as.matrix(exprs)
	if (!is.null(featureData)){
		out <- ExpressionSet( 
				exprs=exprs,
				phenoData = phenoData,
				featureData=featureData, 
				experimentData = experimentData, 
				annotation = annotation)
	}
	else {
		out <- ExpressionSet(
				exprs=exprs,
				phenoData = phenoData,
				experimentData = experimentData, 
				annotation = annotation)
	}
	return(out)
}


#createExprSet <- function(ExprData, PhenoData){
#	require(affy)
#	if(nrow(PhenoData) != ncol(ExprData)){
#		stop("The number of rows in PhenoData is not equal to the number of 
#			columns in the ExprData.\n")
#	}
#	
#	if(all(rownames(PhenoData) != colnames(ExprData))){
#		stop("rownames of PhenoData are not identical to colnames of ExprData.\n")
#	}
#	myExprSet <- new("ExpressionSet", exprs = ExprData)
#	pData(myExprSet) <- PhenoData
#	return(myExprSet)
#}


### utility function to combine two ExpressionSet objects
### 

#' Combine two ExpressionSet objects
#' 
#' Merge two ExpressionSet objects, checking their attributes.
#' @param x An object of class ExpressionSet
#' @param y An object of class ExpressionSet
#' @details 
#' exprs and pData are merged. 
#' Other data (such as MIAME or annotation) are those of x.
#' @return An object of class ExpressionSet
#' @examples 
#' \dontrun{
#' # prepare and combine two ExpressionSet
#' data(data.H2009); data(phenoData.H2009)
#' data(data.SKOV3); data(phenoData.SKOV3)
#' eH2009 <- prepareExpressionSet(exprs = data.H2009, phenoData = phenoData.H2009, changeColumnsNames = TRUE)
#' eSKOV3  <- prepareExpressionSet(exprs = data.SKOV3, phenoData = phenoData.SKOV3, changeColumnsNames = TRUE)
#' newE <- combineTwoExpressionSet(eH2009,eSKOV3)
#' }
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @author Eric Lecoutre
#' @keywords data
#' @importFrom Biobase AnnotatedDataFrame assayData
#' @export
combineTwoExpressionSet <- function(x,y){
# prioritary  keep information from x, append assayData and phylo data from y
	out <- x
	outAssayData <- new.env()
	assign("exprs",cbind(assayData(x)$exprs,assayData(y)$exprs),
			envir=outAssayData)
	assayData(out) <- outAssayData
	outPhenoData <- AnnotatedDataFrame(
			data=rbind(pData(x),pData(y)),
			varMetadata = varMetadata(x)
	) 
	phenoData(out) <- outPhenoData
	return(out)
}
