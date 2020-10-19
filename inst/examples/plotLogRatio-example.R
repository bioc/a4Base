if (require(ALL)){
	data(ALL, package = "ALL")
	ALL <- addGeneInfo(ALL)
	ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
	ALL2 <- ALL[,ALL$BT != 'T1']  # omit subtype T1 as it only contains one sample
	ALL2$BTtype <- as.factor(substr(ALL2$BT,0,1)) # create a vector with only T and B
	
	# Test for differential expression between B and T cells
	tTestResult <- tTest(ALL, "BTtype", probe2gene = FALSE)
	topGenes <- rownames(tTestResult)[1:20]
	
	# plot the log ratios versus subtype B of the top genes 
	LogRatioALL <- computeLogRatio(ALL2, reference=list(var='BT',level='B'))
	a <- plotLogRatio(e=LogRatioALL[topGenes,],openFile=FALSE, tooltipvalues=FALSE, device='pdf',
			colorsColumnsBy=c('BTtype'), main = 'Top 20 genes most differentially between T- and B-cells',
			orderBy = list(rows = "hclust"), probe2gene = TRUE)
\dontrun{		
	a <- plotLogRatio(e=LogRatioALL[topGenes,],openFile=TRUE, tooltipvalues=FALSE, device='pdf',
			colorsColumnsBy=c('BTtype'), main = 'Top 20 genes most differentially between T- and B-cells',
			orderBy = list(rows = "hclust", cols = "sex"), probe2gene = TRUE)
}
}