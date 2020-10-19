\dontrun{
	library(RColorBrewer)
	library(dichromat)
	
	library(Biobase)
	library(grid)
	pdf.directory=getwd()
	
	
	load(file.path(getwd(),"expressionSetRma.Rda"))      #expressionSetRma
	
	
	eset <- expressionSetRma[100:130,pData(phenoData(expressionSetRma))[,"sample"]%in%c(1:10,41:50)] # ARG
	##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	exprs(eset)[1,5] <- 13.8
	exprs(eset)[10,7] <- 0.5
	eset2 <- expressionSetRma[200:250,] # ARG
	eset3 <- expressionSetRma[1000:1009,pData(phenoData(expressionSetRma))[,"sample"]%in%c(1:3,41:46)] # ARG
	eset4 <- expressionSetRma[100:230,pData(phenoData(expressionSetRma))[,"sample"]%in%c(1:20,31:50)] # ARG
	
	eset5 <- expressionSetRma[1:400,] # ARG
	
	# eset <- eset2
	
	pdf(file.path(pdf.directory,"eset.pdf"))
	size <- heatmap.expressionSet(eset,subtitle.main=" ")
	dev.off()
	pdf(file.path(pdf.directory,"eset.pdf"),width=size[1],height=size[2])
	heatmap.expressionSet(eset,subtitle.main=" ")
	dev.off()
	
	
	pdf(file.path(pdf.directory,"eset2.pdf"))
	size <- heatmap.expressionSet(
			eset2,
			colors.nbreaks = 20,
			colors.pergroup=TRUE,
			legend.range="data",
			row.col.groups.display=FALSE,
			cell.gpar=gpar(lwd=0.5),
			legend.height=unit(50,"points"),
			title.just=c("center","center"),
			title.maxlines=2,
			col.groups.sep.width=unit(0,"points"),
			row.labels=featureNames(eset),
			subtitle.main="This is subtitle",
			row.order="hclust",row.groups.hclust=FALSE,
			title.gpar=gpar(cex=2),
			subtitle.gpar=gpar(cex=1.5)
	)
	dev.off()
	pdf(file.path(pdf.directory,"eset2.pdf"),width=size[1],height=size[2])
	size <- heatmap.expressionSet(
			eset2,
			colors.nbreaks = 20,
			colors.pergroup=TRUE,
			legend.range="data",
			row.col.groups.display=FALSE,
			cell.gpar=gpar(lwd=0.5),
			legend.height=unit(50,"points"),
			title.just=c("center","center"),
			title.maxlines=2,
			col.groups.sep.width=unit(0,"points"),
			row.labels=featureNames(eset),
			subtitle.main="This is subtitle",
			row.order="hclust",row.groups.hclust=FALSE,
			title.gpar=gpar(cex=2),
			subtitle.gpar=gpar(cex=1.5)
	
	)
	dev.off()
	
	
	
	
	
	pdf(file.path(pdf.directory,"eset3.pdf"))
	size <- heatmap.expressionSet(
			eset3,
			row.labels.gpar=gpar(cex=0.4,col=c(rep("red",2),rep("black",49))	), # col will correctly be a vector only if no group...
			col.labels.gpar=gpar(cex=0.6),
			colors.nbreaks = 20,
			colors.pergroup=TRUE,
			legend.range="data",
			row.col.groups.display=FALSE,
			cell.gpar=gpar(lwd=0.5),
			legend.height=unit(50,"points"),
			title.just=c("center","center"),
			title.maxlines=2,
			col.groups.sep.width=unit(0,"points"),
			row.labels=featureNames(eset),
			subtitle.main="Essai subtitle",
			row.order="hclust",row.groups.hclust=FALSE,
			interactive=FALSE
	)
	dev.off()
	
	pdf(file.path(pdf.directory,"eset3.pdf"),width=size[1],height=size[2])
	size <- heatmap.expressionSet(
			eset3,
			row.labels.gpar=gpar(cex=0.4,col=c(rep("red",2),rep("black",49))	), # col will correctly be a vector only if no group...
			col.labels.gpar=gpar(cex=0.6),
			colors.nbreaks = 20,
			colors.pergroup=TRUE,
			legend.range="data",
			row.col.groups.display=FALSE,
			cell.gpar=gpar(lwd=0.5),
			legend.height=unit(50,"points"),
			title.just=c("center","center"),
			title.maxlines=2,
			col.groups.sep.width=unit(0,"points"),
			row.labels=featureNames(eset),
			subtitle.main="Essai subtitle",
			row.order="hclust",row.groups.hclust=FALSE,
			interactive=FALSE
	)
	dev.off()
	
	
	
	pdf(file.path(pdf.directory,"eset4.pdf"))
	size <- heatmap.expressionSet(
			eset4,
			legend.range="data",
			colors.palette = dichromat(rich.colors(190)[1:128]),
			row.col.groups.display=TRUE,
			title.just=c("left","top"),
			title.maxlines=2,
			row.labels=featureNames(eset),
			subtitle.main="",
			row.order="hclust",row.groups.hclust=FALSE,
	)
	dev.off()
	
	pdf(file.path(pdf.directory,"eset4.pdf"),width=size[1],height=size[2])
	size <- heatmap.expressionSet(
			eset4,
			legend.range="data",
			colors.palette = dichromat(rich.colors(190)[1:128]),
			row.col.groups.display=TRUE,
			title.just=c("left","top"),
			title.maxlines=2,
			row.labels=featureNames(eset),
			subtitle.main="",
			row.order="hclust",row.groups.hclust=FALSE,
	)
	dev.off()
	
	pdf(file.path(pdf.directory,"eset5.pdf"))
	size <- heatmap.expressionSet(eset5,row.order="hclust",row.groups.hclust=FALSE)
	dev.off()
	
	pdf(file.path(pdf.directory,"eset5.pdf"),width=size[1],height=size[2])
	heatmap.expressionSet(eset5,row.order="hclust",row.groups.hclust=FALSE)
	dev.off()
	
}