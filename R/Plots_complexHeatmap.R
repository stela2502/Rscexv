#' @name complexHeatmap
#' @aliases complexHeatmap,Rscexv-method
#' @rdname complexHeatmap-methods
#' @docType methods
#' @description plot the PCR heatmap using the heatmap.3 function included in this package 
#' @param x the Rscexv object
#' @param ofile the outfile to create in the x@outpath folder
#' @param colGroups columns in the samples table to use to order the data (first == order)
#' @param rowGroups rows in the annotation table to use to color the heatmap rows (first == order)
#' @param colColors a named list of column color vectors
#' @param rowColors a named list of row color vectors
#' @param pdf export as pdf (default = FALSE)
#' @param subpath the subpath for the plots (default = '')
#' @param heapmapCols the color function to calculate the heatmap colours ( default function (x) { c("darkgrey",bluered(x)) } )
#' @title description of function plot.beans
#' @export 
setGeneric('complexHeatmap', ## Name
		function ( x,  ofile=NULL, colGroups=NULL, rowGroups=NULL, colColors=NULL, rowColors=NULL, pdf=FALSE, subpath='', main = '',  heapmapCols= function(x){ c("darkgrey",bluered(x))} ) { ## Argumente der generischen Funktion
			standardGeneric('complexHeatmap') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('complexHeatmap', signature = c ('Rscexv'),
		definition = function ( x,  ofile=NULL, colGroups=NULL, rowGroups=NULL, colColors=NULL, rowColors=NULL, pdf=FALSE, subpath='', main = '' ,  heapmapCols= function(x){ c("darkgrey",bluered(x))} ) {
			
			Rowv = FALSE
			Colv = FALSE
			dendrogram = 'both'
			ColSideColors <- NULL
			RowSideColors <- NULL
			ColSideColorsSize <- 1
			RowSideColorsSize <- 1
			if ( is.null(colColors) ){
				colColors <- x@usedObj[['colorRange']]
			}
			if ( is.null(rowColors) ){
				rowColors <- x@usedObj[['colorRange']]
			}
			if ( ! is.null(colGroups) ) {
				ColSideColorsSize <- length(colGroups)
				x <- reorder.samples(x, colGroups[1] )
				for ( i in colGroups ){
					if ( is.na(match( i, names(colColors))) ){
						stop( paste( "No colours for the grouping", i ) )
					}
					ColSideColors <- cbind(ColSideColors, colColors[[ match( i, names(colColors)) ]][x@samples[, i]] )
				}
				colnames(ColSideColors) = colGroups
				#ColSideColors <- matrix( ColSideColors, ncol= ColSideColorsSize)
				Colv = FALSE
				if ( !is.null(rowGroups)){
					dendrogram = 'none'
				}else{
					dendrogram= 'none'
				}
			}else {
				## probably calculate the clustering??
			}
			if ( ! is.null(rowGroups) ) {
				RowSideColorsSize <- length(rowGroups)
				x <- reorder.genes(x, rowGroups[1] )
				for ( i in rowGroups ){
					RowSideColors <- rbind( RowSideColors,rowColors[[ match( i, names(rowColors)) ]][x@annotation[, i]] )
				}
				rownames(RowSideColors) = rowGroups
				Rowv = FALSE
				#RowSideColors <- matrix( RowSideColors, nrow= RowSideColorsSize)
				if ( !is.null(colGroups)){
					dendrogram = 'none'
				}else{
					dendrogram= 'none'
				}
			}else {
				## probably calculate the clustering??
			}
			data <- as.matrix(t(x@data))
			brks <- unique(as.vector(c(-20,quantile(data[which(data!= -20)],seq(0,1,by=0.1)),max(data))))
			if ( ! is.null(ofile)){
				if ( pdf ) {
					pdf( file=paste(file.path(x@outpath,ofile),'pdf',sep='.'), width=10, height=8)
				}else{
					png( file=paste(file.path(x@outpath,ofile),'png',sep='.'), width=1600, height=800)
				}
			}
			heatmap.3(
					data, breaks=brks,col=heapmapCols(length(brks)-2), Rowv= is.null(RowSideColors), Colv = is.null(ColSideColors),  key=F, symkey=FALSE,
					trace='none', 
					ColSideColors=ColSideColors,ColSideColorsSize=ColSideColorsSize, 
					RowSideColors=RowSideColors,RowSideColorsSize=RowSideColorsSize, 
					cexRow=0.6,cexCol=0.7,main=main, dendrogram=dendrogram, labCol = "", lwid=c(0.5,4), lhei=c(1,4)
			)
			if ( ! is.null(ofile)){
				dev.off()
				pdf( file=paste(file.path(x@outpath,ofile),'_legend_values.pdf',sep='.'), width=8, height=4)
				Z <- as.matrix(1:(length(brks)-2))
				image(Z, col=heapmapCols(length(brks)-2),axes = FALSE, main='color key')
				axis( 1, at=c(0,0.1,1), labels=c('NA','low','high'))
				dev.off()
				for ( cname in c( rowGroups, colGroups) ) {
					plot.legend( x, cname, file= ofile, pdf= pdf )
				}
			}
			
		}
)

