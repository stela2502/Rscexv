#' @name PCR.heatmap
#' @aliases PCR.heatmap,Rscexv-method
#' @rdname PCR.heatmap-methods
#' @docType methods
#' @description This function plots the heatmaps
#' @param dataObj  TEXT MISSING
#' @param ofile  TEXT MISSING
#' @param title  TEXT MISSING default='Heatmap'
#' @param nmax  TEXT MISSING default=4000
#' @param hc.row  TEXT MISSING default=NA
#' @param ColSideColors  TEXT MISSING default=NA
#' @param RowSideColors  TEXT MISSING default=F
#' @param PCR.heatmap  TEXT MISSING
#' @param reorder if yes the data matrix will be reordered to the clusters object (default=F)
#' @title description of function PCR.heatmap
#' @export 
setGeneric('PCR.heatmap', ## Name
		function ( dataObj, ofile,reorder =F,  title='Heatmap', nmax=4000, hc.row=NA, ColSideColors=NA, RowSideColors=F,
				width=6, height=6, margins = c(1,11) ,lwid = c( 1,6), lhei=c(1,5), hclustfun = function(c){hclust( c, method='ward.D')}, distfun = function (x) as.dist( 1- cor(t(x), method='pearson') ), Rowv=T, ... ) {## Argumente der generischen Funktion
			standardGeneric('PCR.heatmap') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('PCR.heatmap', signature = c ('Rscexv'),
		definition = function ( dataObj, ofile,reorder =F,  title='Heatmap', nmax=4000, hc.row=NA, 
				ColSideColors=NA, RowSideColors=F, width=6, height=6, margins = c(1,11),
				lwid = c( 1,6), lhei=c(1,5), hclustfun = function(c){hclust( c, method='ward.D')}, 
				distfun = function (x) as.dist( 1- cor(t(x), method='pearson') ), Rowv=T, ... ) {
			##plot the heatmap as svg image
			if ( nrow(dataObj@data) > nmax ) {
				stop (paste('No plotting for file ',ofile,'- too many genes selected (',nrow(data),')' ))
			}
			if( nrow(dataObj@data) > 2 ){
				data <- as.matrix(t(dataObj@data))
				if ( reorder ){
					data <- data[,order(dataObj@usedObj[['clusters']])]
				}
				brks <- unique(as.vector(c(-20,quantile(data[which(data!= -20)],seq(0,1,by=0.1)),max(data))))
				#rownames( data ) <- paste( dataObj$genes, dataObj$names)
				if ( is.na(hc.row) ){
					hc.row <- hclustfun(distfun(data)) #hclust( as.dist( 1- cor(t(data), method='spearman')), method='ward')
				}
				dendrogram='both'
				if ( length(grep ('color_groups', ofile)) == 0 ) {
					ma <- data[hc.row$order,]
					dendrogram='both'
				}
				else {
					ma <- data
					dendrogram='row'
				}
				if ( ! RowSideColors==F ) {
					ma <- ma[match(geneGroups[order(geneGroups[,3]),1],rownames(ma)),]
					if ( dendrogram=='both'){
						dendrogram='col'
					}else {
						dendrogram='none'
					}
				}
				if (  plotsvg == 1) {
					devSVG( file=paste(ofile,'_Heatmap.svg',sep='') , width=width, height=height)
					if ( ! is.na(ColSideColors) ) {
						if ( RowSideColors != F) {
							heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), key=F, symkey=FALSE,trace='none', 
									cexRow=0.7,cexCol=0.7, main=title,margins = margins, ColSideColors=ColSideColors, RowSideColors=RowSideColors, Rowv=F,dendrogram=dendrogram,lwid = lwid, lhei=lhei, ... )
						}
						else {
							heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), key=F, symkey=FALSE,
									trace='none', cexRow=0.7,cexCol=0.7, main=title,margins = margins, 
									ColSideColors=ColSideColors, hclustfun = hclustfun, distfun = distfun, Rowv=T,dendrogram=dendrogram,lwid = lwid, lhei=lhei, ...)
						}
					}
					else {
						heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), Rowv=F,  key=F, symkey=FALSE,
								trace='none', cexRow=0.7,cexCol=0.7, main=title,margins = margins,
								hclustfun = hclustfun, distfun = distfun, dendrogram=dendrogram,lwid = lwid, lhei=lhei )
					}
					dev.off()
				}
				if ( nrow(data) > 50 ) {
					png( file=paste(ofile,'_Heatmap.png',sep='') , width=width*150, height=height*250 )
				}
				else {
					png( file=paste(ofile,'_Heatmap.png',sep='') , width=width*150, height=height*200 )
				}
				if ( ! is.na(ColSideColors) ) {
					if ( RowSideColors != F) {
						heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), key=F, symkey=FALSE,trace='none', 
								cexRow=2,cexCol=0.7, main=title,margins = margins, ColSideColors=ColSideColors, RowSideColors=RowSideColors, Rowv=F,dendrogram=dendrogram,lwid = lwid, lhei=lhei, ... )
					}
					else {
						heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), key=F, symkey=FALSE,
								trace='none', cexRow=2,cexCol=0.7, main=title,margins = margins, 
								ColSideColors=ColSideColors, hclustfun = hclustfun, distfun = distfun, Rowv=T,dendrogram=dendrogram,lwid = lwid, lhei=lhei, ...)
					}
				}
				else {
					heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), Rowv=F,  key=F, symkey=FALSE,
							trace='none', cexRow=2,cexCol=0.7, main=title,margins = margins,
							hclustfun = hclustfun, distfun = distfun, dendrogram=dendrogram,lwid = lwid, lhei=lhei )
				}
				dev.off()
				write.table( cbind ( 'GeneSymbol' = rownames(ma), 'groupsID' = hc.row$order[hc.row$order], ma),file= paste(ofile,'_data4Genesis.xls', sep=''),sep='\t' )
				write ( rownames(ma),file= paste(ofile,'_Genes_in_order.txt',sep='') ,ncolumns = 1 )
			}
			else {
				print ( paste( 'You have less than two genes for the histogram (',nrow(ma),', ',ofile,') '))
			}
			dataObj@usedObj[['expression.hc.row']] = hc.row
			dataObj
		} 
)

