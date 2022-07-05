#' @name FACS.heatmap
#' @aliases FACS.heatmap,Rscexv-method
#' @rdname FACS.heatmap-methods
#' @docType methods
#' @description Plots the FACS heatmaps
#' @param dataObj  TEXT MISSING
#' @param ofile  TEXT MISSING
#' @param title  TEXT MISSING default='Heatmap'
#' @param nmax  TEXT MISSING default=500
#' @param hc.row  TEXT MISSING default=NA
#' @param ColSideColors  TEXT MISSING default=NA
#' @param RowSideColors  TEXT MISSING default=NA
#' @param FACS.heatmap  TEXT MISSING
#' @param reorder if yes the data matrix will be reordered to the clusters object (default=F)
#' @title description of function FACS.heatmap
#' @export 
setGeneric('FACS.heatmap', ## Name
		function ( dataObj, ofile, title='Heatmap', reorder =F, nmax=500, hc.row=NA, ColSideColors=NA, RowSideColors=NA,width=6, height=6, margins = c(15, 10), hclustfun = function(c){hclust( c, method='ward')}, distfun = function (x) as.dist( 1- cor(t(x), method='pearson') ), ... ){ 
			standardGeneric('FACS.heatmap')
		}
)

setMethod('FACS.heatmap', signature = c ('Rscexv'),
		definition = function ( dataObj, ofile, title='Heatmap', reorder =F, nmax=500, hc.row=NA, ColSideColors=NA, RowSideColors=NA,
				width=6, height=6, margins = c(15, 10), hclustfun = function(c){hclust( c, method='ward.D2')}, distfun = function (x) as.dist( 1- cor(t(x), method='pearson') ), ... ) {
			##plot the heatmap as svg image
			
			if ( ncol(dataObj@facs) > nmax ) {
				stop (paste('No plotting for file ',ofile,'- too many genes selected (',ncol(data),')' ))
			}
			#browser()
			if( ncol(dataObj@facs) > 2 ){
				data <- as.matrix(t(dataObj@facs))
				if ( reorder ){
					data <- data[,order(dataObj@usedObj[['clusters']])]
				}
				if ( length(which(is.na(data))) > 0) {
					v <- rnorm(length(which(is.na(data))))
					v = v + min(v) +1 
					v <- log10(v)
					data[which(is.na(data)) ] = v
				}
				if ( is.na(hc.row) ){
					## this breaks if we have NA values 
					hc.row <- hclustfun(distfun(data)) #hclust( as.dist( 1- cor(t(data), method='spearman')), method='ward')
				}
				ma <- data[hc.row$order,]
				if ( ! is.na(RowSideColors) ) {
					RowSideColors <- RowSideColors[ hc.row$order ]
				} 
				for ( i in 1:2 ){
					#rownames( data ) <- paste( dataObj$genes, dataObj$names)
					if ( i == 1 && plotsvg == 1 ) {
						svglite( file=paste(ofile,'_Heatmap.svg',sep='') , width=width, height=height)
					}
					else {
						png( file=paste(ofile,'_Heatmap.png',sep='') , width=width*150, height=nrow( data ) * 15 + 400 )
					}
					if ( ! is.na(ColSideColors) ) {
						if ( ! is.na(RowSideColors)) {
							heatmap.2(as.matrix(ma), col=bluered, Rowv=F,  key=F, symkey=FALSE,
									trace='none', cexRow=2,cexCol=0.7, main=title,margins = margins, 
									ColSideColors=ColSideColors, RowSideColors=RowSideColors, hclustfun = hclustfun, distfun = distfun,dendrogram='both', ... )
						}
						else {
							heatmap.2(as.matrix(ma), col=bluered, Rowv=T,  key=F, symkey=FALSE,
									trace='none', cexRow=2,cexCol=0.7, main=title,margins = margins, 
									ColSideColors=ColSideColors, hclustfun = hclustfun, distfun = distfun,dendrogram='both',... )
						}
					}
					else {
						heatmap.2(as.matrix(ma), col=bluered, Rowv=F,  key=F, symkey=FALSE,
								trace='none', cexRow=2,cexCol=0.7, main=title,margins = margins,
								hclustfun = hclustfun, distfun = distfun, ... )
					}
					dev.off()
				}
				write.table( cbind ( 'GeneSymbol' = rownames(ma), 'groupsID' = hc.row$order[hc.row$order], ma),file= paste(ofile,'_data4Genesis.txt', sep=''),sep='\t', row.names=F, quote=F  )
				write ( rownames(ma),file= paste(ofile,'_Genes_in_order.txt',sep='') ,ncolumns = 1 )
			}
			else {
				print ( paste( 'You have less than two genes for the histogram (',nrow(ma),', ',ofile,') '))
			}
			
			dataObj@usedObj[['facs.hc.row']] = hc.row
		} 
)

