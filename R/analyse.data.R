#' @name analyse.data
#' @aliases analyse.data,Rscexv-method
#' @rdname analyse.data-methods
#' @docType methods
#' @description This function implements the backend functionallity of the SCExV analyse page.
#' @param obj the Rscexv object
#' @param onwhat cluster on weither Expression of FACS data default='Expression'
#' @param groups.n the number of groups to create
#' @param cmethod clustering method
#' @param clusterby which dataset to cluter on (raw or MDS) default='MDS'
#' @param ctype the clzutering type default = 'hierarchical clust'
#' @param beanplots create beanplots or violinplots
#' @param move.neg should the not expressed value be as set to the minimum expression value?
#' @param plot.neg show the not expressing samples in the bean/violinplots?
#' @param useGrouping do not calculate a new grouping - use this column in the samples table (default=NULL)
#' @param plotsvg create svg figures in addition to the png ones (default=0)
#' @param mds.type which mds algorithm to use (PCA, LLE, ISOMAP or ZIFA) default = PCA
#' @param geneGroups the annotation gene group to use in the heatmaps (default = NULL)
#' @title description of function analyse.data
#' @export 
setGeneric('analyse.data', ## Name
		function (obj,onwhat='Expression',groups.n, cmethod, clusterby='MDS', 
				ctype='hierarchical clust', beanplots=TRUE, move.neg = FALSE, plot.neg=TRUE, useGrouping=NULL, plotsvg = 0, mds.type='PCA', geneGroups=NULL, ...){	## Argumente der generischen Funktion
			standardGeneric('analyse.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('analyse.data', signature = c ('Rscexv'),
		definition = function (obj,onwhat='Expression',groups.n, cmethod, clusterby='MDS', 
				ctype='hierarchical clust', beanplots=TRUE, move.neg = FALSE, plot.neg=TRUE, useGrouping=NULL, plotsvg = 0, mds.type='PCA', geneGroups=NULL, ...){
			
			
			
			if ( ! obj@wFACS ) {
				onwhat="Expression"
			} else if ( ncol(obj@facs)< 4 ) {
				onwhat="Expression"
			}
			
			if ( !is.null(useGrouping) ) {
				if ( useGrouping == 'none' ){
					useGrouping = NULL
					obj@usedObj$usedGrouping = NULL
				}
			}else {
				obj@usedObj$usedGrouping = NULL
			}
			
			obj <- mds.and.clus(obj,onwhat= onwhat,groups.n = groups.n, cmethod=cmethod, 
					clusterby=clusterby,ctype=ctype, useGrouping=useGrouping,mds.type=mds.type)
			
			useGrouping = obj@usedObj$usedGrouping
			cols = obj@usedObj$colorRange[[useGrouping]]
			
			try(plotcoma(obj) )
			
			if ( length(which(obj@usedObj[['clusters']] == 0)) > 0 ){
				obj@usedObj[['clusters']] <- obj@usedObj[['clusters']] + 1
			}
			obj@usedObj[['colors']] <- apply( t(col2rgb( cols ) ), 1, paste,collapse=' ')[obj@usedObj[['clusters']]]
			silent= F
			## plot the mds data
			try(plotDR( obj@usedObj[['mds.proj']][order(obj@usedObj[['clusters']]),], col=cols, labels=obj@usedObj[['clusters']][order(obj@usedObj[['clusters']])] ),silent=silent)
			try(writeWebGL( width=470, height=470, dir = file.path(obj@outpath,'webGL')),silent=silent)
			png(file=file.path(obj@outpath,'webGL', 'MDS_2D.png'), width=800,height=800)
			plotDR( obj@usedObj[['mds.proj']][order(obj@usedObj[['clusters']]),1:2], col=cols, labels=obj@usedObj[['clusters']][order(obj@usedObj[['clusters']])] )
			dev.off()
			#save( obj, file='clusters.RData')
			write.table (obj@usedObj[['mds.proj']][order(obj@usedObj[['clusters']]),1:2], file = file.path(obj@outpath,'2D_data.xls') )
			sample.cols.rgb <-t(col2rgb( cols[obj@usedObj[['clusters']][order(obj@usedObj[['clusters']])]]))
			sample.cols.rgb <- cbind(sample.cols.rgb,  colorname = cols[obj@usedObj[['clusters']][order(obj@usedObj[['clusters']])]] )
			rownames(sample.cols.rgb) <- rownames(obj@data)[order(obj@usedObj[['clusters']])]
			write.table ( sample.cols.rgb , file = file.path(obj@outpath,'2D_data_color.xls') )
			write.table (cbind( names = cols, t(col2rgb( cols))), file=file.path(obj@outpath,'webGL','MDS_2D.png.cols'),sep='\\t',  row.names=F,quote=F )
			
			obj@usedObj[['quality_of_fit']] = quality_of_fit(obj)

			RowV = TRUE
			RowSideColors = FALSE
			
			if ( ! is.null(geneGroups) ){
				geneGroup <- as.numeric(obj@annotation[,geneGroups])
				t <- obj
				obj@data <- obj@data[, order(geneGroup) ]
				RowV = FALSE
				RowSideColors=c(gray.colors(max(geneGroup),start = 0.3, end = 0.9))[geneGroup[order(geneGroup)] ]
			}
			## plot the heatmaps
			
			try( PCR.heatmap ( obj , 
							file.path(obj@outpath,'PCR_color_groups'), 
							title='PCR data', 
							ColSideColors=cols[obj@usedObj[['clusters']]][order(obj@usedObj[['clusters']])],
							RowSideColors=RowSideColors,
							reorder=T,
							width=12,
							height=6, 
							margins = c(1,11), 
							lwid = c( 1,6), lhei=c(1,5),
							Rowv=RowV,
							Colv=F,
							hclustfun = function(c){hclust( c, method=cmethod)}
					), silent=F)
			
#	try( collapsed_heatmaps (obj, what='PCR', functions = c('median', 'mean', 'var', 'quantile70' )), silent=T)
#	try( collapsed_heatmaps (obj, what='FACS', functions = c('median', 'mean', 'var', 'quantile70' )), silent=T)
			try( PCR.heatmap ( obj, 
							file.path(obj@outpath,'PCR'), 
							title='PCR data', 
							ColSideColors=cols[obj@usedObj[['clusters']]],
							width=12,
							height=6, 
							margins = c(1,11), 
							lwid = c( 1,6), lhei=c(1,5),
							hclustfun = function(c){hclust( c, method=cmethod)}
					), silent=F)
			try( FACS.heatmap ( obj, 
							file.path(obj@outpath,'facs'), 
							title='FACS data', 
							ColSideColors=cols[obj@usedObj[['clusters']]],
							width=12,
							height=6, 
							hc.col= obj@usedObj[['hc']],
							margins = c(1,11), 
							lwid = c( 1,6), lhei=c(1,5),
						#	hclustfun = function(c){hclust( c, method=cmethod)}
					), silent=T)
			
			try( FACS.heatmap ( obj, 
							file.path(obj@outpath,'facs_color_groups'), 
							reorder=T,
							title='FACS data', 
							ColSideColors=cols[obj@usedObj[['clusters']]][order(obj@usedObj[['clusters']])],
							width=12,
							height=6, 
							hc.col= obj@usedObj[['hc']],
							margins = c(1,11), 
							lwid = c( 1,6), lhei=c(1,5),
							Colv=F,
						#	hclustfun = function(c){hclust( c, method=cmethod)}
					), silent=T)
			
			ma  <- NULL
			mv <- NULL
			if ( zscoredVioplot == 1 ){
				ma <- as.matrix(obj@data)
				mv <- -20
			}else {
				ma <- as.matrix(obj@snorm)
				mv <- 0
			}
			
			if ( move.neg ){
				neg <- which(ma == mv )
				m <- min(ma[-neg])
				mv <- m -1
				ma[neg] = mv
			}
			if ( beanplots ) {
				plot.funct <-  function (x , ... ) { plot.beans( x, ...) }
			}else{
				plot.funct <-  function (x, ...) { plot.violines( x, ...) }
			}
			
			obj@usedObj[['for.plot']] = ma
			plot.funct( obj, groups.n, clus =  obj@usedObj[['clusters']], plot.neg=plot.neg, mv=mv )
			
			#print ( paste( 'plot.funct( ma , groups.n, clus =  obj@usedObj[["clusters"]], boot = 1000, plot.neg =',plot.neg,', mv =', mv))
			if ( obj@wFACS ){
				write.table( cbind( Samples = rownames(obj@facs), obj@facs ), file=file.path(obj@outpath,'merged_FACS_Table.xls') , row.names=F, sep='	',quote=F )
				all.data <- cbind(obj@data, obj@facs )
				write.table(cbind( Samples = rownames(all.data), all.data ), file=file.path(obj@outpath,'merged_data_Table.xls') , row.names=F, sep='	',quote=F )
			}else {
				write.table( cbind( Samples = rownames(obj@data), obj@data ), file=file.path(obj@outpath,'merged_data_Table.xls') , row.names=F, sep='	',quote=F )
			}
						
			## the lists in one file
			
			write.table( obj@samples,file=file.path(obj@outpath,'Sample_Colors.xls') , row.names=F, sep='	',quote=F )
			
			## now I need to write a groungs info file with all possible groupings to choose from
			exportGroups( obj )
			obj
		} 
)

