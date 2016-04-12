#' @name createDataObj
#' @aliases createDataObj,Rscexv-method
#' @rdname createDataObj-methods
#' @docType methods
#' @description This function is the main function used by SCExV. In a script the steps should be executed one by one.
#' @param PCR  the pcr data file names default=NULL
#' @param FACS the FACS data file names MISSING default=NULL
#' @param max.value an optional maximum value for all failed wells default=40
#' @param ref.genes if the data should be normalized to reference genes put them into this list
#' @param max.ct all sammples with a ct value of more than max.ct in the ref.genes is dropped
#' @param max.control samples are dropped if 1 (0), 2(1), ... etc ref.genes shows a ct value of more than max.ct
#' @param norm.function any of ( "none","mean control genes","max expression","median expression","quantile" )
#' @param negContrGenes which negative control genes to use (samples with expression there are dropped!)
#' @param use_pass_fail whether or not to use the pass_fail, information in the PCR data files.
#' @title description of function createDataObj
#' @export 
setGeneric('createDataObj', ## Name
		function ( PCR=NULL,  FACS=NULL, max.value=40, ref.genes=NULL, max.ct=25, 
				max.control=0,  norm.function='none', negContrGenes=NULL, use_pass_fail = T, ...){ ## Argumente der generischen Funktion
			standardGeneric('createDataObj') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)
setMethod('createDataObj', signature = c ('character'),
		definition = function ( PCR=NULL,  FACS=NULL, max.value=40,	ref.genes=NULL, max.ct=25, max.control=0,  norm.function='none', negContrGenes=NULL, use_pass_fail = T, ... ){
			
			data <- Rscexv( PCR, FACS, use_pass_fail)
			
			data <- kick.expressed.negContr.samples(data, negContrGenes )
			
			data <- plug.999(data, max.value ) ## does nothing for pre-processed data
			
			if ( all ( data@data == 40 ) ){
				system ( 'echo "Please check your filter settings - all samples have been removed from the analysis!" >> R_file_read_error.txt' )
			}
			
			if ( ! is.null(ref.genes)){
				data <- filter.on.controls.no.inv(data,ref.genes,max.ct,max.control)	
			}
			
			## export the unfiltered_not_modified PCR data for publication
			
			write.table( data@data, file= file.path(data@outpath,"InbuiltData.RData"), sep='\t' )
			
			data.filtered <- sd.filter( data )
						
			plot.histograms( data.filtered ) ## this is needed for the web tool
			
			data.filtered <- norm.PCR(data.filtered,norm.function,max.cyc=max.value, ctrl=ref.genes )
			#plot.heatmap( list( data = t(data.filtered$PCR), genes=colnames(data.filtered$PCR)), 'Contr_filtered_inverted_norm', title='SD filtered inverted data', width=12,height=6,Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},distfun = function(x){ 1- cor(t(x),method='spearman')} )
			write.table( data.filtered@data, file=file.path( data.filtered@outpath,"PCR_data_normalized_4_publication.xls"), sep='\t' )
			
			data.filtered <- z.score.PCR.mad(data.filtered)
			
			
			data.filtered
		} 
)


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
#' @title description of function analyse.data
#' @export 
setGeneric('analyse.data', ## Name
		function (obj,onwhat='Expression',groups.n, cmethod, clusterby='MDS', 
				ctype='hierarchical clust', beanplots=TRUE, move.neg = FALSE, plot.neg=TRUE, useGrouping=NULL, plotsvg = 0, ...){	## Argumente der generischen Funktion
			standardGeneric('analyse.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('analyse.data', signature = c ('Rscexv'),
		definition = function (obj,onwhat='Expression',groups.n, cmethod, clusterby='MDS', 
				ctype='hierarchical clust', beanplots=TRUE, move.neg = FALSE, plot.neg=TRUE, useGrouping=NULL, plotsvg = 0, ...){
			
			cols = rainbow( groups.n )
			
			if ( ! obj@wFACS ) {
				onwhat="Expression"
			} else if ( ncol(obj@facs)< 4 ) {
				onwhat="Expression"
			}
			obj <- mds.and.clus(obj,onwhat= onwhat,groups.n = groups.n, cmethod=cmethod, clusterby=clusterby,ctype=ctype, useGrouping=useGrouping)
			
			try(plotcoma(obj) )
			if ( length(which(obj@usedObj[['clusters']] == 0)) > 0 ){
				obj@usedObj[['clusters']] <- obj@usedObj[['clusters']] + 1
			}
			obj@usedObj[['colors']] <- apply( t(col2rgb( cols ) ), 1, paste,collapse=' ')[obj@usedObj[['clusters']]]
			
			## plot the mds data
			try(plotDR( obj@usedObj[['mds.proj']][order(obj@usedObj[['clusters']]),], col=cols, labels=obj@usedObj[['clusters']][order(obj@usedObj[['clusters']])] ),silent=F)
			try(writeWebGL( width=470, height=470, dir = file.path(obj@outpath,'webGL')),silent=F)
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
			
			if ( exists ('geneGroups') ){
				geneGroups$groupID = as.vector(geneGroups$groupID)
				if ( is.vector(as.vector(geneGroups$groupID)) ) {
					t <- obj
					obj$z$PCR <- obj$z$PCR[, order(geneGroups$groupID)]
					RowV = FALSE
					RowSideColors=c(gray.colors(max(geneGroups$groupID),start = 0.3, end = 0.9))[as.numeric(geneGroups$groupID[order(geneGroups$groupID)])]
				}
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
							hclustfun = function(c){hclust( c, method=cmethod)}
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
							hclustfun = function(c){hclust( c, method=cmethod)}
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
			plot.funct( obj, groups.n, clus =  obj@usedObj[['clusters']], boot = 1000, plot.neg=plot.neg, mv=mv )
			
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

#' @name exportGroups
#' @aliases exportGroups,Rscexv-method
#' @rdname exportGroups-methods
#' @docType methods
#' @description This function exports the group names necessary for the SCExV server
#' @param obj the Rscexv object
#' @param file theoutfile default='SCExV_Grps.txt'
#' @title description of function analyse.data
#' @export 
setGeneric('exportGroups', ## Name
		function ( obj, file='SCExV_Grps.txt' ){	## Argumente der generischen Funktion
			standardGeneric('exportGroups') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('exportGroups', signature = c ('Rscexv'),
		definition = function ( obj, file='SCExV_Grps.txt' ){
			v <- NULL
			if ( ! exists('useGrouping') ){
				useGrouping = NULL
			}
			if ( ! is.null(useGrouping)){
				if ( useGrouping == 'ArrayID') {
					v <-'Group by plateID'
				}else {
					v <- useGrouping
				}
			}
			v <- c( v, 'none', 'Group by plateID')
			if ( ncol(obj@samples) > obj@baseSamplesCol ) {
				v <- c( v, colnames(obj@samples)[(obj@baseSamplesCol+1):ncol(obj@samples)] )
			}
			write( v, file= file.path(obj@outpath,file))
			## and the important Sample_Colors.xls file
			write.table( cbind( 
							Samples = obj@samples[,1], 
							ArrayID = obj@samples[,2], 
							Cluster =  obj@usedObj[['clusters']], 
							'color.[rgb]' =  obj@usedObj[['colors']] 
						  ),
					file='Sample_Colors.xls' , row.names=F, sep='\t',quote=F )
			write.table( obj@samples, file='Sample_complete_Data.xls' , row.names=F, sep='\t',quote=F )
		}
)

#' @name saveObj
#' @aliases saveObj,Rscexv-method
#' @rdname saveObj-methods
#' @docType methods
#' @description This function saves the object either as analysis.RData or norm_data.RData if the analysi.RData has not been produced before
#' @param obj the Rscexv object
#' @param file theoutfile default='SCExV_Grps.txt'
#' @title description of function analyse.data
#' @export 
setGeneric('saveObj', ## Name
		function ( data, file='analysis.RData' ){	## Argumente der generischen Funktion
			standardGeneric('saveObj') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('saveObj', signature = c ('Rscexv'),
		definition = function ( data, file='analysis.RData' ){
			exportGroups(data)
			if ( file.exists( file.path(data@outpath, file) ) ){
				print ( 'data exported to analysis.RData')
				save(data , file=file.path(data@outpath, file) )
			}else {
				data.filtered <- data
				save(data , file= file.path(data@outpath, 'norm_data.RData') )
			}
		}
)


#' @name plotDensity
#' @aliases plotDensity,Rscexv-method
#' @rdname plotDensity-methods
#' @docType methods
#' @description This function creates the density plot for the analysis page
#' @param obj the Rscexv object
#' @title description of function analyse.data
#' @export 
setGeneric('plotDensity', ## Name
		function ( obj ){	## Argumente der generischen Funktion
			standardGeneric('plotDensity') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plotDensity', signature = c ('Rscexv'),
		definition = function ( obj ){
			usable <- is.na(match( obj@usedObj[['clusters']], which(table(as.factor(obj@usedObj[['clusters']])) < 4 ) )) == T
			
			use <- obj
			use@usedObj[['clusters']] <- obj@usedObj[['clusters']][usable]
			use@usedObj[['mds.proj']] <- obj@usedObj[['mds.proj']][usable,]
			cols <- rainbow(max(as.numeric(obj@usedObj[['clusters']])))
			H <- Hkda( use@usedObj[['mds.proj']], use@usedObj[['clusters']], bw='plugin')
			kda.fhat <- kda( use@usedObj[['mds.proj']], use@usedObj[['clusters']],Hs=H, compute.cont=TRUE)
			try(plot(kda.fhat, size=0.001, colors = cols[as.numeric(names(table(use@usedObj[['clusters']])))] ),silent=F)
			try( writeWebGL(
							dir = file.path(obj@outpath,'densityWebGL'), 
							width=470, height=470, prefix='K', 
							template= system.file(file.path("densityWebGL.html"), 
									package = "Rscexv") 
							),
					silent=F 
			)
		}
)
