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
			write.table( data@data, file= paste( data@outpath, "/PCR_data_4_publication.xls",sep=''), sep='\t' )
			
			data.filtered <- sd.filter( data )
						
			plot.histograms( data.filtered ) ## this is needed for the web tool
			
			data.filtered <- norm.PCR(data.filtered,norm.function,max.cyc=max.value, ctrl=ref.genes )
			#plot.heatmap( list( data = t(data.filtered$PCR), genes=colnames(data.filtered$PCR)), 'Contr_filtered_inverted_norm', title='SD filtered inverted data', width=12,height=6,Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},distfun = function(x){ 1- cor(t(x),method='spearman')} )
			write.table( data.filtered@data, file=paste( data@outpath,"/PCR_data_normalized_4_publication.xls",sep=''), sep='\t' )
			
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
#' @title description of function analyse.data
#' @export 
setGeneric('analyse.data', ## Name
		function (obj,onwhat='Expression',groups.n, cmethod, clusterby='MDS', 
				ctype='hierarchical clust', beanplots=TRUE, move.neg = FALSE, plot.neg=TRUE, useGrouping=NULL, ...){	## Argumente der generischen Funktion
			standardGeneric('analyse.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('analyse.data', signature = c ('Rscexv'),
		definition = function (obj,onwhat='Expression',groups.n, cmethod, clusterby='MDS', 
				ctype='hierarchical clust', beanplots=TRUE, move.neg = FALSE, plot.neg=TRUE, ...){
			
		#	if ( exists( 'userGroups' )) {
		#		userGroups <- checkGrouping(userGroups, obj)
		#	}
			cols = rainbow( groups.n )
			
			if ( is.null(obj@facs)) {
				onwhat="Expression"
			} else if ( ncol(obj@facs)< 4 ) {
				onwhat="Expression"
			}
			obj <- mds.and.clus(obj,onwhat= onwhat,groups.n = groups.n, cmethod=cmethod, clusterby=clusterby,ctype=ctype)
			
			plotcoma(obj)
			if ( length(which(obj@usedObj[['clusters']] == 0)) > 0 ){
				obj@usedObj[['clusters']] <- obj@usedObj[['clusters']] + 1
			}
			obj@usedObj[['colors']] <- apply( t(col2rgb( cols ) ), 1, paste,collapse=' ')[obj@usedObj[['clusters']]]
			
			## plot the mds data
			try(plotDR( obj@usedObj[['mds.proj']][order(obj@usedObj[['clusters']]),], col=cols, labels=obj@usedObj[['clusters']][order(obj@usedObj[['clusters']])] ),silent=F)
			
			try(writeWebGL( width=470, height=470 ),silent=F)
			png(file='./webGL/MDS_2D.png', width=800,height=800)
			plotDR( obj@usedObj[['mds.proj']][order(obj@usedObj[['clusters']]),1:2], col=cols, labels=obj@usedObj[['clusters']][order(obj@usedObj[['clusters']])] )
			dev.off()
			#save( obj, file='clusters.RData')
			write.table (obj@usedObj[['mds.proj']][order(obj@usedObj[['clusters']]),1:2], file = './2D_data.xls' )
			sample.cols.rgb <-t(col2rgb( cols[obj@usedObj[['clusters']][order(obj@usedObj[['clusters']])]]))
			sample.cols.rgb <- cbind(sample.cols.rgb,  colorname = cols[obj@usedObj[['clusters']][order(obj@usedObj[['clusters']])]] )
			rownames(sample.cols.rgb) <- rownames(obj@data)[order(obj@usedObj[['clusters']])]
			write.table ( sample.cols.rgb , file = './2D_data_color.xls' )
			write.table (cbind( names = cols, t(col2rgb( cols))), file='./webGL/MDS_2D.png.cols', sep='\\t',  row.names=F,quote=F )
			
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
							'./PCR_color_groups', 
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
					), silent=T)
			
#	try( collapsed_heatmaps (obj, what='PCR', functions = c('median', 'mean', 'var', 'quantile70' )), silent=T)
#	try( collapsed_heatmaps (obj, what='FACS', functions = c('median', 'mean', 'var', 'quantile70' )), silent=T)
			try( PCR.heatmap ( obj, 
							'./PCR', 
							title='PCR data', 
							ColSideColors=cols[obj@usedObj[['clusters']]],
							width=12,
							height=6, 
							margins = c(1,11), 
							lwid = c( 1,6), lhei=c(1,5),
							hclustfun = function(c){hclust( c, method=cmethod)}
					), silent=T)
			try( FACS.heatmap ( obj, 
							'./facs', 
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
							'./facs_color_groups', 
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
				plot.funct <-  plot.beans
			}else{
				plot.funct <-  plot.violines
			}
			
			## plot the violoines
			if ( ! is.null(obj$FACS)){
				plot.funct( obj$FACS, groups.n, clus =  obj@usedObj[['clusters']], boot = 1000, plot.neg=plot.neg, mv=mv )
			}
			print ( paste( 'plot.funct( ma , groups.n, clus =  obj@usedObj[["clusters"]], boot = 1000, plot.neg =',plot.neg,', mv =', mv))
			plot.funct( ma , groups.n, clus =  obj@usedObj[['clusters']], boot = 1000, plot.neg=plot.neg, mv = mv  )
			
			obj
		} )





