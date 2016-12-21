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
			
			write.table( t(data@data), file= file.path(data@outpath,"PCR_data_RawExpression_4_GEO.xls"), sep='\t' )
			
			data.filtered <- sd.filter( data )
						
			plot.histograms( data.filtered ) ## this is needed for the web tool
			
			data.filtered <- norm.PCR(data.filtered,norm.function,max.cyc=max.value, ctrl=ref.genes )
			#plot.heatmap( list( data = t(data.filtered$PCR), genes=colnames(data.filtered$PCR)), 'Contr_filtered_inverted_norm', title='SD filtered inverted data', width=12,height=6,Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},distfun = function(x){ 1- cor(t(x),method='spearman')} )
			write.table( t(data.filtered@data), file=file.path( data.filtered@outpath,"PCR_data_normalized_4_GEO.xls"), sep='\t' )
			
			data.filtered <- z.score.PCR.mad(data.filtered)
			
			write.table( t(data.filtered@data), file=file.path( data.filtered@outpath,"PCR_data_zscored_4_GEO.xls"), sep='\t' )
			
			colnames(data.filtered@annotation) <- c('Gene Name')
			
			data.filtered
		} 
)


