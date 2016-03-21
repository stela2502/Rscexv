
#' @name Rscexv
#' @aliases Rscexv,Rscexv-method
#' @rdname Rscexv-methods
#' @docType methods
#' @description create a new Rscexv object from the data files.
#' @param PCR  the pcr data file names default=NULL
#' @param FACS the FACS data file names MISSING default=NULL
#' @param use_pass_fail whether or not to use the pass_fail, information in the PCR data files.
#' @title description of function createDataObj
#' @export 
setGeneric('Rscexv', ## Name
		function  ( PCR=NULL,  FACS=NULL, use_pass_fail = T ){ ## Argumente der generischen Funktion
			standardGeneric('Rscexv') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('Rscexv', signature = c ('character'),
		definition = function ( PCR=NULL,  FACS=NULL, use_pass_fail = T ){
			
			
			PCR <- read.PCR.set ( PCR, use_pass_fail )
		
			system ( 'rm filtered*' )
			
			try ( system ( 'rm R_file_read_error.txt' ) )
			
			if ( is.null(PCR) ){
				system ( 'echo "You need to upload at least one PCR data set in order to start the analysis" > R.error')
				stop("You need to upload at least one PCR data set in order to start the analysis")
			}
			
			colnames(PCR$data) <- str_replace_all( colnames(PCR$data), '/', '_' )
			wFACS=F
			if ( ! is.null(FACS)){
				FACS <- read.FACS.set ( FACS)
				wFACS=T
				## black magick with the FACS gene names
				colnames(FACS$data) <- str_replace_all( colnames(FACS$data), '^P[0-9]+.', '' )
				colnames(FACS$data) <- str_replace_all( colnames(FACS$data), '.Min$', '' )
			}
			
			## check if the samples do overlapp
			
			
			data <- check.dataObj(PCR, FACS)
			
			## now create the object and done
			res <- new('Rscexv', data=data.frame(data$PCR), 
					facs=data.frame(data$FACS), samples=data$samples, 
					annotation=data$annotation, wFACS=wFACS )
			
			res
		}
)



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
#' @param norm.function any of ( )
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
			
			data@PCR <- plug.999(data$PCR, max.value ) ## does nothing for pre-processed data
			
			if ( all ( data@data == 40 ) ){
				system ( 'echo "Please check your filter settings - all samples have been removed from the analysis!" >> R_file_read_error.txt' )
			}
			
			if ( ! is.null(ref.genes)){
				data <- filter.on.controls.no.inv(data,ref.genes,max.ct,max.control)	
			}
			
			## export the unfiltered_not_modified PCR data for publication
			write.table( data@data, file="../PCR_data_4_publication.xls", sep='\t' )
			
			data.filtered <- sd.filter( data )
			
			#plot.heatmap( list( data = t(data.filtered$PCR), genes=colnames(data.filtered$PCR)), 'SD_filtered_not_norm', title='SD filtered RAW data', width=12,height=6,Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},distfun = function(x){ 1- cor(t(x),method='spearman')} )
			
			plot.histograms( data.filtered )
			
			data.filtered$PCR <- norm.PCR(data.filtered$PCR,norm.function,max.cyc=40, ctrl=ref.genes )
			#plot.heatmap( list( data = t(data.filtered$PCR), genes=colnames(data.filtered$PCR)), 'Contr_filtered_inverted_norm', title='SD filtered inverted data', width=12,height=6,Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},distfun = function(x){ 1- cor(t(x),method='spearman')} )
			write.table( data.filtered$z$PCR, file="../PCR_data_normalized_4_publication.xls", sep='\t' )
			
			data.filtered <- z.score.PCR.mad(data.filtered)
			#data.filtered$z$PCR <- data.filtered$PCR
			
			arrays <- arrays <- max(data.filtered$ArrayID)
			cols <- rainbow(arrays)
			cmethod <- 'ward'
			try(PCR.heatmap( 
							list( data = t(data.filtered$PCR), genes=colnames(data.filtered$PCR)),
							'Contr_filtered_inverted_norm', 
							title='SD filtered and normlizied data', 
							width=12,height=6,
							Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},
							distfun = function(x){ as.dist(1- cor(t(x),method='spearman'))}, 
							margins = c(0,15), 
							lwid = c( 0.05, 0.45),
							ColSideColors=cols[data.filtered$ArrayID] )
			)
			
			try(PCR.heatmap( 
							list( data = t(data.filtered$z$PCR), genes=colnames(data.filtered$z$PCR)),
							'Contr_filtered_inverted_norm_Zscore', 
							title='SD filtered normalized data and Z scored', 
							width=12,height=6,
							Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},
							distfun = function(x){as.dist( 1- cor(t(x),method='spearman'))}, 
							margins = c(0,15), 
							lwid = c( 0.05, 0.45),
							ColSideColors=cols[data.filtered$ArrayID] )
			)
			
			devSVG( file="boxplot_filtered_samples.svg", height =8, width =18 )
			par(mar=c(12,5,2,2))
			tmp <- t(data.filtered$PCR)
			tmp [ which(tmp == 0 ) ] <- NA
			boxplot( tmp , las=2,cex.lab=0.5, main ="Normalized expression in all used Samples")
			for ( i in 1:arrays){
				da <- tmp[,which(data.filtered$ArrayID == i ) ]
				abline( h=median(da[which(! is.na(da) ) ]), col = cols[i], lwd=3 )
			}
			dev.off()
			
			devSVG( file="boxplot_filtered_zscored_samples.svg", height =8, width =18 )
			par(mar=c(12,5,2,2))
			tmp <- t(data.filtered$z$PCR)
			tmp [ which(tmp == -20 ) ] <- NA
			boxplot( tmp , las=2,cex.lab=0.5, main="Normalized and z-scored")
			for ( i in 1:arrays){
				da <- tmp[,which(data.filtered$ArrayID == i ) ]
				abline( h=median(da[which(! is.na(da) ) ]), col = cols[i], lwd=3 )
			}
			dev.off()
			
			data.filtered
		} )