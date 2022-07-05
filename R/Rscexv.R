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
		function  ( PCR=NULL,  FACS=NULL, use_pass_fail = T ){ 
			standardGeneric('Rscexv')
		}
)

setMethod('Rscexv', signature = c ('character'),
		definition = function ( PCR=NULL,  FACS=NULL, use_pass_fail = T ){
			
			fnamesPCR <- PCR
			fnamesFACS <- FACS
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
			
			data$samples$fnamesPCR <- fnamesPCR[data$samples$ArrayID]
			if ( wFACS ) {
				data$samples$fnamesFACS <- fnamesFACS[data$samples$ArrayID]
			}
			missing <- setdiff( 1:length(fnamesPCR), unique(data$samples$ArrayID) )
			if ( length(missing) > 0 ){
				stop( paste( length(missing),"files failed to load into the analysis:", paste( collapse= ', ', fnamesPCR[missing] ),". Most likely due to gene name missmatches." ))
			}
			## now create the object and done
			res <- new('Rscexv', data=data.frame(data$PCR), 
					facs=data.frame(data$FACS), samples=data$samples, 
					annotation=data$annotation, wFACS=wFACS, outpath=pwd(), baseSamplesCol=ncol(data$samples) )
			colnames(res@annotation) <- c('Gene Name' )
			res
		}
)

