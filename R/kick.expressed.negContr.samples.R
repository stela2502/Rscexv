#' @name kick.expressed.negContr.samples
#' @aliases kick.expressed.negContr.samples,Rscexv-method
#' @rdname kick.expressed.negContr.samples-methods
#' @docType methods
#' @description IF you are sure, that a list of genes must not be expressed in cells of interest you can fore to exclude these ceells using this function.
#' @description But this is a quite bad absolutely not recommended workflow!
#' @description This function drops samples with expression different from 999 - therefore you must call this function before removing the 999 values.
#' @param dataObj the Rscexv object
#' @param negContrGenes which genes should be checked default=NULL
#' @title description of function kick.expressed.negContr.samples
#' @export 
setGeneric('kick.expressed.negContr.samples', ## Name
		function ( dataObj, negContrGenes=NULL ) { 
			standardGeneric('kick.expressed.negContr.samples')
		}
)

setMethod('kick.expressed.negContr.samples', signature = c ('Rscexv'),
		definition = function ( dataObj, negContrGenes=NULL ) {
			if ( ! is.null(negContrGenes) ){
				rem.samples <- NULL
				rem.samples <- unique( 
						apply( dataObj@data[,negContrGenes], 2, function(x) { 
									which(t(dataObj@data[,i]) < 999 )
								} ) 
				)
				dataObj <- remove.samples(dataObj, rem.samples )
			}
			dataObj }
)

