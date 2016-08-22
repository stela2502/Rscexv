#' @name reorder.samples
#' @aliases reorder.samples,Rscexv-method
#' @rdname reorder.samples-methods
#' @docType methods
#' @description this function reorderes the Rscexv object based on a column in the samples table (e.g. for plotting)
#' @param dataObj the Rscexv object
#' @param column the samples column to reorder on
#' @title description of function remove.genes
#' @export 
setGeneric('reorder.samples', ## Name
		function ( dataObj, column ) { ## Argumente der generischen Funktion
			standardGeneric('reorder.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('reorder.samples', signature = c ('Rscexv'),
		definition = function ( dataObj, column ) {
			dataObj@data <- dataObj@data[ order( dataObj@samples[,column]),]
			dataObj@samples <- dataObj@samples[order( dataObj@samples[,column]),]
			dataObj
		}
)

