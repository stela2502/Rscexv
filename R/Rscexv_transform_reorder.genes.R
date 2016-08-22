#' @name reorder.genes
#' @aliases reorder.genes,Rscexv-method
#' @rdname reorder.genes-methods
#' @docType methods
#' @description this function reorderes the Rscexv object based on a column in the annotation table (e.g. for plotting)
#' @param dataObj the Rscexv object
#' @param column the annotation column to reorder on
#' @title description of function remove.genes
#' @export 
setGeneric('reorder.genes', ## Name
		function ( dataObj, column ) { ## Argumente der generischen Funktion
			standardGeneric('reorder.genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('reorder.genes', signature = c ('Rscexv'),
		definition = function ( dataObj, column ) {
			dataObj@data <- dataObj@data[ , order( dataObj@annotation[,column])]
			dataObj@annotation <- dataObj@annotation[order( dataObj@annotation[,column]),]
			dataObj
		}
)

