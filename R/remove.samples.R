#' @name remove.samples
#' @aliases remove.samples,Rscexv-method
#' @rdname remove.samples-methods
#' @docType methods
#' @description Remove samples by id. 
#' @param dataObj the Rscexv object
#' @param ids which samples (ids!) to remove
#' @title description of function remove.samples
#' @export 
setGeneric('remove.samples', ## Name
		function ( dataObj, ids ) { ## Argumente der generischen Funktion
			standardGeneric('remove.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('remove.samples', signature = c ('Rscexv'),
		definition = function ( dataObj, ids ) {
			
			if ( length(ids) > 0 ){
				write ( rownames(dataObj@data)[ids], file="./filtered_samples.txt",ncolumn=1, append=T )
				dataObj@data <- dataObj@data[-ids,]
				if ( dataObj@wFACS ){
					dataObj@facs <- dataObj@facs[-ids,]
				}
				if ( ncol(dataObj@snorm) > 0 ){
					dataObj@snorm <- dataObj@snorm[-ids,]
				}
				if ( ncol(dataObj@raw) > 0 ){
					dataObj@raw <- dataObj@raw[-ids,]
				}
				dataObj@samples <- dataObj@samples[-ids,]
			}
			else {
				print ( "No samples to filter out!" )
			}
			dataObj	
		} 
)

