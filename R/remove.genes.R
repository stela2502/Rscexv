#' @name remove.genes
#' @aliases remove.genes,Rscexv-method
#' @rdname remove.genes-methods
#' @docType methods
#' @description this function removes a list of PCR genes from the analysis
#' @param dataObj the Rscexv object
#' @param ids the gene ids to remove
#' @title description of function remove.genes
#' @export 
setGeneric('remove.genes', ## Name
		function ( dataObj, ids ) { ## Argumente der generischen Funktion
			standardGeneric('remove.genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('remove.genes', signature = c ('Rscexv'),
		definition = function ( dataObj, ids ) {
			if ( length(ids) > 0 ){
				write ( colnames(dataObj@data)[ids], file="./filtered_genes.txt",ncolumn=1, append=T )
				dataObj@data <- dataObj@data[,-ids]
				dataObj@annotation <- data.frame(dataObj@annotation[-ids,])
				if ( ncol(dataObj@snorm) > 0 ){
					dataObj@snorm <- dataObj@snorm[,-ids]
				}
				if ( ncol(dataObj@raw) > 0 ){
					dataObj@raw <- dataObj@raw[,-ids]
				}
			}
			else {
				print ( "No genes to filter out!" )
			}
			dataObj	
		} 
)

