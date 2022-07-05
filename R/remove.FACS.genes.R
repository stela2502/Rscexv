#' @name remove.FACS.genes
#' @aliases remove.FACS.genes,Rscexv-method
#' @rdname remove.FACS.genes-methods
#' @docType methods
#' @description this function removes a list of FACS genes from the analysis
#' @param dataObj the Rscexv object
#' @param ids the gene ids to remove
#' @title description of function remove.FACS.genes
#' @export 
setGeneric('remove.FACS.genes', ## Name
		function ( dataObj, ids ) { 
			standardGeneric('remove.FACS.genes')
		}
)

setMethod('remove.FACS.genes', signature = c ('Rscexv'),
		definition = function ( dataObj, ids ) {
			if ( length(ids) > 0 ){
				write ( colnames(dataObj@facs)[ids], file="./filtered_genes.txt",ncolumn=1, append=T )
				dataObj@facs <- dataObj@facs[,-ids]
			}
			else {
				print ( "No genes to filter out!" )
			}
			dataObj	
		} 
)

