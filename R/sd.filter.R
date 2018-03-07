#' @name sd.filter
#' @aliases sd.filter,Rscexv-method
#' @rdname sd.filter-methods
#' @docType methods
#' @description This function is removing cells and genes that do not show any variation in the data.
#' @param dataObj the Rscexv object
#' @title description of function sd.filter
#' @export 
setGeneric('sd.filter', ## Name
		function (dataObj) { ## Argumente der generischen Funktion
			standardGeneric('sd.filter') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('sd.filter', signature = c ('Rscexv'),
		definition = function (dataObj) {
			
			sds1 <- NULL
			sds2 <- NULL

			sds1 <- apply(dataObj@data,1,sd)
			sds2 <- apply(dataObj@data,2,sd)
			
			cuto1 <- which(sds1==0)
			cuto2 <- which(sds2==0)
			
			if(length(cuto1) >0){
				write ( colnames(dataObj@data)[cuto2], file="./filtered_samples.txt",ncolumn=1, append=T )
				dataObj <- remove.samples(dataObj, cuto1 )
			}
			
			if(length(cuto2) >0){
				write ( colnames(dataObj@data)[cuto2], file="./filtered_genes.txt",ncolumn=1, append=T )
				dataObj <- remove.genes(dataObj, cuto2 )
			}
			
			dataObj			
		} 
)

