#' @name filter.on.controls.no.inv
#' @aliases filter.on.controls.no.inv,Rscexv-method
#' @rdname filter.on.controls.no.inv-methods
#' @docType methods
#' @description This function filteres out samples that have an expression of less than thresh in there ref.nms(gene names)
#' @param dataObj the Rscexv object
#' @param ref.nms the gene names used as positive control
#' @param thresh the minimum expression in each ref gene
#' @param howm drop if one sample has howm genes that fail the test
#' @title description of function filter.on.controls.no.inv
#' @export 
setGeneric('filter.on.controls.no.inv', ## Name
		function (dataObj,ref.nms,thresh,howm) { ## Argumente der generischen Funktion
			standardGeneric('filter.on.controls.no.inv') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('filter.on.controls.no.inv', signature = c ('Rscexv'),
		definition = function (dataObj,ref.nms,thresh,howm) {
			
			conts <- as.matrix(dataObj@data[,ref.nms ])
			conts[which(conts < thresh)] <- 0
			conts[conts>0] <- 1
			
			if ( is.null(dim(conts)) ){
				nums <- as.vector(conts)
			}
			else {
				nums <- apply(conts,1,sum)
			}
			
			dataObj <- remove.samples(dataObj, which(nums>howm) )
			dataObj
		} 
)

