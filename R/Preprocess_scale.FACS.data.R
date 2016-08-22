#' @name scale.FACS.data
#' @aliases scale.FACS.data,Rscexv-method
#' @rdname scale.FACS.data-methods
#' @docType methods
#' @description This function applies a simple log10 onto all expression data.
#' @param dataObj  TEXT MISSING
#' @title description of function scale.FACS.data
#' @export 
setGeneric('scale.FACS.data', ## Name
		function (dataObj) { ## Argumente der generischen Funktion
			standardGeneric('scale.FACS.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('scale.FACS.data', signature = c ('data.frame'),
		definition = function (dataObj) {
			rown <- rownames( dataObj )
			dataObj <-  apply( dataObj,2, as.numeric)			
			dataObj[which(dataObj < 1)] <- 1
			dataObj <- log10(dataObj)
			rownames(dataObj) <- rown
			dataObj
		} 
)

