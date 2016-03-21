
#' @name plug.999
#' @aliases plug.999,Rscexv-method
#' @rdname plug.999-methods
#' @docType methods
#' @description this function replaces the 999 (failed) values by the plug value.
#' @param tab the data.frame containge the read PCR table
#' @param plug the replacement value (default=45)
#' @title description of function plug.999
#' @export 
setGeneric('plug.999', ## Name
		function (tab,plug=45) { ## Argumente der generischen Funktion
			standardGeneric('plug.999') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plug.999', signature = c ('data.frame'),
		definition = function (tab,plug=45) {
			
			ind <- which(tab==999)
			tab[ind] <- plug
			tab
			
		} 
)

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



