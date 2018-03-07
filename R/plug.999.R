#' @name plug.999
#' @aliases plug.999,Rscexv-method
#' @rdname plug.999-methods
#' @docType methods
#' @description this function replaces the 999 (failed) values by the plug value.
#' @param x the Rscexv object
#' @param plug the replacement value (default=45)
#' @title description of function plug.999
#' @export 
setGeneric('plug.999', ## Name
		function (x,plug=45) { ## Argumente der generischen Funktion
			standardGeneric('plug.999') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plug.999', signature = c ('Rscexv'),
		definition = function (x,plug=45) {
			ma <- as.matrix(x@data)
			ind <- which(ma==999)
			ma[ind] <- plug
			x@data <- data.frame(ma)
			x
		} 
)

