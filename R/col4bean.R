#' @name col4bean
#' @aliases col4bean,Rscexv-method
#' @rdname col4bean-methods
#' @docType methods
#' @description get the color information for the beans
#' @param x  TEXT MISSING
#' @param tic  TEXT MISSING default='black'
#' @title description of function col4bean
#' @export 
setGeneric('col4bean', ## Name
		function (x, tic='black') { ## Argumente der generischen Funktion
			standardGeneric('col4bean') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('col4bean', signature = c ('character'),
		definition = function (x, tic='black') {
			ret <- list()
			for ( i in 1:length(x) ){
				ret[[i]] <- c(x[i], tic )
			}
			ret
		}
)


