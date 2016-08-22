#' @name checkGrouping
#' @aliases checkGrouping,Rscexv-method
#' @rdname checkGrouping-methods
#' @docType methods
#' @description This function checks that there is no 0 group - probably more later on
#' @param userGroups the grouping vector
#' @param data the optional Rscexv object (at the moment unused)
#' @title description of function checkGrouping
#' @export 
setGeneric('checkGrouping', ## Name
		function ( userGroups,  data=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('checkGrouping') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('checkGrouping', signature = c ('data.frame'),
		definition = function ( userGroups,  data=NULL ) {
			## there can not be any sample problems!
			if ( max(userGroups[,3]) != 0){ 
				if ( min(userGroups[,3]) == 0 ){
					userGroups[,3] <- userGroups[,3] +1
				}
			}
			as.vector(t(userGroups[,3]))
		} 
)


