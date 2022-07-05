#' @name quality_of_fit
#' @aliases quality_of_fit,Rscexv-method
#' @rdname quality_of_fit-methods
#' @docType methods
#' @description Calculates a quality of fit
#' @param obj  TEXT MISSING
#' @title description of function quality_of_fit
#' @export 
setGeneric('quality_of_fit', ## Name
		function ( obj ) { 
			standardGeneric('quality_of_fit')
		}
)

setMethod('quality_of_fit', signature = c ('Rscexv'),
		definition = function ( obj ) {
			test <- as.matrix(obj@data)
			rem <- which(test ==  -20 )
			if (length (rem) == 0 ) {
				# that is not possible in a single cell dataset!
				rem <- which(test ==  0 )
			}
			test[which(test ==  -20 ) ] = NA
			ret <- list ( 'per_expression' = apply(test,2, difference, obj ) )
			ret$Expression = round(sum(ret$per_expression))
			if ( obj@wFACS ) {
				test <- obj@facs
				ret$per_FACS = apply(test,2, difference, obj ) 
				ret$FACS = round(sum(ret$per_FACS))
			}
			else {
				ret$per_FACS <- NA
				ret$FACS <- NA
			}
			ret
		} 
)


