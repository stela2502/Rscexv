#' @name difference
#' @aliases difference,Rscexv-method
#' @rdname difference-methods
#' @docType methods
#' @description This function calculates the 
#' @param x  TEXT MISSING
#' @param obj  TEXT MISSING
#' @title description of function difference
#' @export 
setGeneric('difference', ## Name
		function ( x, obj ) { 
			standardGeneric('difference')
		}
)

setMethod('difference', signature = c ('numeric'),
		definition = function ( x, obj ) {
			ret = 0 
			for ( i in 1:groups.n  ) {
				a <- x[which( obj@usedObj[['clusters']] == i)]
				a <- a[- (is.na(a))==F]
				if ( length(a) > 1 ) {  ret = ret + sum( (a- mean(a) )^2 ) }
			}
			ret
		} 
)