#' @name force.unique.sample
#' @aliases force.unique.sample,Rscexv-method
#' @rdname force.unique.sample-methods
#' @docType methods
#' @description This function will just make sure, that all samples have uniaue sample names, even if they come from different arrays.
#' @param x  TEXT MISSING
#' @title description of function force.unique.sample
#' @export 
setGeneric('force.unique.sample', ## Name
	function ( x ) { 
		standardGeneric('force.unique.sample')
	}
)

setMethod('force.unique.sample', signature = c ('character'),
	definition = function ( x ) {
	ret <- NULL
	last <- ""
	use <- NULL
	id <- 1
	repl <- vector(length=length(x))
	for ( i in 1:length(x) ){
		if ( is.null(ret) ){
			last = x[i]
			ret <- c( last )
			use  <- last
		}
		else if ( x[i] != last ){
			last = x[i]
			if ( ! is.na(match( last, ret )) ){
				use <- paste(last,"_",id, sep = '')
				ret <- c(ret , use )
				id = id + 1
			}else {
				use  <- last
				ret <- c(ret,  last )
			}
		}
		repl[i] <- use
	}
	l <- list( 'unique' = ret, 'replacement' = repl )
	l
} )

