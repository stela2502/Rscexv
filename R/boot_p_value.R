#' @name boot_p_value
#' @aliases boot_p_value,Rscexv-method
#' @rdname boot_p_value-methods
#' @docType methods
#' @description calculates the p value from the boot strap distribution for the lin_lang analysis
#' @param cmp the boot strap list of correlation values
#' @param real_val the correlation value
#' @param i the numer of bootstrap runs
#' @title description of function boot_p_value
#' @export 
setGeneric('boot_p_value', ## Name
		function (cmp, real_val, i ) { 
			standardGeneric('boot_p_value')
		}
)

setMethod('boot_p_value', signature = c ('numeric'),
		definition = function (cmp, real_val, i ) {
			a <- length(which(cmp > real_val))
			if ( a == 0 ) { a<-1}
			p_value = a / i
			p_value
		} 
)

