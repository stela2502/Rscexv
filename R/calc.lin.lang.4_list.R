#' @name calc.lin.lang.4_list
#' @aliases calc.lin.lang.4_list,Rscexv-method
#' @rdname calc.lin.lang.4_list-methods
#' @docType methods
#' @description This function calculates the correlation between the median expression value of a group and the fraction of expression cells in the group.
#' @param l the list of values
#' @title description of function calc.lin.lang.4_list
#' @export 
setGeneric('calc.lin.lang.4_list', ## Name
		function ( l =list() ) { ## Argumente der generischen Funktion
			standardGeneric('calc.lin.lang.4_list') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('calc.lin.lang.4_list', signature = c ('list'),
		definition = function ( l =list() ) {
			n = length(l)
			or <- NULL
			medians <- NULL
			groupIDs<- NULL
			weight <- NULL
			for (i in 1:n) {
				data <- l[[i]][ is.na(l[[i]]) ==F ]
				if ( length(data) == 0) {
					data <- c(0)
				}
				data.min <- min(data)
				data.max <- max(data)
				if ( data.min != data.max ) {
					or <- c( or, length(data) / length(l[[i]] ))
					weight <- c( weight,  length(l[[i]] ) )
					medians <- c( medians,  median(data) )
					groupIDs <- c( groupIDs, i )
				}
			}
			ret <- list (  'weight' = weight, medians = medians, groupIDs = groupIDs, or = or , cor = NULL)
			if( length(weight) > 2 ){
				ret$cor <- corr( cbind(or, medians), w = weight / sum(weight) )
			}
			
			ret
		} 
)


