#' @name this.color
#' @aliases this.color,Rscexv-method
#' @rdname this.color-methods
#' @docType methods
#' @description Get the actual color range from the object
#' @param x the Rscexv object
#' @param name the name of data column or empty if you want the last auto color
#' @title description of function createRFgrouping_col
#' @export 
setGeneric('this.color', ## Name
		function ( x, name=NULL  ) { 
			standardGeneric('this.color')
		}
)

setMethod('this.color', signature = c ('Rscexv'),
		definition = function ( x, name=NULL ) {
			if ( is.null(name) ){
				name = paste( 'auto_clusters', x@usedObj[['auto_clusters']] ,sep='.')
			}
			if (is.null( x@usedObj$colorRange[[name]] )) {
				x <- colors_4( x, name )	
			}
			x@usedObj$usedGrouping <- name
			x
		}
)