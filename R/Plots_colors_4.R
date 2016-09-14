#' @name colors_4
#' @aliases colors_4,Rscexv-method
#' @rdname colors_4-methods
#' @docType methods
#' @description Create the colour information for a samples or annotation column
#' @param x the Rscexv object
#' @param name the name of data column to colour
#' @param colFunc a colour function like default = function(x) {rainbow(x)}
#' @title description of function createRFgrouping_col
#' @export 
setGeneric('colors_4', ## Name
		function ( x, name,  colFunc = NULL  ) { ## Argumente der generischen Funktion
			standardGeneric('colors_4') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('colors_4', signature = c ('Rscexv'),
		definition = function ( x, name, colFunc = NULL ) {
			if ( is.null(x@usedObj[['colorRange']])){
				x@usedObj[['colorRange']] = list()
			}
			if ( is.null(colFunc) ){
				colFunc = function(x) {rainbow(x)}
			}
			else if ( ! is.na( match( name, names(x@usedObj[['colorRange']]))) ){
				x@usedObj[['colorRange']][[name]] <- NULL
			}
			mkF <- function( x, name, w=1 ) {
				if ( w == 1) {
					if ( is.na(match(class(x@samples[, name ]), 'factor') )) {
						x@samples[, name ] <- factor(x@samples[, name ], levels = unique(x@samples[, name ]) ) 
					}
				}else {
					if ( is.na(match(class(x@annotation[, name ]), 'factor') )) {
						x@annotation[, name ] <- factor(x@annotation[, name ], levels = unique(x@annotation[, name ]) ) 
					}
				}
				 
				x 
			}
			if ( is.na( match( name, names(x@usedObj[['colorRange']])))) {
				if ( !is.na( match(name, colnames(x@samples)))){
					x <- mkF( x, name, 1 )
					x@usedObj[['colorRange']][[name]] <- 
							colFunc( length(levels( x@samples[, name ])))
				}
				else if ( !is.na( match(name, colnames(x@annotation)))){
					x <- mkF( x, name, 2 )
					x@usedObj[['colorRange']][[name]] <- 
							colFunc( length(levels( x@annotation[, name ])))
				}
				else {
					stop( "Sorry this column is nether defined in the samples nor in the annotation table!" )
				}
				
			}
			x
		}
)

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
		function ( x, name=NULL  ) { ## Argumente der generischen Funktion
			standardGeneric('this.color') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
			x@usedObj$colorRange[[name]]
		}
)
