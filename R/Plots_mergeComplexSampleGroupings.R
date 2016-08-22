#' @name mergeComplexSampleGroupings
#' @aliases mergeComplexSampleGroupings,Rscexv-method
#' @rdname mergeComplexSampleGroupings-methods
#' @docType methods
#' @description Create a new group and create shaded rainbow color shemata stored in the x@usedObj[['colorRange']] variable
#' @param x the Rscexv object
#' @param g1 the first grouping name (can be complex)
#' @param g2 the second grouping name (can be complex too)
#' @param newName the name of the new variable to create
#' @title description of function mergeSampleGroupings
#' @export 
setGeneric('mergeComplexSampleGroupings', ## Name
		function ( x,  g1, g2, newName ) { ## Argumente der generischen Funktion
			standardGeneric('mergeComplexSampleGroupings') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('mergeComplexSampleGroupings', signature = c ('Rscexv'),
		definition = function ( x, g1, g2, newName ){
			if(is.null( x@usedObj[['colorRange']] ) ){
				x@usedObj[['colorRange']] <- list()
			}
			newG <- factor(paste( as.vector(x@samples[, g1]), as.vector( x@samples[, g2]) , sep='_' ))
			
			if ( is.na(match( newName, colnames(x@samples))) ){
				x@samples =cbind( x@samples, newName=newG )
				colnames(x@samples)[ncol(x@samples)] = newName
			}else{
				x@samples[,newName] = newG
			}
			if ( is.na( match( newName, names(x@usedObj[['colorRange']]) ) ) ) {
				x@usedObj[['colorRange']][[length(x@usedObj[['colorRange']])+1]] <- 1
				names(x@usedObj[['colorRange']])[[length(x@usedObj[['colorRange']])]] = newName
			}
			x@usedObj[['colorRange']][[ match( newName, names(x@usedObj[['colorRange']])) ]] <-  rainbow( max(as.numeric(newG)) )
			x
		}
)

