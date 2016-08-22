#' @name mergeSampleGroupings
#' @aliases mergeSampleGroupings,Rscexv-method
#' @rdname mergeSampleGroupings-methods
#' @docType methods
#' @description Create a new group and create shaded rainbow color shemata stored in the x@usedObj[['colorRange']] variable
#' @param x the Rscexv object
#' @param g1 the first grouping name (can be complex)
#' @param g2 the second grouping name (has to be simple yes/no!)
#' @param newName the name of the new variable to create
#' @title description of function mergeSampleGroupings
#' @export 
setGeneric('mergeSampleGroupings', ## Name
		function ( x,  g1, g2, newName ) { ## Argumente der generischen Funktion
			standardGeneric('mergeSampleGroupings') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('mergeSampleGroupings', signature = c ('Rscexv'),
		definition = function ( x, g1, g2, newName ){
	if (length( table( x@samples[,g2] ) ) > 2) {
		stop( "This function can only merge one complex grouping with a yes/no information")
	}
	if ( ! all.equal( names( table(x@samples[,g2] )), c('no','yes')) ){
		stop( "This function can only merge one complex grouping with a yes/no information")
	}
	if(is.null( x@usedObj[['colorRange']] ) ){
		x@usedObj[['colorRange']] <- list()
	}
	interleave <- function(v1,v2)
	{
		ord1 <- 2*(1:length(v1))-1
		ord2 <- 2*(1:length(v2))
		c(v1,v2)[order(c(ord1,ord2))]
	}
	newG <- (as.numeric(t(x@samples[, g1])) *2)-1
	newG[ which( x@samples[,g2] == 'yes') ] = newG[ which( x@samples[,g2] == 'yes') ] +1
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
	x@usedObj[['colorRange']][[ match( newName, names(x@usedObj[['colorRange']])) ]] <- interleave ( rainbow( max(as.numeric(x@samples[,g1]))), rainbow( max(as.numeric(x@samples[,g1])), s=0.5, v=0.5) )
	x
	}
)

