


#' @name reorder_on_correlation
#' @aliases reorder_on_correlation,Rscexv-method
#' @rdname reorder_on_correlation-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param order  TEXT MISSING default='hclust'
#' @param hclust.method  TEXT MISSING default= 'single' 
#' @param ...  TEXT MISSING
#' @title description of function reorder_on_correlation
#' @export 
setGeneric('reorder_on_correlation', ## Name
	function ( dataObj, order='hclust', hclust.method= 'single' ,... ) { ## Argumente der generischen Funktion
		standardGeneric('reorder_on_correlation') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('reorder_on_correlation', signature = c ('Tool_Coexpression'),
	definition = function ( dataObj, order='hclust', hclust.method= 'single' ,... ) {
	if ( exists ( 'oldclusters', where=dataObj)) {
		return ( dataObj )
	}
	if ( ! exists ( 'coma', where =dataObj ) ) {
		dataObj <- coexpressionMatrix ( dataObj )
	}
	newOrder <- corrMatOrder(dataObj$coma, order=order, hclust.method=  hclust.method, ...)
	#rg <- vector('list', max(newOrder))
	#names(rg) <- 1:max(newOrder)
	dataObj$oldclusters <- dataObj$clusters
	for ( i in 1:max(dataObj$clusters) ) {
		#rg[[newOrder[i]]] <- as.vector(rownames(dataObj$PCR)[ which( dataObj$clusters == i )] )
		dataObj$clusters[ which( dataObj$oldclusters == i )] <- newOrder[i]
	}
	#regroup( dataObj, group2sample= rg )
	dataObj
} )
