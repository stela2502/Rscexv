#' @name implyCloseOrder
#' @aliases implyCloseOrder,Rscexv-method
#' @rdname implyCloseOrder-methods
#' @docType methods
#' @description This function reorderes a grouping using hclust. The grouing will simply be converted into a factor.
#' @param dataObj the Rscexv data object
#' @param groupName the name of the grouping that has to be reordered
#' @title description of function implyCloseOrder
#' @export 
setGeneric('implyCloseOrder', ## Name
		function ( dataObj,  groupName=NULL, hclust.method= 'complete' ) { ## Argumente der generischen Funktion
			standardGeneric('implyCloseOrder') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('implyCloseOrder', signature = c ('Rscexv'),
		definition = function ( dataObj,  groupName=NULL, hclust.method= 'ward.D2' ) {
			
			dataObj@usedObj[['clusters']] <- as.numeric( dataObj@samples[,groupName] )
			dataObj <- coexpressionMatrix ( dataObj )
			hc <- hclust ( as.dist( 1- dataObj@usedObj[['coma']]), method= hclust.method )
			if ( is.factor(dataObj@samples[,groupName]) ) {
				interm <-factor( dataObj@samples[,groupName], levels=levels(data@samples[, 5])[hc$order] )
				dataObj@samples[,groupName] <- factor( as.numeric(interm) )
			}else {
				interm <- factor( dataObj@samples[,groupName], levels=	hc$order )
				dataObj@samples[,groupName] <- factor( as.numeric(interm) )
			}
			dataObj@usedObj$clusters <- as.numeric( dataObj@samples[,groupName] )
			dataObj@usedObj$coma <- NULL
			dataObj
		} 
)

