#' @name coexpressGenes
#' @aliases coexpressGenes,Rscexv-method
#' @rdname coexpressGenes-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @title description of function coexpressGenes
#' @export 
setGeneric('coexpressGenes', ## Name
	function ( dataObj ) { ## Argumente der generischen Funktion
		standardGeneric('coexpressGenes') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('coexpressGenes', signature = c ('Tool_Coexpression'),
	definition = function ( dataObj ) {
	
	cor.funct <- function ( ma ){
		ma <- ma[, which( apply( ma, 2, function ( x) { length(which( x != 0)) }) > 9 )]
		if ( ncol(ma) < 2 ) {
			NULL
		}
		else {
			cor.t <- cor( ma , method='spearman')
			cor.p <- cor.t
			diag(cor.p) <- 1
			for ( i in 1:(ncol(ma)-1) ) {
				for (a in (i+1):ncol(ma) ) {
					if ( length( as.vector(ma[,i]) ) != length(as.vector(ma[,a]))){
						browser()
					}
					cor.p[a,i] <- cor.test( as.vector(ma[,i]), as.vector(ma[,a]),method='spearman')$p.value
				}
			}
			cor.t.m <- melt(cor.t)
			cor.p.m <- melt(cor.p)
			cor.t.m <- cbind(cor.t.m, cor.p.m[,3])
			cor.t.m <- cor.t.m[which(cor.t.m[,4] < 0.05), ]
			cor.t.m
		}
	}
	ret <- NULL
	for (i in 1:max(dataObj$clusters)){
		t <- cor.funct ( dataObj$PCR[which(dataObj$clusters == i),] )
		if ( ! is.null(t)){
				if ( nrow(t) > 0  ){
					t[,5] <- i
					ret <- rbind(ret, t)
				}
		}
	}
	colnames(ret) <- c('Source.Node','Target.Node', 'rho', 'p.value','Group' )
	ret
} )


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
