#' @name write.stats
#' @aliases write.stats,Rscexv-method
#' @rdname write.stats-methods
#' @docType methods
#' @description write a statistics table from the lin lang list
#' @param stats the lin lang list default= NULL
#' @param file the outfile default='lin_lang_stats.xls'
#' @title description of function write.stats
setGeneric('write.stats', ## Name
		function ( stats = NULL, file='lin_lang_stats.xls' ) { ## Argumente der generischen Funktion
			standardGeneric('write.stats') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('write.stats', signature = c ('list'),
		definition = function ( stats = NULL, file='lin_lang_stats.xls' ) {
			
			groupL <- function(x) {
				if ( ! is.vector(x$medians)){ x$medians = c(-1,-2) }
				if ( ! is.vector(x$groupIDs)){ x$groupIDs = c(-1,-2) }
				if ( ! is.vector(x$weight)){ x$weight = c(-1,-2) }
				c( x$cor, x$p_value, 
						paste(x$groupIDs[order(x$medians)], collapse =', '), 
						paste(x$medians[order(x$medians)], collapse =', '), 
						paste(x$weight[order(x$medians)], collapse =', ') 
				) }
			ma <- NULL
			if ( ! is.null(stats) ) {
				ma <- t(as.data.frame(lapply(stats, groupL )))
				rownames(ma) <- names(stats)
				colnames(ma)<- c('Correlation', 'p value', 'groups in order', 'median expression in group', 'weight of group' )
				write.table( ma, file=file ,  sep='\t',quote=F ) 
			}
			else {
				print ( "No starts to print!" )
			}
			ma
		} 
)


