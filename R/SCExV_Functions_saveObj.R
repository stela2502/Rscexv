#' @name saveObj
#' @aliases saveObj,Rscexv-method
#' @rdname saveObj-methods
#' @docType methods
#' @description This function saves the object either as analysis.RData or norm_data.RData if the analysi.RData has not been produced before
#' @param obj the Rscexv object
#' @param file theoutfile default='SCExV_Grps.txt'
#' @title description of function analyse.data
#' @export 
setGeneric('saveObj', ## Name
		function ( data, file='analysis.RData' ){	## Argumente der generischen Funktion
			standardGeneric('saveObj') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('saveObj', signature = c ('Rscexv'),
		definition = function ( data, file='analysis.RData' ){
			exportGroups(data)
			print ( paste('data exported to', file ) )
			save(data , file=file.path(data@outpath, file) )
			if ( locked(file)){
				release.lock( file )
			}
			if ( ! is.null( data@usedObj$usedGrouping )) {
				write( data@usedObj$usedGrouping , file= file.path(data@outpath,'usedGrouping.txt') )
			}
			
		}
)


