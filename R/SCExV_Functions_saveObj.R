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
			if ( file.exists( file.path(data@outpath, file) ) ){
				print ( 'data exported to analysis.RData')
				save(data , file=file.path(data@outpath, file) )
			}else {
				data.filtered <- data
				save(data , file= file.path(data@outpath, 'norm_data.RData') )
			}
		}
)


