#' @name exportGroups
#' @aliases exportGroups,Rscexv-method
#' @rdname exportGroups-methods
#' @docType methods
#' @description This function exports the group names necessary for the SCExV server
#' @param obj the Rscexv object
#' @param file theoutfile default='SCExV_Grps.txt'
#' @title description of function analyse.data
#' @export 
setGeneric('exportGroups', ## Name
		function ( obj, file='SCExV_Grps.txt' ){	## Argumente der generischen Funktion
			standardGeneric('exportGroups') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('exportGroups', signature = c ('Rscexv'),
		definition = function ( obj, file='SCExV_Grps.txt' ){
			v <- NULL
			if ( ! exists('useGrouping') ){
				useGrouping = NULL
			}
			if ( ! is.null(useGrouping)){
				if ( useGrouping == 'ArrayID') {
					v <-'Group by plateID'
				}else {
					v <- useGrouping
				}
			}
			v <- c( v, 'none', 'Group by plateID')
			if ( ncol(obj@samples) > obj@baseSamplesCol ) {
				v <- c( v, colnames(obj@samples)[(obj@baseSamplesCol+1):ncol(obj@samples)] )
			}
			write( v, file= file.path(obj@outpath,file))
			## and the important Sample_Colors.xls file
			write.table( cbind( 
							Samples = obj@samples[,1], 
							ArrayID = obj@samples[,2], 
							Cluster =  obj@usedObj[['clusters']], 
							'color.[rgb]' =  obj@usedObj[['colors']] 
						  ),
					file='Sample_Colors.xls' , row.names=F, sep='\t',quote=F )
			write.table( obj@samples, file='Sample_complete_Data.xls' , row.names=F, sep='\t',quote=F )
		}
)

