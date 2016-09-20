#' @name exportGroups
#' @aliases exportGroups,Rscexv-method
#' @rdname exportGroups-methods
#' @docType methods
#' @description This function exports the group names necessary for the SCExV server
#' @param obj the Rscexv object
#' @param sample.file the outfile for the sample groups; default='SCExV_Grps.txt'
#' @param gene.file the outfile for the gene groups; default ='GeneGroupings.txt'
#' @title description of function analyse.data
#' @export 
setGeneric('exportGroups', ## Name
		function ( obj, sample.file='SCExV_Grps.txt', gene.file='GeneGroupings.txt' ){	## Argumente der generischen Funktion
			standardGeneric('exportGroups') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('exportGroups', signature = c ('Rscexv'),
		definition = function ( obj, sample.file='SCExV_Grps.txt',  gene.file='GeneGroupings.txt'){
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
			write( v, file= file.path(obj@outpath,sample.file))
			if ( ! is.null(obj@usedObj[['clusters']])) {
			## the gene grouping
			g <- c( 'none' )
			if ( ncol(obj@annotation) > 1 ){
				g <- c( g, colnames(obj@annotation)[2:ncol(obj@annotation)])
			}
			write( g, file= file.path(obj@outpath,gene.file))
			
			## and the important Sample_Colors.xls file
			write.table( cbind( obj@samples[,c(1,2)],
							grouping =  obj@usedObj[['clusters']], 
							'colors' =  obj@usedObj[['colors']] 
						  ),
					file='Sample_Colors.xls' , row.names=F, sep='\t',quote=F )
			write.table( obj@samples, file='Sample_complete_Data.xls' , row.names=F, sep='\t',quote=F )
		}
			
		}
)



