#' @name group_1D
#' @aliases group_1D,Rscexv-method
#' @rdname group_1D-methods
#' @docType methods
#' @description Create grouping based on the expression of one gene
#' @param dataObj the Rscexv object
#' @param gene the gene the groups are based on
#' @param ranges the ranges for the rgoups
#' @title description of function group_1D
#' @export 
setGeneric('group_1D', ## Name
		function (dataObj, gene, ranges) { 
			standardGeneric('group_1D')
		}
)

setMethod('group_1D', signature = c ('Rscexv'),
		definition = function (dataObj, gene, ranges) {
			name= paste(gene, '1D Group')
			userGroups <- group_1D_worker ( dataObj@data, gene, ranges)
			if ( max(userGroups) == 0 ){
				userGroups <- group_1D_worker ( dataObj@facs, gene, ranges)
			}
			if ( is.na(match( name, colnames(dataObj@samples))) ){
				dataObj@samples = cbind(dataObj@samples, userGroups)
				colnames(dataObj@samples)[ncol(dataObj@samples)] = name
			}else{
				dataObj@samples[,name] = userGroups
			}
			dataObj
		} 
)
