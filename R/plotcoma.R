#' @name plotcoma
#' @aliases plotcoma,Rscexv-method
#' @rdname plotcoma-methods
#' @docType methods
#' @description This function plots the correlation plot for the grouping. 
#' @description It will create an svg file if the global variable plotsvg ==1.
#' @param dataObj  TEXT MISSING
#' @param fname  TEXT MISSING default='CorrelationPlot'
#' @title description of function plotcoma
#' @export 
setGeneric('plotcoma', ## Name
		function ( dataObj, fname='/CorrelationPlot' ) { 
			standardGeneric('plotcoma')
		}
)

setMethod('plotcoma', signature = c ('Rscexv'),
		definition = function ( dataObj, fname='/CorrelationPlot' ) {
			if ( is.null(dataObj@usedObj[['coma']] ) ) {
				if ( max(dataObj@usedObj[['clusters']]) > 1){
					dataObj <- coexpressionMatrix ( dataObj )
				}else {
					print ("I can not calculate a co-expression dataset on 1 groups!" )
					return (0)
				}
			}
			png ( file=paste(dataObj@outpath,fname,'.png',sep=''), width=800, height=800 )
			corrplot(dataObj@usedObj[['coma']], order = "hclust", method = "square", hclust.method ='single' )
			dev.off()
			if ( plotsvg == 1 ) {
				svglite ( file=paste(fname,'.svg',sep=''), width=6, height=6 )
				corrplot(dataObj@usedObj[['coma']], order = "hclust", method = "square", hclust.method ='single' )
				dev.off()
			}
		}  
)
