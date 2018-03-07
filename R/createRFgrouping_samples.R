#' @name createRFgrouping_samples
#' @aliases createRFgrouping_samples,Rscexv-method
#' @rdname createRFgrouping_samples-methods
#' @docType methods
#' @description Create a sample grouping data from one RFclust.SGE object
#' @param x the Rscexv object
#' @param RFname the name of the RFclust.SGE object in the Rscexv data. This object has to be populized with data!
#' @param k the number of wanted groups ( default = 10)
#' @param single_res_col the new column in the samples table default= paste('RFgrouping', RFname)
#' @title description of function createRFgrouping_samples
#' @export 
setGeneric('createRFgrouping_samples', ## Name
		function ( x, RFname, k=10, single_res_col = paste('RFgrouping',RFname)) { ## Argumente der generischen Funktion
			standardGeneric('createRFgrouping_samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('createRFgrouping_samples', signature = c ('Rscexv'),
		definition = function ( x, RFname, k=10, single_res_col = paste('RFgrouping',RFname)) {
			if ( is.na( match( RFname, names(x@usedObj[['rfObj']])))){
				stop( paste("the RFname",RFname,"is not defined in this object; defined grouings are:",paste(names(x@usedObj[['rfObj']]), collapse=" ",sep=', ') ) )
			}
			groups <- createGroups( x@usedObj[['rfObj']][[RFname]], k=k, name=RFname )
			x@usedObj[['rfExpressionSets']][[RFname]]@samples <- cbind ( x@usedObj[['rfExpressionSets']][[RFname]]@samples, groups[,3:(2+length(k))] )
			le <- ncol(x@usedObj[['rfExpressionSets']][[RFname]]@samples)
			colnames(x@usedObj[['rfExpressionSets']][[RFname]]@samples)[(le-length(k)+1):le] <- paste('group n=',k)
			m <- max(k)
			## create the predictive random forest object
			x@usedObj[['rfExpressionSets']][[RFname]] <- bestGrouping( x@usedObj[['rfExpressionSets']][[RFname]], group=paste('group n=', m), bestColname = paste('OptimalGrouping',m ,name) )
			x@samples[, paste( single_res_col) ] <-
					predict( x@usedObj[['rfExpressionSets']][[RFname]]@usedObj[[paste( 'predictive RFobj group n=',m) ]], as.matrix(x@data) )
			if ( is.null(x@usedObj[['colorRange']])){
				x@usedObj[['colorRange']] = list()
			}
			x@usedObj[['colorRange']][[paste( single_res_col)]] <- rainbow( length(levels( x@samples[, paste( single_res_col) ])))
			x
		} 
)


