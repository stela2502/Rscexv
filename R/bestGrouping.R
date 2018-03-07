#' @name bestGrouping
#' @aliases 'bestGrouping,Rscexv-method
#' @docType methods
#' @description The function is using a randomForest classifier with 2000 trees to classify the given data using the given grooping
#' @description All groups that fail to be prediceted using the random forest are deemed ungrouped.
#' @description All groups where less than 50 percent of the total samples geting classified as being from that group fail.
#' @param x the single cells ngs object
#' @param group a vector of sample columns that should be checked (the most complex is used only)
#' @param bestColname the column name to store the best grouping in
#' @param cutoff the cutoff percentage where all groups showing less than this percentacge of remapped samples are dropped
#' @title description of function randomForest
#' @return a distRF object to be analyzed by pamNew
#' @export 
setGeneric('bestGrouping',
		function ( x, group , bestColname='QualifiedGrouping', cutoff=0.5){
			standardGeneric('bestGrouping')
		}
)
setMethod('bestGrouping', signature = c ('Rscexv'),
		definition = function (x, group, bestColname='QualifiedGrouping' , cutoff=0.5) {
			uObj <- paste( 'predictive RFobj', group )
			rf <- NULL
			if (  is.null( x@usedObj[[uObj]])){
				x@usedObj[[uObj]] <- randomForest( x= as.matrix(x@data), y=factor(x@samples[, group]),ntree=2000 )
			}
			if ( FALSE){
			t <- table( observed= x@samples[,group ], predicted = x@usedObj[[uObj]]$predicted )
			i <- 0
			r <- vector('numeric', ncol(t))
			names(r) <- colnames(t)
			for (i in 1:nrow(t)) {
				if ( which(t[i,] == max(t[i,])) == i) {
					r[i]= max(t[i,]) / sum(t[i,])
				}
				else {
					r[i]= 0
				}
			}
			BAD <- which(r < cutoff )
			## remove an optional 'gr. ' from the group ids
			x@samples[,bestColname] <- as.numeric(str_replace_all( x@samples[, group], 'gr. ', ''))
			for ( b in BAD ) {
				x@samples[ which(x@samples[,bestColname] == b), bestColname] <- 0
			}
			
			for (i in 0:(length(table(x@samples[,bestColname]))-1)){
				modify <- which(x@samples[,bestColname] >= i )
				if ( length(modify) == 0 ) { break}
				while( length(which(x@samples[modify,bestColname] == i)) == 0 ){
					x@samples[modify,bestColname] = x@samples[modify, bestColname] -1
				}
			}
			x@samples[,bestColname] <- paste( 'gr.', x@samples[,bestColname])
		}
			x
		}
)

