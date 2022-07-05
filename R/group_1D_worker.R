#' @name group_1D_worker
#' @aliases group_1D_worker,Rscexv-method
#' @rdname group_1D_worker-methods
#' @docType methods
#' @description grouping worker for the 1D grouping.
#' @param ma one data matrix from the Rscexv object
#' @param gene the gene to analyze
#' @param ranges the cut offs to create
#' @title description of function group_1D_worker
#' @export 
setGeneric('group_1D_worker', ## Name
		function (ma, gene, ranges ) { 
			standardGeneric('group_1D_worker')
		}
)

setMethod('group_1D_worker', signature = c ('data.frame'),
		definition = function (ma, gene, ranges ) {
			
			position <- which ( colnames(ma) == gene )
			
			userGroups <- data.frame( cellName = rownames(ma), userInput = rep.int(0, nrow(ma)), groupID = rep.int(0, nrow(ma)) )
			
			if ( length(position) > 0 ){
				min <- min(ma[,position])
				max <- max(ma[,position])+1
				
				ranges = ranges[order(ranges)]
				minor = 0
				now <- as.vector( which( ma[,position] >= min & ma[,position] < ranges[1] ))
				userGroups$userInput[now] = paste ('min <= x <',ranges[1] )
				userGroups$groupID[now] = 1
				for ( i in 2:length(ranges) ) {
					now <- as.vector( which( ma[,position] >= ranges[i-1] & ma[,position] < ranges[i] ))
					userGroups$userInput[now] = paste(ranges[i-1],'<= x <',ranges[i])
					if ( length(now) > 0 ){
						userGroups$groupID[now] = i
					}
					else {
						minor = minor + 1
					}
				}
				now <- as.vector( which( ma[,position] >= ranges[length(ranges)] & ma[,position] < max ))
				userGroups$userInput[now] = paste(ranges[length(ranges)],'<= x < max')
				userGroups$groupID[now] = length(ranges) +1
			}
			checkGrouping ( userGroups )
		} 
)


