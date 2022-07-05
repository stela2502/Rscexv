#' @name regroup
#' @aliases regroup,Rscexv-method
#' @rdname regroup-methods
#' @docType methods
#' @description This function is used to re-group the data from within SCExV - there should be easier ways in R.
#' @param dataObj  TEXT MISSING
#' @param group2sample  TEXT MISSING default= list ( '1' = c( 'Sample1'
#' @param groupName which group to re-group
#' @title description of function regroup
#' @export 
setGeneric('regroup', ## Name
		function ( dataObj, group2sample = list ( '1' = c( 'Sample1', 'Sample2' ) ), name=NULL ) { 
			standardGeneric('regroup')
		}
)

setMethod('regroup', signature = c ('Rscexv'),
		definition = function ( dataObj, group2sample = list ( '1' = c( 'Sample1', 'Sample2' ) ), name=NULL ) {
			userGroups <-  data.frame( cellName = rownames(dataObj@data), userInput = rep.int(0, nrow(dataObj@data)), groupID = rep.int(0, nrow(dataObj@data)) )
			n <- names(group2sample)
			n <- n[order( n )]
			minor = 0
			for ( i in 1:length(n) ){
				if ( sum(is.na(match(group2sample[[i]], userGroups$cellName))==F) == 0 ){
					minor = minor +1
				}
				else if (  sum(is.na(match(group2sample[[i]], userGroups$cellName))==F) < length( group2sample[[i]] ) ) {
					rows <- match(group2sample[[i]], userGroups$cellName)
					rows <- rows[-which(is.na(rows))]
					userGroups[ rows ,3] <- i
				}
				else { ## the whole group has failed??
					userGroups[ match(group2sample[[i]], userGroups$cellName) ,3] = i - minor
				}
			}
			gr <- checkGrouping ( userGroups )
			dataObj@samples[,name] = gr
			dataObj
		} 
)

