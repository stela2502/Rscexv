#' @name checkGrouping
#' @aliases checkGrouping,Rscexv-method
#' @rdname checkGrouping-methods
#' @docType methods
#' @description This function checks that there is no 0 group - probably more later on
#' @param userGroups the grouping vector
#' @param data the optional Rscexv object (at the moment unused)
#' @title description of function checkGrouping
#' @export 
setGeneric('checkGrouping', ## Name
		function ( userGroups,  data=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('checkGrouping') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('checkGrouping', signature = c ('numeric'),
		definition = function ( userGroups,  data=NULL ) {
			## there can not be any sample problems!
			if ( max(userGroups) != 0){ 
				if ( min(userGroups) == 0 ){
					userGroups <- userGroups +1
				}
			}
			userGroups
		} 
)


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
		function (dataObj, gene, ranges) { ## Argumente der generischen Funktion
			standardGeneric('group_1D') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
		function (ma, gene, ranges ) { ## Argumente der generischen Funktion
			standardGeneric('group_1D_worker') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
			checkGrouping ( as.vector(t(userGroups[,3])) )
		} 
)


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
		function ( dataObj, group2sample = list ( '1' = c( 'Sample1', 'Sample2' ) ), name=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('regroup') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
			gr <- checkGrouping ( userGroups[,3] )
			dataObj@samples[,name] = gr
			dataObj
		} 
)

#' @name group_on_strings
#' @aliases group_on_strings,Rscexv-method
#' @rdname group_on_strings-methods
#' @docType methods
#' @description group_an_strings allows the user to group based on a string in the sample names.
#' @param dataObj the Rscexv object
#' @param strings the grouping strings as vector default= c()
#' @param name the new group name
#' @title description of function group_on_strings
#' @export 
setGeneric('group_on_strings', ## Name
		function (dataObj, strings = c(), name='GroupOnStrings' ) { ## Argumente der generischen Funktion
			standardGeneric('group_on_strings') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('group_on_strings', signature = c ('Rscexv'),
		definition = function (dataObj, strings = c(), name='GroupOnStrings' ) {
			userGroups <-  data.frame( cellName = rownames(dataObj@data), userInput = rep.int(0, nrow(dataObj@data)), groupID = rep.int(0, nrow(dataObj@data)) )
			minor = 0	
			for ( i in 1:length(strings) ) {
				g <- grep(strings[i], userGroups$cellName)
				if ( length(g) == 0 ){
					system (paste ('echo "The group name',strings[i] ,'did not match to any sample" > Grouping_R_Error.txt', collaps=" ") )
					minor = minor +1
				}
				else {
					userGroups[g ,3] = i - minor
					userGroups[g ,2] = strings[i]
				}
			}
			gr <- checkGrouping ( userGroups[,3], dataObj )
			dataObj@samples[,name] = gr
			dataObj
		} 
)



