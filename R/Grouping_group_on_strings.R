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
			gr <- checkGrouping ( userGroups, dataObj )
			dataObj@samples[,name] = gr
			dataObj
		} 
)



