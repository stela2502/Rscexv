#' @name checkGrouping
#' @aliases checkGrouping,Rscexv-method
#' @rdname checkGrouping-methods
#' @docType methods
#' @description This function should become obsolete.
#' @param userGroups  TEXT MISSING
#' @param data  TEXT MISSING default=NULL
#' @title description of function checkGrouping
#' @export 
setGeneric('checkGrouping', ## Name
		function ( userGroups,  data=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('checkGrouping') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('checkGrouping', signature = c ('data.frame'),
		definition = function ( userGroups,  data=NULL ) {
			userGroups$groupID <- as.vector( as.numeric( userGroups$groupID ))
			if ( !is.null(data) ){
				if ( nrow(data$PCR) != nrow(userGroups) ) {
					### CRAP - rebuilt the grouping information - the data files have been re-created!
					rn <- rownames(data$PCR)
					for ( i in 1:length(rn) ){
						rownames(userGroups) <- userGroups[,1]
						userGroups2 <- as.matrix(userGroups[ rownames(data$PCR), ])
						missing <- which(is.na(userGroups2[,1]))
						userGroups2[missing,1] <- rn[missing]
						userGroups2[missing,2] <- 'previousely dropped'
						userGroups2[missing,3] <- 0
						userGroups2[, 3] <- as.numeric(as.vector(userGroups2[, 3])) +1
						userGroups2 <- as.data.frame(userGroups2)
						userGroups2[,3] <- as.numeric(userGroups2[,3])
						userGroups <- userGroups2
					} 
				}
			}else {
				if ( length(which(userGroups$groupID == 0)) > 0 ){
					userGroups$groupID = userGroups$groupID + 1
				}
				ta <-table(userGroups$groupID)
				exp <- 1:max(as.numeric(userGroups$groupID))
				miss <- exp[(exp %in% names(ta)) == F]
				for ( i in 1:length(miss) ){
					miss[i] = miss[i] -(i -1)
					userGroups$groupID[which(userGroups$groupID > miss[i] )] = userGroups$groupID[which(userGroups$groupID > miss[i] )] -1 
					
				}
			}
			if ( min(userGroups$groupID) > 1 ){
				userGroups$groupID <- userGroups$groupID - ( min(userGroups$groupID) -1 )
			}
			userGroups
		} 
)