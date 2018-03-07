#' @name qualityTest
#' @aliases qualityTest,Rscexv-method
#' @rdname qualityTest-methods
#' @docType methods
#' @description This function calculates an anova test for the groupCol and all data.
#' @return A list containing the updated 'x' Rscexv object and res - a vector of pasted gene names that turned out to be significantly different.
#' @param x The Rscexv object
#' @param groups the group name to clister the data (default=NULL)
#' @param cut the p value cut off (BenjaminHochberg corrected; default=0.05)
#' @param numbers if true return amount of significant genes, not the names
#' @title description of function qualityTest
#' @export 
setGeneric('qualityTest', ## Name
		function (x, groups=NULL, cut=0.05, numbers=F  ) { ## Argumente der generischen Funktion
			standardGeneric('qualityTest') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('qualityTest', signature = c ('Rscexv'),
		definition = function (x, groups=NULL, cut=0.05, numbers=F ) {
			if (is.null(groups) ) {
				stop( "please give me a vector of columns to check!" )
			}
			res <- vector( 'numeric', length(groups))
			names(res) <- groups
			for ( i in 1:length(groups) ) {
				x <- simpleAnova( x,groups[i] )
				grname <- paste('simpleAnova', groups[i])
				g <-which(x@stats[[grname]][,3] < cut)
				
				if ( numbers ){
					res[i] <- length( x@stats[[grname]][which(x@stats[[grname]][,3]< cut),1] )
				}
				else {
					if ( length( g ) > 0 ) {
						res[i] <- paste(collapse=" ",as.vector(x@stats[[grname]][which(x@stats[[grname]][,3] < cut),1] ))
					}
					else{
						res[i] <- ''
					}
				}
			}
			list ( res = res, x=x)
		} )

