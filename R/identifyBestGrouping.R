#' @name identifyBestGrouping
#' @aliases identifyBestGrouping,Rscexv-method
#' @rdname identifyBestGrouping-methods
#' @docType methods
#' @description This function compares several groupings to each other. 
#' @description In addition the groupings are combined in a way, that if a sample end up in group 1 in grouing A but in group 2 in groupung B
#' @description it might be something else than a sample ending up in group 1 for both.
#' @description The quality of the groupings will be accessed in two stated (1) there should not be too many small groups and
#' @description (2) the genes should be differentially expressed in these groups. The differential expression is accessed 
#' @description using a straight forward anova approach (excluding the not called genes for a Rscexv object).
#' @param x the Rscexv
#' @param groups the colnames (samples) that contain the grouping information (default)
#' @param namePrefix a common name prefix for this analysis
#' @param cut define the p value cut off for the test ( if too view signfican genes are detected at 0.05)
#' @return A list containing the number of significant genes in each possible combination of groups and the Rscexv object containg all annotations and stats.
#' @title description of function identifyBestGrouping
#' @export 
setGeneric('identifyBestGrouping', ## Name
		function (x, groups, namePrefix='identifyBestGrouping', cut=0.05) { 
			standardGeneric('identifyBestGrouping')
		}
)

setMethod('identifyBestGrouping', signature = c ('Rscexv'),
		definition = function (x, groups, namePrefix='identifyBestGrouping', cut=0.05) {
			## I do not want to loose all so get me the most out of this gouping
			ret <- list( 'x' = 0, groups = 0 )
			groupLengthT <- function ( a, g ) {
				sum (which(t) < 10 ) < 20 && length(t) < 20
			}
			groupPaste <- function ( a, g, name ) {
				if ( ! is.na(match (name,colnames(a@samples)) ) ) {
					stop( paste( "The column",name,'already exists - STOP') )
				}
				a@samples[, name ] <- apply( a@samples[, g ],1, function (x) {paste(x,collapse= ' ') } )
				a
			}

			## first oder the group cols by the output of identifyBestGrouping			
			te <- qualityTest ( x, groups, numbers=T , cut=cut)
			x <- te$x
			groups <- groups[ order( te$res, decreasing =T ) ]
			names = c(paste( namePrefix, 'All (',length(groups),')' ))
			x <- groupPaste (x, groups, names)
			for ( i in 2:length(groups) ) {
				## then paste them best to worst together and check where you get better stats
				x <- groupPaste( x, groups[1:i], paste( namePrefix, i,'/',length(groups) ) )
				names<- c(names, paste( namePrefix, i,'/',length(groups) ))
			}
			ret <- qualityTest ( x, names, cut=cut )
			ret$names <- te$res
			ret
} )






