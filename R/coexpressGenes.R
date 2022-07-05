#' @name coexpressGenes
#' @aliases coexpressGenes,Rscexv-method
#' @rdname coexpressGenes-methods
#' @docType methods
#' @description  calculates the coexpression for all genes in all groups in the data set
#' @param dataObj the Rscexv object
#' @param grouping the column in the samples table that describes the grouping to use
#' @param pcutoff the minimum p value to report the co-expression for
#' @param file an optional filename to store the correlation to. This file can be read by Cytoscape (default = NULL)
#' @param keepCrap sometimes you want to collect each and every correlation value - even the useless ones (default = 0)
#' @title description of function coexpressGenes
#' @return A table that contains correlation for EVERY gene gene combination in all groups + all data.
#' @return falied correlations will have a NA value but still a p value of 1.
#' @export 
setGeneric('coexpressGenes', ## Name
		function ( dataObj, grouping=NULL, pcutoff= 0.05, file=NULL, keepCrap=0  ) { 
			standardGeneric('coexpressGenes')
		}
)

setMethod('coexpressGenes', signature = c ('Rscexv'),
		definition = function ( dataObj, grouping=NULL, pcutoff= 0.05, file=NULL, keepCrap=0 ) {
			if ( keepCrap ==1 ){
				pcutoff = 1
			}
			cor.funct <- function ( ma ){
				#ma <- ma[, which( apply( ma, 2, function ( x) { length(which( x != 0)) }) > 9 )]
				if ( ncol(ma) < 2 ) {
					NULL
				}
				else {
					cor.t <- cor( ma , method='spearman')
					diag(cor.t) <- 1
					cor.p <- cor.t
					cor.p[] <- 1 
					diag(cor.p) <- 0
					for ( i in 1:(ncol(ma)-1) ) {
						for (a in (i+1):ncol(ma) ) {
							cor.p[a,i] <- cor.test( as.vector(ma[,i]), as.vector(ma[,a]),method='spearman')$p.value
						}
					}
					cor.t.m <- melt(cor.t)
					cor.p.m <- melt(cor.p)
					cor.t.m <- cbind(cor.t.m, cor.p.m[,3])
					if ( keepCrap==1) {
						cor.t.m[which(is.na(cor.t.m[,4])),4] <- 1
					}
					cor.t.m <- cor.t.m[which(cor.t.m[,4] <= pcutoff), ]
					cor.t.m
				}
			}
			ret <- NULL
			for (i in unique(as.vector(dataObj@samples[,grouping] ))){
				t <- cor.funct ( dataObj@snorm[which(dataObj@samples[,grouping]== i),] )
				if ( ! is.null(t)){
					if ( nrow(t) > 0  ){
						t[,5] <- i
						ret <- rbind(ret, t)
					}
				}
			}
			i <- 0
			t <- cor.funct ( dataObj@snorm )
			if ( ! is.null(t)){
				if ( nrow(t) > 0  ){
					t[,5] <- i
					ret <- rbind(ret, t)
				}
			}
			colnames(ret) <- c('Source.Node','Target.Node', 'rho', 'p.value','Group' )
			if ( ! is.null(file) ) {
				write.table( ret, file , sep=" ",quote=F, row.names=F )
			}
			ret
		} 
)

