#' @name plot.beans
#' @aliases plot.beans,Rscexv-method
#' @rdname plot.beans-methods
#' @docType methods
#' @description plot the beanplots
#' @param x the Rscexv object
#' @param groups.n  TEXT MISSING
#' @param clus  TEXT MISSING
#' @param plot.neg  TEXT MISSING default=TRUE
#' @param mv  TEXT MISSING default=-20
#' @param subpath the subpath for the plots (default = '')
#' @param names a vector of grou names (default NULL)
#' @param col a colour vector for colouring the beans (default=NULL -> rainbow )
#' @title description of function plot.beans
#' @export 
setGeneric('plot.beans', ## Name
		function ( x, groups.n, clus,  plot.neg=TRUE, mv=-20, subpath='', names=NULL, col=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('plot.beans') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plot.beans', signature = c ('Rscexv'),
		definition = function ( x, groups.n, clus, plot.neg=TRUE, mv=-20, subpath='', names=NULL, col=NULL ) {
			ma <- NULL
			if ( x@wFACS ){
				ma <- t(cbind( x@usedObj[['for.plot']], x@facs ))
			}else{
				ma <- t(x@usedObj[['for.plot']])
			}
			n <- rownames(ma)
			if ( is.null(col)){
				col = this.color(x,useGrouping)
			}
			cols = col4bean(col)	
			s <-  split(seq(ncol(ma)), clus)
			if ( subpath != '' ){
				system( paste('mkdir',file.path(x@outpath,subpath ) ) )
			}
			for ( i in 1:nrow( ma ) ) {
				#print (paste( 'plot.beans working on gene', n[i] ) )
				fn <- file.path(x@outpath,subpath, n[i] )
				png( file=paste(fn,'.png',sep=''), width=800,height=800)
				lila <- vector('list', groups.n)
				lila$names <- NULL
				for( a in 1:groups.n){
					lila[[a]]=ma[i,which(clus == a)]
					lila$names <- c( lila$names, paste(length(which(lila[[a]] != mv)), length(lila[[a]]) ,sep='/' ) )
					if ( ! plot.neg ){
						lila[[a]][which(lila[[a]] == mv)] <- NA
						if ( sum(is.na(lila[[a]]) ) == length( lila[[a]]) ){
							lila[[a]] = c(0)
						}
					}
				}
				if ( ! is.null(names) && length(names) == length(lila$names) ){
					lila$names= paste( names, lila$names)
				}
				lila$main <- n[i]
				lila$what <- c(1,1,0,1) ## not plot medain line
				lila$col <- cols
				try(do.call(beanplot,lila), silent=F )
				dev.off()
				if ( plotsvg == 1 ) {
					devSVG( file=paste(fn,'.svg',sep=''), width=6,height=6)
					#lila$cex.axis=0.5
					try(do.call(beanplot,lila), silent=T )
					dev.off()
				}
			}
		} 
)

