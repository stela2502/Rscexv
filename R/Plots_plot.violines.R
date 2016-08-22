#' @name plot.violines
#' @aliases plot.violines,Rscexv-method
#' @rdname plot.violines-methods
#' @docType methods
#' @description This fucntion converts the table into a format, that can be fead to the vioplot function
#' @param x the Rscexv object
#' @param groups.n  TEXT MISSING
#' @param clus  TEXT MISSING
#' @param plot.neg  TEXT MISSING default=FALSE
#' @param mv  TEXT MISSING default=-20
#' @param subpath the path to write the figures to (default = '')
#' @param col a colour vector for colouring the violines (default=NULL  -> rainbow )
#' @title description of function plot.violines
#' @export 
setGeneric('plot.violines', ## Name
		function ( x, groups.n, clus, plot.neg=FALSE, mv=-20, subpath='', names=NULL, col=NULL) { ## Argumente der generischen Funktion
			standardGeneric('plot.violines') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plot.violines', signature = c ('Rscexv'),
		definition = function ( x, groups.n, clus,  plot.neg=FALSE, mv=-20, subpath='', names=NULL, col=NULL) {
			ma <- NULL
			if ( x@wFACS ){
				ma <- t(cbind( x@usedObj[['for.plot']], x@facs ))
			}else{
				ma <- t(x@usedObj[['for.plot']])
			}
			if ( subpath != '' ){
				system( paste('mkdir',file.path(x@outpath,subpath ) ) )
			}
			n <- rownames(ma)
			if ( is.null(col)){
				col = rainbow( groups.n )
			}
			cols = col
			
			s <-  split(seq(ncol(ma)), clus)
			for ( i in 1:nrow( ma ) ) {
				fn <- file.path(x@outpath,subpath, n[i] )
			#	print (paste( 'plot.violines working on gene', n[i] ) )
				png( file=paste(fn,'.png',sep=''), width=800,height=800)
				#create color info
				lila <- lapply(s ,function(x) { ma[ i, x] } )
				lila$names <- as.vector(unlist(lapply( lila, 
							function(x) { 
								lila$names <- c( lila$names, paste(length(which(x != mv)), length(x) ,sep='/' ) )
							} ) ))
				if ( ! plot.neg ){
					lila = lapply( lila, function(x) { x[which(x != mv)] = NA; x } )
				}
				if ( ! is.null(names) && length(names) == length(lila$names) ){
					lila$names= paste( names, lila$names)
				}
				names(lila)[1]= 'x'
				lila$col= cols
				lila$main=n[i]
				if ( ! plot.neg ) {
					lila$neg = mv
					lila$drawRect = FALSE
				}else{
					lila$neg = NULL
				}
				lila$h = 0.3
				try(do.call(vioplot,lila), silent=F )
				dev.off()
				if ( plotsvg == 1 ) {
					devSVG( file=paste(fn,'.svg',sep=''), width=6,height=6)
					lila$cex.axis=0.5
					try(do.call(vioplot,lila), silent=F )
					dev.off()
				}
			}
		} 
)

