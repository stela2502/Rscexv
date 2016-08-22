#' @name plot.histograms
#' @aliases plot.histograms,Rscexv-method
#' @rdname plot.histograms-methods
#' @docType methods
#' @description This function plots all genes as histograms to check whether there ar4e clear expression differences in different plates.
#' @param dataObj the Rscexv object
#' @param cuts the cuts are used for the 1D gene groups default=vector('list',1)
#' @param subpath the subpath to plot to ( default = preprocess)
#' @title description of function plot.histograms
#' @export 
setGeneric('plot.histograms', ## Name
		function ( dataObj, cuts=vector('list',1), subpath='preprocess' ) { ## Argumente der generischen Funktion
			standardGeneric('plot.histograms') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plot.histograms', signature = c ('Rscexv'),
		definition = function ( dataObj, cuts=vector('list',1), subpath='preprocess' ) {
			ma <- t(dataObj@data)
			if ( dataObj@wFACS ){
				ma <- rbind( ma, t( dataObj@facs) )
			}
			n <- rownames(ma)
			arrays <- max(dataObj@samples$ArrayID)
			cols <- rainbow(arrays)
			n.cuts <- names(cuts)
			opath = file.path(dataObj@outpath,subpath )
			dir.create(opath, showWarnings = FALSE)
			for ( i in 1:nrow(ma) ) {
				png( file=file.path( opath, paste(n[i],'png',sep='.')),width=800, height=800 )
				temp <- vector('list',arrays)
				m <- NULL
				for (a in 1:arrays ) {
					temp[[a]] <- density(ma[i,which(dataObj@samples$ArrayID == a )])
					m <- c(m,max(temp[[a]]$y))
				}
				#h <- hist(ma[i,],main=n[i], xlab='expression values [raw]', freq=F, col=rgb(0, 1, 0, 0.5), cex.lab = 1.5, breaks = 15, ylim=c(0,max(m)) )
				h <- hist(ma[i,], breaks = 15,plot=F)
				m <- c(m, max(h$density) )
				plot( h, freq=F,main=n[i], col=rgb(0, 1, 0, 0.5), xlab="Ct", cex.lab = 1.5, breaks = 15, ylim=c(0,max(m)) )
				for (a in 1:arrays ) {
					lines( temp[[a]] , col=cols[a], lwd=2)
				}
				pos <- which( n.cuts == n[i] )
				if ( length(pos) > 0 ){
					for (c in 1:length(cuts[[pos]]) ) {
						abline( v= cuts[[pos]][c], col='black', lwd = 3, lty = 2 )
					}
				}
				dev.off()
			}
		} 
)

