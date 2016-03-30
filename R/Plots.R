#' @name plot.histograms
#' @aliases plot.histograms,Rscexv-method
#' @rdname plot.histograms-methods
#' @docType methods
#' @description This function plots all genes as histograms to check whether there ar4e clear expression differences in different plates.
#' @param dataObj the Rscexv object
#' @param cuts the cuts are used for the 1D gene groups default=vector('list',1)
#' @title description of function plot.histograms
#' @export 
setGeneric('plot.histograms', ## Name
		function ( dataObj, cuts=vector('list',1) ) { ## Argumente der generischen Funktion
			standardGeneric('plot.histograms') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plot.histograms', signature = c ('Rscexv'),
		definition = function ( dataObj, cuts=vector('list',1) ) {
			ma <- t(dataObj@data)
			if ( ! is.null(dataObj@facs)){
				ma <- rbind( ma, t( dataObj@facs) )
			}
			n <- rownames(ma)
			arrays <- max(dataObj@samples$ArrayID)
			cols <- rainbow(arrays)
			n.cuts <- names(cuts)
			opath = file.path(dataObj@outpath, 'preprocess')
			dir.create(opath, showWarnings = FALSE)
			for ( i in 1:nrow(ma) ) {
				png( file=paste(opath,"/",n[i],'.png',sep=''),width=800, height=800 )
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

#' @name PCR.heatmap
#' @aliases PCR.heatmap,Rscexv-method
#' @rdname PCR.heatmap-methods
#' @docType methods
#' @description This function plots the heatmaps
#' @param dataObj  TEXT MISSING
#' @param ofile  TEXT MISSING
#' @param title  TEXT MISSING default='Heatmap'
#' @param nmax  TEXT MISSING default=500
#' @param hc.row  TEXT MISSING default=NA
#' @param ColSideColors  TEXT MISSING default=NA
#' @param RowSideColors  TEXT MISSING default=F
#' @param PCR.heatmap  TEXT MISSING
#' @param reorder if yes the data matrix will be reordered to the clusters object (default=F)
#' @title description of function PCR.heatmap
#' @export 
setGeneric('PCR.heatmap', ## Name
		function ( dataObj, ofile,reorder =F,  title='Heatmap', nmax=500, hc.row=NA, ColSideColors=NA, RowSideColors=F,
				width=6, height=6, margins = c(1,11) ,lwid = c( 1,6), lhei=c(1,5), hclustfun = function(c){hclust( c, method='ward.D')}, distfun = function (x) as.dist( 1- cor(t(x), method='pearson') ), Rowv=T, ... ) {## Argumente der generischen Funktion
			standardGeneric('PCR.heatmap') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('PCR.heatmap', signature = c ('Rscexv'),
		definition = function ( dataObj, ofile,reorder =F,  title='Heatmap', nmax=500, hc.row=NA, 
				ColSideColors=NA, RowSideColors=F, width=6, height=6, margins = c(1,11),
				lwid = c( 1,6), lhei=c(1,5), hclustfun = function(c){hclust( c, method='ward.D')}, 
				distfun = function (x) as.dist( 1- cor(t(x), method='pearson') ), Rowv=T, ... ) {
			##plot the heatmap as svg image
			if ( nrow(dataObj@data) > nmax ) {
				stop (paste('No plotting for file ',ofile,'- too many genes selected (',nrow(data),')' ))
			}
			if( nrow(dataObj@data) > 2 ){
				data <- as.matrix(t(dataObj@data))
				if ( reorder ){
					data <- data[,order(dataObj@usedObj[['clusters']])]
				}
				brks <- as.vector(c(-20,quantile(data[which(data!= -20)],seq(0,1,by=0.1)),max(data)))
				#rownames( data ) <- paste( dataObj$genes, dataObj$names)
				if ( is.na(hc.row) ){
					hc.row <- hclustfun(distfun(data)) #hclust( as.dist( 1- cor(t(data), method='spearman')), method='ward')
				}
				dendrogram='both'
				if ( length(grep ('color_groups', ofile)) == 0 ) {
					ma <- data[hc.row$order,]
					dendrogram='both'
				}
				else {
					ma <- data
					dendrogram='row'
				}
				if ( ! RowSideColors==F ) {
					ma <- ma[match(geneGroups[order(geneGroups[,3]),1],rownames(ma)),]
					if ( dendrogram=='both'){
						dendrogram='col'
					}else {
						dendrogram='none'
					}
				}
				if (  plotsvg == 1) {
					devSVG( file=paste(ofile,'_Heatmap.svg',sep='') , width=width, height=height)
					if ( ! is.na(ColSideColors) ) {
						if ( RowSideColors != F) {
							heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), key=F, symkey=FALSE,trace='none', 
									cexRow=0.7,cexCol=0.7, main=title,margins = margins, ColSideColors=ColSideColors, RowSideColors=RowSideColors, Rowv=F,dendrogram=dendrogram,lwid = lwid, lhei=lhei, ... )
						}
						else {
							heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), key=F, symkey=FALSE,
									trace='none', cexRow=0.7,cexCol=0.7, main=title,margins = margins, 
									ColSideColors=ColSideColors, hclustfun = hclustfun, distfun = distfun, Rowv=T,dendrogram=dendrogram,lwid = lwid, lhei=lhei, ...)
						}
					}
					else {
						heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), Rowv=F,  key=F, symkey=FALSE,
								trace='none', cexRow=0.7,cexCol=0.7, main=title,margins = margins,
								hclustfun = hclustfun, distfun = distfun, dendrogram=dendrogram,lwid = lwid, lhei=lhei )
					}
					dev.off()
				}
				if ( nrow(data) > 50 ) {
					png( file=paste(ofile,'_Heatmap.png',sep='') , width=width*150, height=height*250 )
				}
				else {
					png( file=paste(ofile,'_Heatmap.png',sep='') , width=width*150, height=height*200 )
				}
				if ( ! is.na(ColSideColors) ) {
					if ( RowSideColors != F) {
						heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), key=F, symkey=FALSE,trace='none', 
								cexRow=2,cexCol=0.7, main=title,margins = margins, ColSideColors=ColSideColors, RowSideColors=RowSideColors, Rowv=F,dendrogram=dendrogram,lwid = lwid, lhei=lhei, ... )
					}
					else {
						heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), key=F, symkey=FALSE,
								trace='none', cexRow=2,cexCol=0.7, main=title,margins = margins, 
								ColSideColors=ColSideColors, hclustfun = hclustfun, distfun = distfun, Rowv=T,dendrogram=dendrogram,lwid = lwid, lhei=lhei, ...)
					}
				}
				else {
					heatmap.2(as.matrix(ma), breaks=brks,col=c("darkgrey",bluered(length(brks)-2)), Rowv=F,  key=F, symkey=FALSE,
							trace='none', cexRow=2,cexCol=0.7, main=title,margins = margins,
							hclustfun = hclustfun, distfun = distfun, dendrogram=dendrogram,lwid = lwid, lhei=lhei )
				}
				dev.off()
				write.table( cbind ( 'GeneSymbol' = rownames(ma), 'groupsID' = hc.row$order[hc.row$order], ma),file= paste(ofile,'_data4Genesis.xls', sep=''),sep='\t' )
				write ( rownames(ma),file= paste(ofile,'_Genes_in_order.txt',sep='') ,ncolumns = 1 )
			}
			else {
				print ( paste( 'You have less than two genes for the histogram (',nrow(ma),', ',ofile,') '))
			}
			dataObj@usedObj[['expression.hc.row']] = hc.row
			dataObj
		} 
)

#' @name FACS.heatmap
#' @aliases FACS.heatmap,Rscexv-method
#' @rdname FACS.heatmap-methods
#' @docType methods
#' @description Plots the FACS heatmaps
#' @param dataObj  TEXT MISSING
#' @param ofile  TEXT MISSING
#' @param title  TEXT MISSING default='Heatmap'
#' @param nmax  TEXT MISSING default=500
#' @param hc.row  TEXT MISSING default=NA
#' @param ColSideColors  TEXT MISSING default=NA
#' @param RowSideColors  TEXT MISSING default=NA
#' @param FACS.heatmap  TEXT MISSING
#' @param reorder if yes the data matrix will be reordered to the clusters object (default=F)
#' @title description of function FACS.heatmap
#' @export 
setGeneric('FACS.heatmap', ## Name
		function ( dataObj, ofile, title='Heatmap', reorder =F, nmax=500, hc.row=NA, ColSideColors=NA, RowSideColors=NA,width=6, height=6, margins = c(15, 10), hclustfun = function(c){hclust( c, method='ward')}, distfun = function (x) as.dist( 1- cor(t(x), method='pearson') ), ... ){ ## Argumente der generischen Funktion
			standardGeneric('FACS.heatmap') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('FACS.heatmap', signature = c ('Rscexv'),
		definition = function ( dataObj, ofile, title='Heatmap', reorder =F, nmax=500, hc.row=NA, ColSideColors=NA, RowSideColors=NA,
				width=6, height=6, margins = c(15, 10), hclustfun = function(c){hclust( c, method='ward')}, distfun = function (x) as.dist( 1- cor(t(x), method='pearson') ), ... ) {
			##plot the heatmap as svg image
			if ( nrow(dataObj@facs) > nmax ) {
				stop (paste('No plotting for file ',ofile,'- too many genes selected (',nrow(data),')' ))
			}
			if( nrow(dataObj@facs) > 2 ){
				data <- as.matrix(t(dataObj@facs))
				if ( reorder ){
					data <- data[,order(dataObj@usedObj[['clusters']])]
				}
				for ( i in 1:2 ){
					#rownames( data ) <- paste( dataObj$genes, dataObj$names)
					if ( i == 1 && plotsvg == 1 ) {
						devSVG( file=paste(ofile,'_Heatmap.svg',sep='') , width=width, height=height)
					}
					else {
						png( file=paste(ofile,'_Heatmap.png',sep='') , width=width*150, height=nrow( data ) * 15 + 400 )
					}
					if ( is.na(hc.row) ){
						
						hc.row <- hclustfun(distfun(data)) #hclust( as.dist( 1- cor(t(data), method='spearman')), method='ward')
					}
					ma <- data[hc.row$order,]
					if ( ! is.na(RowSideColors) ) {
						RowSideColors <- RowSideColors[ hc.row$order ]
					} 
					if ( ! is.na(ColSideColors) ) {
						if ( ! is.na(RowSideColors)) {
							heatmap.2(as.matrix(ma), col=bluered, Rowv=F,  key=F, symkey=FALSE,
									trace='none', cexRow=2,cexCol=0.7, main=title,margins = margins, 
									ColSideColors=ColSideColors, RowSideColors=RowSideColors, hclustfun = hclustfun, distfun = distfun,dendrogram='both', ... )
						}
						else {
							heatmap.2(as.matrix(ma), col=bluered, Rowv=T,  key=F, symkey=FALSE,
									trace='none', cexRow=2,cexCol=0.7, main=title,margins = margins, 
									ColSideColors=ColSideColors, hclustfun = hclustfun, distfun = distfun,dendrogram='both',... )
						}
					}
					else {
						heatmap.2(as.matrix(ma), col=bluered, Rowv=F,  key=F, symkey=FALSE,
								trace='none', cexRow=2,cexCol=0.7, main=title,margins = margins,
								hclustfun = hclustfun, distfun = distfun, ... )
					}
					dev.off()
				}
				write.table( cbind ( 'GeneSymbol' = rownames(ma), 'groupsID' = hc.row$order[hc.row$order], ma),file= paste(ofile,'_data4Genesis.txt', sep=''),sep='\t', row.names=F, quote=F  )
				write ( rownames(ma),file= paste(ofile,'_Genes_in_order.txt',sep='') ,ncolumns = 1 )
			}
			else {
				print ( paste( 'You have less than two genes for the histogram (',nrow(ma),', ',ofile,') '))
			}
			
			dataObj@usedObj[['facs.hc.row']] = hc.row
		} 
)

#' @name plot.violines
#' @aliases plot.violines,Rscexv-method
#' @rdname plot.violines-methods
#' @docType methods
#' @description This fucntion converts the table into a format, that can be fead to the vioplot function
#' @param x the Rscexv object
#' @param groups.n  TEXT MISSING
#' @param clus  TEXT MISSING
#' @param boot  TEXT MISSING default= 1000
#' @param plot.neg  TEXT MISSING default=FALSE
#' @param mv  TEXT MISSING default=-20
#' @title description of function plot.violines
#' @export 
setGeneric('plot.violines', ## Name
		function ( x, groups.n, clus, boot = 1000, plot.neg=FALSE, mv=-20) { ## Argumente der generischen Funktion
			standardGeneric('plot.violines') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plot.violines', signature = c ('Rscexv'),
		definition = function ( x, groups.n, clus, boot = 1000, plot.neg=FALSE, mv=-20) {
			ma <- NULL
			if ( x@wFACS ){
				ma <- t(cbind( x@usedObj[['for.plot']], x@facs ))
			}else{
				ma <- t(x@usedObj[['for.plot']])
			}
			
			n <- rownames(ma)
			cols = rainbow( groups.n )
			s <-  split(seq(ncol(ma)), clus)
			for ( i in 1:nrow( ma ) ) {
			#	print (paste( 'plot.violines working on gene', n[i] ) )
				png( file=paste(n[i],'.png',sep=''), width=800,height=800)
				#create color info
				lila <- lapply(s ,function(x) { ma[ i, x] } )
				lila$names <- as.vector(unlist(lapply( lila, 
							function(x) { 
								lila$names <- c( lila$names, paste(length(which(x != mv)), length(x) ,sep='/' ) )
							} ) ))
				if ( ! plot.neg ){
					lila = lapply( lila, function(x) { x[which(x != mv)] = NA; x } )
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
					devSVG( file=paste(n[i],'.svg',sep=''), width=6,height=6)
					lila$cex.axis=0.5
					try(do.call(vioplot,lila), silent=F )
					dev.off()
				}
			}
		} 
)

#' @name vioplot
#' @aliases vioplot,Rscexv-method
#' @rdname vioplot-methods
#' @docType methods
#' @description the vioplot function patched to allow different colors for the different violines.
#' @param x  TEXT MISSING
#' @param ...  TEXT MISSING
#' @param range  TEXT MISSING default= 1.5
#' @param h  TEXT MISSING default= NULL
#' @param ylim  TEXT MISSING default= NULL
#' @param names  TEXT MISSING default= NULL
#' @param vioplot  TEXT MISSING
#' @title description of function vioplot
#' @export 
setGeneric('vioplot', ## Name
		function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
				horizontal = FALSE, col = 'magenta', border = 'black', lty = 1, 
				lwd = 1, rectCol = 'black', colMed = 'white', pchMed = 19, 
				at, add = FALSE, wex = 1, drawRect = TRUE, main=NULL, cex.axis=1, neg=NULL) { ## Argumente der generischen Funktion
			standardGeneric('vioplot') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)
setMethod('vioplot', signature = c ('numeric'),
		definition = function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
		horizontal = FALSE, col = 'magenta', border = 'black', lty = 1, 
		lwd = 1, rectCol = 'black', colMed = 'white', pchMed = 19, 
		at, add = FALSE, wex = 1, drawRect = TRUE, main=NULL, cex.axis=1, neg=NULL) 
{
	datas <- list(x, ...)
	n <- length(datas)
	if (missing(at)) 
		at <- 1:n
	upper <- vector(mode = 'numeric', length = n)
	lower <- vector(mode = 'numeric', length = n)
	q1 <- vector(mode = 'numeric', length = n)
	q3 <- vector(mode = 'numeric', length = n)
	med <- vector(mode = 'numeric', length = n)
	base <- vector(mode = 'list', length = n)
	height <- vector(mode = 'list', length = n)
	baserange <- c(Inf, -Inf)
	args <- list(display = 'none')
	if (!(is.null(h))) 
		args <- c(args, h = h)
	names.2 <- NULL
	for (i in 1:n) {
		data <- datas[[i]][ is.na(datas[[i]]) ==F ]
		if ( ! is.null(neg)) {
			names.2 <- c ( names.2, paste( length(which( datas[[i]] != neg )),"/",length(data),sep='') )
		}else {
			names.2 <- c ( names.2, paste( length(data),"/",length(datas[[i]]),sep='') )
		}
		if ( length(data) == 0) {
			data <- c(0)
		}
		data.min <- min(data)
		data.max <- max(data)
		
		if ( data.min == data.max ) {
			next;
		}
		q1[i] <- quantile(data, 0.25)
		q3[i] <- quantile(data, 0.75)
		med[i] <- median(data)
		
		iqd <- q3[i] - q1[i]
		upper[i] <- min(q3[i] + range * iqd, data.max)
		lower[i] <- max(q1[i] - range * iqd, data.min)
		est.xlim <- c(min(lower[i], data.min), max(upper[i], 
						data.max))
		smout <- do.call('sm.density', c(list(data, xlim = est.xlim), 
						args))
		hscale <- 0.4/max(smout$estimate) * wex
		base[[i]] <- smout$eval.points
		height[[i]] <- smout$estimate * hscale
		t <- range(base[[i]])
		baserange[1] <- min(baserange[1], t[1])
		baserange[2] <- max(baserange[2], t[2])
	}
	if (!add) {
		xlim <- if (n == 1) 
					at + c(-0.5, 0.5)
				else range(at) + min(diff(at))/2 * c(-1, 1)
		if (is.null(ylim)) {
			ylim <- baserange
		}
	}
	if ( ! is.null(names)) {
		label <- names
	}
	else if (is.null(names.2)) {
		label <- 1:n
	}
	else {
		label <- names.2
	}
	boxwidth <- 0.05 * wex
	if ( length( col ) == 1 ){
		col = rep(col,n)
	}
	if ( length(col ) < n ) {
		stop(paste(length(col),'colors are too view to color',n,'data sets'))
	} 
	if (!add) 
		plot.new()
	if (!horizontal) {
		if (!add) {
			plot.window(xlim = xlim, ylim = ylim)
			axis(2, cex.axis=cex.axis)
			axis(1, at = at, label = label, cex.axis=cex.axis)
		}
		box()
		for (i in 1:n) {
			polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
					c(base[[i]], rev(base[[i]])), col = col[i], border = border, 
					lty = lty, lwd = lwd)
			if (drawRect) {
				lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
						lty = lty)
				rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
						q3[i], col = rectCol)
				points(at[i], med[i], pch = pchMed, col = colMed)
			}
			else{
				lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
						lty = lty)
				lines(c(at[i]- boxwidth/2, at[i] + boxwidth/2), c(lower[i], lower[i]), lwd = lwd, 
						lty = lty)
				lines( c(at[i]- boxwidth/2, at[i] + boxwidth/2), c(upper[i], upper[i]), lwd = lwd, 
						lty = lty)
				points(at[i], med[i], pch = pchMed, col = colMed, cex=2)
			}
		}
		
	}
	else {
		if (!add) {
			plot.window(xlim = ylim, ylim = xlim)
			axis(1,cex.axis =cex.axis)
			axis(2, at = at, label = label, cex.axis=cex.axis)
		}
		box()
		for (i in 1:n) {
			polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
							rev(at[i] + height[[i]])), col = col[i], border = border, 
					lty = lty, lwd = lwd)
			if (drawRect) {
				lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
						lty = lty)
				rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
								boxwidth/2, col = rectCol)
				points(med[i], at[i], pch = pchMed, col = colMed)
			}
		}
	}
	if ( ! is.null(main) ){
		title( main, cex.main = 2)
	}
	invisible(list(upper = upper, lower = lower, median = med, 
					q1 = q1, q3 = q3))
}

)

#' @name col4bean
#' @aliases col4bean,Rscexv-method
#' @rdname col4bean-methods
#' @docType methods
#' @description get the color information for the beans
#' @param x  TEXT MISSING
#' @param tic  TEXT MISSING default='black'
#' @title description of function col4bean
#' @export 
setGeneric('col4bean', ## Name
		function (x, tic='black') { ## Argumente der generischen Funktion
			standardGeneric('col4bean') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('col4bean', signature = c ('character'),
		definition = function (x, tic='black') {
			ret <- list()
			for ( i in 1:length(x) ){
				ret[[i]] <- c(x[i], tic )
			}
			ret
		}
)


#' @name plot.beans
#' @aliases plot.beans,Rscexv-method
#' @rdname plot.beans-methods
#' @docType methods
#' @description plot the beanplots
#' @param x the Rscexv object
#' @param groups.n  TEXT MISSING
#' @param clus  TEXT MISSING
#' @param boot  TEXT MISSING default= 1000
#' @param plot.neg  TEXT MISSING default=TRUE
#' @param mv  TEXT MISSING default=-20
#' @title description of function plot.beans
#' @export 
setGeneric('plot.beans', ## Name
		function ( x, groups.n, clus, boot = 1000, plot.neg=TRUE, mv=-20 ) { ## Argumente der generischen Funktion
			standardGeneric('plot.beans') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plot.beans', signature = c ('Rscexv'),
		definition = function ( x, groups.n, clus, boot = 1000, plot.neg=TRUE, mv=-20 ) {
			ma <- NULL
			if ( x@wFACS ){
				ma <- t(cbind( x@usedObj[['for.plot']], x@facs ))
			}else{
				ma <- t(x@usedObj[['for.plot']])
			}
			n <- rownames(ma)
			cols = col4bean(rainbow( groups.n ))	
			s <-  split(seq(ncol(ma)), clus)
			for ( i in 1:nrow( ma ) ) {
				#print (paste( 'plot.beans working on gene', n[i] ) )
				png( file=paste(n[i],'.png',sep=''), width=800,height=800)
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
				lila$main <- n[i]
				lila$what <- c(1,1,0,1) ## not plot medain line
				lila$col <- cols
				try(do.call(beanplot,lila), silent=F )
				dev.off()
				if ( plotsvg == 1 ) {
					devSVG( file=paste(n[i],'.svg',sep=''), width=6,height=6)
					#lila$cex.axis=0.5
					try(do.call(beanplot,lila), silent=F )
					dev.off()
				}
			}
		} 
)


