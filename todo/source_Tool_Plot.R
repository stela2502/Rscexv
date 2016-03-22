
#' @name collapsData
#' @aliases collapsData,Rscexv-method
#' @rdname collapsData-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param method  TEXT MISSING default='median'
#' @title description of function collapsData
#' @export 
setGeneric('collapsData', ## Name
	function ( dataObj, method='median' ) { ## Argumente der generischen Funktion
		standardGeneric('collapsData') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('collapsData', signature = c ('Tool_Plot'),
	definition = function ( dataObj, method='median' ) {
	ret <- matrix ( ncol= max(dataObj$clusters),nrow= ncol(data$z$PCR))
	colnames(ret) <- paste('Cluster', 1:ncol(ret))
	rownames(ret) <- colnames(dataObj$z$PCR)
	for ( genecol in 1:nrow(ret) ) { ## genes
		for ( cluster in 1:ncol(ret) ){
			if ( method == 'median' ){
				ret[genecol,cluster] = median(dataObj$z$PCR[which(dataObj$clusters == cluster),genecol ] )
			}else if ( method == 'mean' ){
				ret[genecol,cluster] = mean(dataObj$z$PCR[which(dataObj$clusters == cluster),genecol ] )
			}else if ( method == 'var' ){
				ret[genecol,cluster] = var(dataObj$z$PCR[which(dataObj$clusters == cluster),genecol ] )
			}else if ( method == 'quantile70' ){
				ret[genecol,cluster] = as.vector(quantile(dataObj$z$PCR[which(dataObj$clusters == cluster),genecol ],0.7 ))
			}
			
			else{
				stop('method not implemented!')
			}
			
		}
	}
	if ( length( which(apply(ret,1,var) == 0))> 0 ){
		ret <- ret[-which(apply(ret,1,var) == 0),]
	}
	ret
} )
#' @name collapsed_heatmaps
#' @aliases collapsed_heatmaps,Rscexv-method
#' @rdname collapsed_heatmaps-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param what  TEXT MISSING default='PCR'
#' @param functions  TEXT MISSING default= c('median'
#' @param 'mean'  TEXT MISSING
#' @param 'var'  TEXT MISSING
#' @param 'quantile70')  TEXT MISSING
#' @title description of function collapsed_heatmaps
#' @export 
setGeneric('collapsed_heatmaps', ## Name
	function ( dataObj, what='PCR', functions = c('median', 'mean', 'var', 'quantile70' ) ) { ## Argumente der generischen Funktion
		standardGeneric('collapsed_heatmaps') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('collapsed_heatmaps', signature = c ('Tool_Plot'),
	definition = function ( dataObj, what='PCR', functions = c('median', 'mean', 'var', 'quantile70' ) ) {
	if ( ! is.vector(functions) ){
		functions = c( functions )
	}
	data <- NULL
	if ( what == 'PCR' ){
		data = dataObj$z$PCR
	}else if ( what =='FACS' ){
		data = dataObj$FACS
	}else {
		stop('collapsed_heatmaps can only collaps PCR or FACS data' )
	}
	if ( !is.null(data)){
		for( i in 1:length(functions)) {
			reduced.filtered <- collapsData( data ,method=functions[i])
			PCR.heatmap ( dataObj= list( data= reduced.filtered, genes= rownames(reduced.filtered)), ofile=paste(what,functions[i],sep="_") , margins=c(3,10),ColSideColors=rainbow(max(data$clusters)), Colv=F, Rowv=F, title=functions[i],RowSideColors=1)
		}
	}
} )

#' @name vioplot
#' @aliases vioplot,Rscexv-method
#' @rdname vioplot-methods
#' @docType methods
#' @description 
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
			at, add = FALSE, wex = 1, drawRect = TRUE, main=NULL, cex.axis=1, neg=NULL)  ) { ## Argumente der generischen Funktion
		standardGeneric('vioplot') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('vioplot', signature = c ('Tool_Plot'),
	definition = function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL,
		horizontal = FALSE, col = 'magenta', border = 'black', lty = 1, 
		lwd = 1, rectCol = 'black', colMed = 'white', pchMed = 19, 
		at, add = FALSE, wex = 1, drawRect = TRUE, main=NULL, cex.axis=1, neg=NULL)  )
#' @name calc.ann
#' @aliases calc.ann,Rscexv-method
#' @rdname calc.ann-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @param groups  TEXT MISSING
#' @title description of function calc.ann
#' @export 
setGeneric('calc.ann', ## Name
	function (x, groups ) { ## Argumente der generischen Funktion
		standardGeneric('calc.ann') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('calc.ann', signature = c ('Tool_Plot'),
	definition = function (x, groups ) {
	a <- min(TukeyHSD(aov(formula = as.vector(x) ~ as.factor(groups)))$"as.factor(groups)"[,4] )
	if ( a < 1e-26 ) {
		a = 1e-26
	}
	a
} )
#' @name get.GOI
#' @aliases get.GOI,Rscexv-method
#' @rdname get.GOI-methods
#' @docType methods
#' @description 
#' @param ma  TEXT MISSING
#' @param group  TEXT MISSING
#' @param exclude  TEXT MISSING default= NULL
#' @title description of function get.GOI
#' @export 
setGeneric('get.GOI', ## Name
	function ( ma, group, exclude = NULL ) { ## Argumente der generischen Funktion
		standardGeneric('get.GOI') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('get.GOI', signature = c ('Tool_Plot'),
	definition = function ( ma, group, exclude = NULL ) {
	d <- apply( ma, 2, calc.ann, groups=group )
	d <- d*ncol(ma) * max(group)
	ret <- d[which(d< 0.05 )]
	ret
} )

#' @name col4bean
#' @aliases col4bean,Rscexv-method
#' @rdname col4bean-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @param tic  TEXT MISSING default='black'
#' @title description of function col4bean
#' @export 
setGeneric('col4bean', ## Name
	function (x, tic='black') { ## Argumente der generischen Funktion
		standardGeneric('col4bean') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('col4bean', signature = c ('Tool_Plot'),
	definition = function (x, tic='black') {
	ret <- list()
	for ( i in 1:length(x) ){
		ret[[i]] <- c(x[i], tic )
	}
	ret
} )
#' @name plot.beans
#' @aliases plot.beans,Rscexv-method
#' @rdname plot.beans-methods
#' @docType methods
#' @description 
#' @param ma  TEXT MISSING
#' @param groups.n  TEXT MISSING
#' @param clus  TEXT MISSING
#' @param boot  TEXT MISSING default= 1000
#' @param plot.neg  TEXT MISSING default=TRUE
#' @param mv  TEXT MISSING default=-20
#' @title description of function plot.beans
#' @export 
setGeneric('plot.beans', ## Name
	function ( ma, groups.n, clus, boot = 1000, plot.neg=TRUE, mv=-20 ) { ## Argumente der generischen Funktion
		standardGeneric('plot.beans') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('plot.beans', signature = c ('Tool_Plot'),
	definition = function ( ma, groups.n, clus, boot = 1000, plot.neg=TRUE, mv=-20 ) {
	ma <- t(ma)
	n <- rownames(ma)
	cols = col4bean(rainbow( groups.n ))
	for ( i in 1:nrow( ma ) ) {
		print (paste( 'plot.beans working on gene', n[i] ) )
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
} )
#' @name plot.violines
#' @aliases plot.violines,Rscexv-method
#' @rdname plot.violines-methods
#' @docType methods
#' @description 
#' @param ma  TEXT MISSING
#' @param groups.n  TEXT MISSING
#' @param clus  TEXT MISSING
#' @param boot  TEXT MISSING default= 1000
#' @param plot.neg  TEXT MISSING default=FALSE
#' @param mv  TEXT MISSING default=-20
#' @title description of function plot.violines
#' @export 
setGeneric('plot.violines', ## Name
	function ( ma, groups.n, clus, boot = 1000, plot.neg=FALSE, mv=-20) { ## Argumente der generischen Funktion
		standardGeneric('plot.violines') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('plot.violines', signature = c ('Tool_Plot'),
	definition = function ( ma, groups.n, clus, boot = 1000, plot.neg=FALSE, mv=-20) {
	ma <- t(ma)
	n <- rownames(ma)
	cols = rainbow( groups.n )
	for ( i in 1:nrow( ma ) ) {
		print (paste( 'plot.violines working on gene', n[i] ) )
		png( file=paste(n[i],'.png',sep=''), width=800,height=800)
		#create color info
		lila <- vector('list', groups.n)
		lila$names <- NULL
		for( a in 1:groups.n){
			lila[[a]]=ma[i,which(clus == a)]
			lila$names <- c( lila$names, paste(length(which(lila[[a]] != mv)), length(lila[[a]]) ,sep='/' ) )
			if ( ! plot.neg ){
				lila[[a]][which(lila[[a]] == mv)] <- NA
			}
		}
		names(lila)[1]= 'x'
		lila$col= cols
		lila$main=n[i]
		if ( ! plot.neg ) {
			lila$neg = mv
		}else{
			lila$neg = NULL
		}
		
		lila$h = 0.3
		if ( ! is.null(plot.neg) ){
			lila$drawRect = FALSE
		}
		try(do.call(vioplot,lila), silent=F )
		dev.off()
		if ( plotsvg == 1 ) {
			devSVG( file=paste(n[i],'.svg',sep=''), width=6,height=6)
			lila$cex.axis=0.5
			try(do.call(vioplot,lila), silent=F )
			dev.off()
		}
	}
} )
