
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



