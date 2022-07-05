

#' @name file.position
#' @aliases file.position,Rscexv-method
#' @rdname file.position-methods
#' @docType methods
#' @description 
#' @param fname  TEXT MISSING
#' @param previous  TEXT MISSING default=NULL
#' @param keyword  TEXT MISSING default='Array'
#' @title description of function file.position
#' @export 
setGeneric('file.position', ## Name
	function ( fname, previous=NULL, keyword='Array' ) { 
		standardGeneric('file.position')
	}
)

setMethod('file.position', signature = c ('Tool_Pipe'),
	definition = function ( fname, previous=NULL, keyword='Array' ) {
	match <- regexpr( paste("(",keyword,"\\d+)", sep=""), fname , perl = TRUE )
	key = substr(fname,match,  attr(match, "match.length") + match -1 )
	if ( is.null(previous) ){
		previous <- list( 'x' = c( 1, key ) )
		names(previous)[1] = key
	}
	else {
		pos = length(names(previous)) +1
		previous[[pos]] = c( pos, paste(key))
		names(previous)[pos] = key
	}
	previous[[which(names(previous) == key ) ]]
} )



#' @name hist.and.report.range
#' @aliases hist.and.report.range,Rscexv-method
#' @rdname hist.and.report.range-methods
#' @docType methods
#' @description 
#' @param tab  TEXT MISSING
#' @title description of function hist.and.report.range
#' @export 
setGeneric('hist.and.report.range', ## Name
	function (tab) { 
		standardGeneric('hist.and.report.range')
	}
)

setMethod('hist.and.report.range', signature = c ('Tool_Pipe'),
	definition = function (tab) {
	
	tab.r <- tab[-which(tab==999)]
	hist(tab.r)
	range(tab.r)
	
} )
#' @name filter.on.controls
#' @aliases filter.on.controls,Rscexv-method
#' @rdname filter.on.controls-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param ref.nms  TEXT MISSING
#' @param thresh  TEXT MISSING
#' @param howm  TEXT MISSING
#' @title description of function filter.on.controls
#' @export 
setGeneric('filter.on.controls', ## Name
	function (dataObj,ref.nms,thresh,howm) { 
		standardGeneric('filter.on.controls')
	}
)

setMethod('filter.on.controls', signature = c ('Tool_Pipe'),
	definition = function (dataObj,ref.nms,thresh,howm) {
	
	conts <- dataObj$PCR[,ref.nms]
	conts[which(conts<=thresh)] <- 0
	conts[conts>0] <- 1
	
	nums <- apply(conts,1,sum)
	
	dataObj <- remove.samples(dataObj, which(nums<howm) )
	dataObj
	
} )

#' @name auto_problem
#' @aliases auto_problem,Rscexv-method
#' @rdname auto_problem-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @title description of function auto_problem
#' @export 
setGeneric('auto_problem', ## Name
	function ( dataObj ) { 
		standardGeneric('auto_problem')
	}
)

setMethod('auto_problem', signature = c ('Tool_Pipe'),
	definition = function ( dataObj ) {
	arrays <- max(dataObj$ArrayID )
	ret <- matrix ( ncol=arrays, nrow=ncol(dataObj$PCR) )
	calculate_acf <- function (x ) { max(as.numeric( fft(acf( x , lag.max = 80 , plot=F )$acf) ) )  }
	for ( i in 1:arrays ){
		ret[,i] <- apply( as.matrix(data.filtered$PCR )[which(data.filtered$ArrayID == i ),],2, calculate_acf)
	}
	rownames(ret)<-colnames(dataObj$PCR)
	ret
} )
#' @name clus.PCR
#' @aliases clus.PCR,Rscexv-method
#' @rdname clus.PCR-methods
#' @docType methods
#' @description 
#' @param tab  TEXT MISSING
#' @title description of function clus.PCR
#' @export 
setGeneric('clus.PCR', ## Name
	function (tab) { 
		standardGeneric('clus.PCR')
	}
)

setMethod('clus.PCR', signature = c ('Tool_Pipe'),
	definition = function (tab) {
	
	gene.c <- hclust(as.dist(1-cor(tab,method="spearman")),method="ward.D2")$order
	cell.c <- hclust(as.dist(1-cor(t(tab),method="spearman")),method="ward.D2")$order
	
	tabc <- tab[cell.c, gene.c]
	
	brks <- as.vector(c(0,quantile(tab[which(tab!=0)],seq(0,1,by=0.05)),max(tab)))
	
	print(brks)
	print(length(brks))
	
	image.plot(tabc,col=c("darkgrey",bluered((length(brks)-2))),breaks=brks)
} )
#' @name clus.z.PCR
#' @aliases clus.z.PCR,Rscexv-method
#' @rdname clus.z.PCR-methods
#' @docType methods
#' @description 
#' @param tab  TEXT MISSING
#' @param k  TEXT MISSING
#' @title description of function clus.z.PCR
#' @export 
setGeneric('clus.z.PCR', ## Name
	function (tab,k) { 
		standardGeneric('clus.z.PCR')
	}
)

setMethod('clus.z.PCR', signature = c ('Tool_Pipe'),
	definition = function (tab,k) {
	
	gene.c <- hclust(as.dist(1-cor(tab,method="spearman")),method="ward.D2")
	cell.c <- hclust(as.dist(1-cor(t(tab),method="spearman")),method="ward.D2")
	
	cell.k <- cutree(cell.c,k)
	
	tabc <- tab[cell.c$order, gene.c$order]
	
	brks <- as.vector(c(-20,quantile(tab[which(tab!= -20)],seq(0,1,by=0.1)),max(tab)))
	
	
	heatmap.2(t(tabc),breaks=brks,col=c("darkgrey",bluered(length(brks)-2)),trace="none",scale="none")
	cell.k
} )


#' @name update.FACS.rownames
#' @aliases update.FACS.rownames,Rscexv-method
#' @rdname update.FACS.rownames-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param old.names  TEXT MISSING
#' @param new.names  TEXT MISSING
#' @title description of function update.FACS.rownames
#' @export 
setGeneric('update.FACS.rownames', ## Name
	function (dataObj, old.names, new.names ) { 
		standardGeneric('update.FACS.rownames')
	}
)

setMethod('update.FACS.rownames', signature = c ('Tool_Pipe'),
	definition = function (dataObj, old.names, new.names ) {
	if ( ! is.null(dataObj$FACS ) ){
		rownames(dataObj$FACS$data)[ match( old.names, rownames(dataObj$FACS$data)) ] <- new.names
	}
	dataObj
} )
#' @name files.sorted
#' @aliases files.sorted,Rscexv-method
#' @rdname files.sorted-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @param keyword  TEXT MISSING default= 'Array'
#' @title description of function files.sorted
#' @export 
setGeneric('files.sorted', ## Name
	function ( x, keyword = 'Array' ) { 
		standardGeneric('files.sorted')
	}
)

setMethod('files.sorted', signature = c ('Tool_Pipe'),
	definition = function ( x, keyword = 'Array' ) {
	x[order( regmatches(x, regexpr( paste(keyword,"\\d+",sep=''),x ,perl=T)) )]
} )


