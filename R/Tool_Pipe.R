#' @name force.unique.sample
#' @aliases force.unique.sample,Rscexv-method
#' @rdname force.unique.sample-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @title description of function force.unique.sample
#' @export 
setGeneric('force.unique.sample', ## Name
	function ( x ) { ## Argumente der generischen Funktion
		standardGeneric('force.unique.sample') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('force.unique.sample', signature = c ('Tool_Pipe'),
	definition = function ( x ) {
	ret <- NULL
	last <- ""
	use <- NULL
	id <- 1
	repl <- vector(length=length(x))
	for ( i in 1:length(x) ){
		if ( is.null(ret) ){
			last = x[i]
			ret <- c( last )
			use  <- last
		}
		else if ( x[i] != last ){
			last = x[i]
			if ( ! is.na(match( last, ret )) ){
				use <- paste(last,"_",id, sep = '')
				ret <- c(ret , use )
				id = id + 1
			}else {
				use  <- last
				ret <- c(ret,  last )
			}
		}
		repl[i] <- use
	}
	l <- list( 'unique' = ret, 'replacement' = repl )
	l
} )
#' @name force.absolute.unique.sample
#' @aliases force.absolute.unique.sample,Rscexv-method
#' @rdname force.absolute.unique.sample-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @title description of function force.absolute.unique.sample
#' @export 
setGeneric('force.absolute.unique.sample', ## Name
	function ( x) { ## Argumente der generischen Funktion
		standardGeneric('force.absolute.unique.sample') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('force.absolute.unique.sample', signature = c ('Tool_Pipe'),
	definition = function ( x) {
	last = ''
	ret <- vector(length=length(x))
	
	for ( i in 1:length(x) ){
		if ( is.null(ret) ){
			last = x[i]
			ret[i] <- last
		}
		else{
			last = x[i]
			if ( ! is.na(match( last, ret )) ){
				last <- paste(last,"_",sum( ! is.na(match( x[1:i], last )))-1, sep = '')
			}
			ret[i] <- last
		}
	}
	ret
} )
#' @name read.PCR
#' @aliases read.PCR,Rscexv-method
#' @rdname read.PCR-methods
#' @docType methods
#' @description 
#' @param fname  TEXT MISSING
#' @param use_pass_fail  TEXT MISSING default=T
#' @title description of function read.PCR
#' @export 
setGeneric('read.PCR', ## Name
	function (fname,use_pass_fail=T) { ## Argumente der generischen Funktion
		standardGeneric('read.PCR') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('read.PCR', signature = c ('Tool_Pipe'),
	definition = function (fname,use_pass_fail=T) {
	
	## we also support tab separated files containing only genes as columns and samples as rows.
	
	top20 <- readLines(fname,20)
	if ( length(grep(";", top20 )) > 0 ) {
		### CRAP!!! this is the file output from the windows software!
		print ( paste ( "sed -e 's/,/./g' ", fname, " | sed -e 's/;/,/g' -- > ",fname,".tmp",sep='' ) )
		print ( paste ( "mv ",fname,".tmp ",fname , sep='' ) )
		system ( paste ( "sed -e 's/,/./g' ", fname, " | sed -e 's/;/,/g' -- > ",fname,".tmp",sep='' ) )
		system ( paste ( "mv ",fname,".tmp ",fname , sep='' ) )
		top20 <- readLines(fname,20)
	}
	
	av <- grep("Application Version",top20)
	
	ftab <- NULL
	
	line.start <- grep("^ID",top20)
	options(show.error.messages = FALSE)
	test <- matrix(ncol=0,nrow=0)
	try (test <- read.delim(fname,header=T))
	options(show.error.messages = TRUE)
	if ( ncol(test) == 0 ){
		options(show.error.messages = FALSE)
		test <- matrix(ncol=0,nrow=0)
		try (test <- read.delim(fname,header=T,sep=','))
		options(show.error.messages = TRUE)
	}
	if ( ncol(test) > 1 ){
		#rownames(test) <- test[,1]
		#test <- test[,-1]
		rownames(test) <- force.unique.sample(as.vector(test[,1]))$unique
		ftab <- as.matrix(test[,-1])
	}
	else if(length(av)>0 && length(line.start)!= 0 ){
		
		
		tab <- read.delim(fname,skip=(line.start),header=F,sep=",")
		colnames(tab) <- make.names(paste( unlist(strsplit(top20[line.start-1],',')), unlist(strsplit(top20[line.start],',') ), sep='_' ),unique=T)
		#
		#browser()
		fin.wells <- force.unique.sample ( as.vector(tab$Sample_Name))
		#fin.wells <- unique(as.vector(tab$Name))
		tab$Sample_Name <- fin.wells$replacement
		fin.wells <- fin.wells$unique
		gnameID <- grep ('Name', colnames(tab))[2]
		for(i in 1:length(fin.wells)){
			
			rq.ind <- which(tab[,2]==fin.wells[i])
			
			cts <- matrix(tab$Ct_Value[rq.ind],nrow=1)
			gname <- as.vector(tab[rq.ind,gnameID])
			if ( use_pass_fail ){
				cts[,which(tab$Ct_Call[rq.ind] == 'Fail' )] <- 999
			}
			colnames(cts) <- gname
			ftab <- rbind(ftab,cts)
			
		}
		
		rownames(ftab) <- fin.wells
		dups <- colnames(ftab)[ which(duplicated(colnames(ftab))==T)]
		if ( length(dups) > 0 ){
			for(i in 1:length(dups)){
				dup.ind <- which(colnames(ftab)==dups[i])
				colnames(ftab)[dup.ind] <- paste(colnames(ftab)[dup.ind],1:length(dup.ind),sep="")
				
			}
		}
		
	}
	else if ( is.null(ftab) && length(av)>0 ){
		ftab <- read.PCR.heatmap( fname , use_pass_fail )
	}
	else{
		ftab <- read.delim(fname,header=T,row.names=1)
	}
	for ( i in 1:ncol(ftab) ) {
		ftab[,i] <- as.numeric(as.vector(ftab[,i]))
	}
	ftab
} )
#' @name read.PCR.heatmap
#' @aliases read.PCR.heatmap,Rscexv-method
#' @rdname read.PCR.heatmap-methods
#' @docType methods
#' @description 
#' @param fname  TEXT MISSING
#' @param use_pass_fail  TEXT MISSING default=T
#' @title description of function read.PCR.heatmap
#' @export 
setGeneric('read.PCR.heatmap', ## Name
	function (fname,use_pass_fail=T) { ## Argumente der generischen Funktion
		standardGeneric('read.PCR.heatmap') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('read.PCR.heatmap', signature = c ('Tool_Pipe'),
	definition = function (fname,use_pass_fail=T) {
	
	top20 <- readLines(fname)
	
	ftab <- NULL
	
	line.start <- grep(",,1,2,3,4,5,6,7",top20)
	
	if ( length(line.start) == 0){
		system( paste('echo "File format error: file',fname, 'not readable\n" >> R_file_read_error.txt' ) )
	}
	else{
		tmp <- read.delim(fname,sep=',')
		line.end <- grep("Quality Results,,,,,,,",top20)-1[1]
		while (tmp[line.end,2] == '' ){
			line.end = line.end -1 
		}
		ftab <- tmp[(line.start[1]+1):(line.end),3:ncol(tmp)]
		rownames(ftab) <- tmp[(line.start[1]+1):(line.end),2]
		colnames(ftab) <- as.vector(t(tmp[line.start[1],3:ncol(tmp)]))
		for ( i in 1:ncol(ftab) ) {
			ftab[,i] <- as.numeric(as.vector(ftab[,i]))
		}
		if ( use_pass_fail ){
			## now we read the quality control data
			line.end <-nrow(tmp)
			while ( length(grep( 'Pass',t(tmp[line.end,3:ncol(tmp)]) ) ) == 0 ){ 
				line.end = line.end -1 
			}
			qval <- tmp[(line.start[2]+1):line.end,3:ncol(tmp)]
			rownames(qval) <- tmp[(line.start[2]+1):line.end,2]
			colnames(qval) <- as.vector(t(tmp[line.start[2],3:ncol(tmp)]))
			for ( i in 1:ncol(ftab) ) {
				ftab[which(qval[,i]=='Fail'),i] <- 999
			}
		}
		as.matrix(ftab)
	}
} )
#' @name read.FACS
#' @aliases read.FACS,Rscexv-method
#' @rdname read.FACS-methods
#' @docType methods
#' @description 
#' @param fname  TEXT MISSING
#' @param use_pass_fail  TEXT MISSING default=T
#' @title description of function read.FACS
#' @export 
setGeneric('read.FACS', ## Name
	function (fname,use_pass_fail=T) { ## Argumente der generischen Funktion
		standardGeneric('read.FACS') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('read.FACS', signature = c ('Tool_Pipe'),
	definition = function (fname,use_pass_fail=T) {
	ftab <- matrix(ncol=0, nrow=1)
	try( ftab <- read.PCR  ( fname, use_pass_fail ), silent=T )
	if ( ncol(ftab) == 0 ){
		
		top20 <- readLines(fname)
		
		ftab <- NULL
		
		line.start <- grep("^Well,",top20)[1]
		tab.pre <- read.delim(fname,skip=(line.start-1),header=T,sep=",",row.names=1)
		if ( length(grep ('Min$', colnames( tab.pre))) > 0 ){
			ftab <- tab.pre[, grep ('Min$', colnames( tab.pre)) ]
		}
		else { ## I suppose the file contains only data columns!
			ftab <- tab.pre
		}
		if (length(grep ('All.Events.', colnames( ftab ))) > 0 ) {
			ftab <- ftab[, grep ('All.Events.', colnames( ftab )) ]
			colnames( ftab ) <- str_replace_all( colnames( ftab ), 'All.Events.', '' )
		}
	}
	as.matrix(ftab)
} )
#' @name read.FACS.set
#' @aliases read.FACS.set,Rscexv-method
#' @rdname read.FACS.set-methods
#' @docType methods
#' @description 
#' @param fnames  TEXT MISSING
#' @title description of function read.FACS.set
#' @export 
setGeneric('read.FACS.set', ## Name
	function (fnames) { ## Argumente der generischen Funktion
		standardGeneric('read.FACS.set') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('read.FACS.set', signature = c ('Tool_Pipe'),
	definition = function (fnames) {
	etab <-NULL
	order <- NULL
	if ( length(fnames) > 0) {
		for(i in 1:length(fnames)){
			
			ttab <- read.FACS(fnames[i])
			ttab <- ttab[,order(colnames(ttab))]
			if ( length( grep( "P\\d+$", rownames(ttab))) == 0 ) {
				rownames(ttab) <- paste(rownames(ttab),".P",i-1,sep="")
			}
			## check whether the gene names are axactly the same
			if ( is.null(etab)){
				etab <-ttab
				order <- c(order, rep(i, nrow(ttab) ) )
			}
			else {
				if ( ! identical ( colnames(etab), colnames(ttab))) {
					system ( paste('echo "file', fnames[i],'was rejected due to Gene Symbol mismatches!" >> R_file_read_error.txt') )
				}
				else {
					etab <- rbind(etab,ttab)
					order <- c(order, rep(i, nrow(ttab) ) )
				}
			}
		}
	}
	list ( data = scale.FACS.data(etab), order = order )
} )
#' @name scale.FACS.data
#' @aliases scale.FACS.data,Rscexv-method
#' @rdname scale.FACS.data-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @title description of function scale.FACS.data
#' @export 
setGeneric('scale.FACS.data', ## Name
	function (dataObj) { ## Argumente der generischen Funktion
		standardGeneric('scale.FACS.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('scale.FACS.data', signature = c ('Tool_Pipe'),
	definition = function (dataObj) {
	dataObj <- as.matrix( dataObj )
	rown <- rownames( dataObj )
	dataObj <-  apply( dataObj,2, as.numeric)
	
#	neg <- which ( dataObj < 0)
#	pos <- which ( dataObj > 0) # == 0 stay 0
#	dataObj[neg] <- -log10( - dataObj[neg] )
#	dataObj[pos] <- log10( dataObj[pos] )
	dataObj[which(dataObj < 1)] <- 1
	dataObj <- log10(dataObj)
	rownames(dataObj) <- rown
	dataObj
} )
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
	function ( fname, previous=NULL, keyword='Array' ) { ## Argumente der generischen Funktion
		standardGeneric('file.position') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
#' @name read.PCR.set
#' @aliases read.PCR.set,Rscexv-method
#' @rdname read.PCR.set-methods
#' @docType methods
#' @description 
#' @param fnames  TEXT MISSING
#' @param use_pass_fail  TEXT MISSING
#' @title description of function read.PCR.set
#' @export 
setGeneric('read.PCR.set', ## Name
	function (fnames, use_pass_fail) { ## Argumente der generischen Funktion
		standardGeneric('read.PCR.set') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('read.PCR.set', signature = c ('Tool_Pipe'),
	definition = function (fnames, use_pass_fail) {
	
	etab <-NULL
	order <- NULL
	if ( length(fnames) > 0) {
		for(i in 1:length(fnames)){
			if ( ! fnames[i] == '../---' ){
				ttab <- read.PCR(fnames[i], use_pass_fail)
				ttab <- ttab[,order(colnames(ttab))]
				if ( length( grep( "P\\d+$", rownames(ttab))) == 0 ) {
					rownames(ttab) <- paste(rownames(ttab),".P",i-1,sep="")
				}
				## check whether the gene names are axactly the same
				if ( is.null(etab)){
					etab <-ttab
					order <- c(order, rep(i, nrow(ttab) ) )
				}
				else {
					if ( ! identical ( colnames(etab), colnames(ttab))) {
						system ( paste('echo "file', fnames[i],'not readable due to Gene Symbol mismatches!" >> R_file_read_error.txt') )
					}
					else {
						etab <- rbind(etab,ttab)
						order <- c(order, rep(i, nrow(ttab) ) )
					}
				}
			}
		}
	}
	list ( data = etab, order = order )
} )
#' @name match.sample.names
#' @aliases match.sample.names,Rscexv-method
#' @rdname match.sample.names-methods
#' @docType methods
#' @description 
#' @param PCR.s  TEXT MISSING
#' @param FACS.s  TEXT MISSING
#' @title description of function match.sample.names
#' @export 
setGeneric('match.sample.names', ## Name
	function ( PCR.s, FACS.s ) { ## Argumente der generischen Funktion
		standardGeneric('match.sample.names') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('match.sample.names', signature = c ('Tool_Pipe'),
	definition = function ( PCR.s, FACS.s ) {
#	system ( paste ('echo "PCR samples c(\'',paste(PCR.s,collapse="', '" ),'\')" >> R_file_read_error.txt', collapse=' ' ) )
#	system ( paste ('echo "FACS samples c(\'',paste(FACS.s,collapse="', '" ),'\')" >> R_file_read_error.txt', collapse=' ' ) )
	
	FACS.s <- str_replace_all( FACS.s, '.P[0-9]+', '' )
	
#	system ( paste ('echo "FACS replaced samples',paste(FACS.s,collapse=' ' ),'" >> R_file_read_error.txt', collapse=' ' ) )
	
	pos.in.FACS <- vector( length=length(PCR.s) )
	for ( i in 1:length(PCR.s) ) {
		pos <- which ((is.na(match (PCR.s, FACS.s[i] )) == FALSE) == T)
		if ( length(pos) == 1 ) {
#			print ( paste ( FACS.s[i],i, 'matches to (which)',pos, PCR.s[pos] ) )
			pos.in.FACS[pos] = i
			next
		}
		pos <- grep( paste(FACS.s[i], '$',sep=""), PCR.s, perl=T )
		if ( length(pos) == 1 ) {
#			print ( paste ( FACS.s[i],i, 'grep $',pos, PCR.s[pos] ) )
			pos.in.FACS[pos] = i
			next
		}
		pos <- grep( paste(FACS.s[i], '[\\s_\\-\\.]',sep=""), PCR.s, perl=T )
		if ( length(pos) == 1 ) {
#			print ( paste ( FACS.s[i], i,'grep with space after',pos, PCR.s[pos] ) )
			pos.in.FACS[pos] = i
			next
		}
		pos <- grep( FACS.s[i], PCR.s, perl=T )
		if ( length(pos) == 1 ) {
#			print ( paste ( FACS.s[i],i, 'grep string',pos, PCR.s[pos] ) )
			pos.in.FACS[pos] = i
			next
		}
	}
	system ( paste ('echo "unmatched PCR in the FACS data',paste( PCR.s[which(pos.in.FACS==0)], collapse=', '),'" >> R_file_read_error.txt', collapse=' ' ) )
	pos.in.FACS
} )
#' @name check.dataObj
#' @aliases check.dataObj,Rscexv-method
#' @rdname check.dataObj-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @title description of function check.dataObj
#' @export 
setGeneric('check.dataObj', ## Name
	function ( dataObj ) { ## Argumente der generischen Funktion
		standardGeneric('check.dataObj') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('check.dataObj', signature = c ('Tool_Pipe'),
	definition = function ( dataObj ) {
	colnames( dataObj$PCR$data ) = make.names( colnames( dataObj$PCR$data ),unique=T)
	if ( ! is.null(dataObj$FACS ) ){
		colnames( dataObj$FACS$data ) = make.names( colnames( dataObj$FACS$data ),unique=T)
		reject = 0
		if ( length(unique(dataObj$PCR$order)) != length(unique(dataObj$FACS$order))){
			reject = 1
		}else {
			new.facs <- matrix( nrow = nrow(dataObj$PCR$data), ncol=ncol(dataObj$FACS$data), 0 )
			colnames( new.facs ) = colnames ( dataObj$FACS$data )
			for ( i in 1:length(unique(dataObj$PCR$order))){
				## all samples existant in PCR but not in FACS (putative control wells) get a random negative expression
				f1 <- which( dataObj$PCR$order == i )
				F <- which( dataObj$FACS$order == i )
				map <- match.sample.names ( rownames(dataObj$PCR$data)[f1], rownames(dataObj$FACS$data)[F] )
				new.facs[f1[which(map > 0) ],] <- dataObj$FACS$data[F[ map[which(map > 0) ]], ]
				missing = which( map == 0 )
				if ( length ( missing) > 0) {
					all.f1 <- rownames(dataObj$PCR$data)[f1]
					system ( paste ('echo "missing cell in FACS data',paste(all.f1[missing],collapse=', '),'" >> R_file_read_warn.txt' ) )
					for ( a in 1:length(missing) ){
						#print ( paste("Problematic a =", a, "?", missing[a]) )
						if ( length(missing[a]) == 1 ){
							new.facs[f1[ missing[a] ], ] = log10(abs(rnorm ( ncol(dataObj$FACS$data), mean = 5 , sd = 1 )))
						}
					}
				}
			}
			rownames(new.facs) <- rownames(dataObj$PCR$data)
			dataObj$FACS$olddata <- dataObj$FACS$data
			## there might be empty lines in the FACS data!
			missing <- which (apply( new.facs,1,sd) == 0)
			if ( length ( missing ) > 0) {
				for ( a in 1:length(missing) ){
					new.facs[missing[a],] = log10(abs(rnorm ( ncol(dataObj$FACS$data), mean = 5 , sd = 1 )))
				}
			}
			dataObj$FACS$data <- new.facs
		}
	}
	dataObj
} )
#' @name plug.999
#' @aliases plug.999,Rscexv-method
#' @rdname plug.999-methods
#' @docType methods
#' @description 
#' @param tab  TEXT MISSING
#' @param plug  TEXT MISSING default=45
#' @title description of function plug.999
#' @export 
setGeneric('plug.999', ## Name
	function (tab,plug=45) { ## Argumente der generischen Funktion
		standardGeneric('plug.999') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('plug.999', signature = c ('Tool_Pipe'),
	definition = function (tab,plug=45) {
	
	ind <- which(tab==999)
	tab[ind] <- plug
	tab
	
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
	function (tab) { ## Argumente der generischen Funktion
		standardGeneric('hist.and.report.range') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
	function (dataObj,ref.nms,thresh,howm) { ## Argumente der generischen Funktion
		standardGeneric('filter.on.controls') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
#' @name remove.samples
#' @aliases remove.samples,Rscexv-method
#' @rdname remove.samples-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param ids  TEXT MISSING
#' @title description of function remove.samples
#' @export 
setGeneric('remove.samples', ## Name
	function ( dataObj, ids ) { ## Argumente der generischen Funktion
		standardGeneric('remove.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('remove.samples', signature = c ('Tool_Pipe'),
	definition = function ( dataObj, ids ) {
	
	if ( length(ids) > 0 ){
		write ( rownames(dataObj$PCR)[ids], file="./filtered_samples.txt",ncolumn=1, append=T )
		dataObj$PCR <- dataObj$PCR[-ids,]
		if ( ! is.null( dataObj$FACS) ){
			dataObj$FACS <- dataObj$FACS[-ids,]
		}
		dataObj$ArrayID <- dataObj$ArrayID[-ids]
		## work on all other possible objects
		if ( with( dataObj, exists('z') ) ) {
			dataObj$z$PCR <- dataObj$z$PCR[-ids,]
		}
	}
	else {
		print ( "No samples to filter out!" )
	}
	dataObj	
} )
#' @name remove.FACS.genes
#' @aliases remove.FACS.genes,Rscexv-method
#' @rdname remove.FACS.genes-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param ids  TEXT MISSING
#' @title description of function remove.FACS.genes
#' @export 
setGeneric('remove.FACS.genes', ## Name
	function ( dataObj, ids ) { ## Argumente der generischen Funktion
		standardGeneric('remove.FACS.genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('remove.FACS.genes', signature = c ('Tool_Pipe'),
	definition = function ( dataObj, ids ) {
	if ( length(ids) > 0 ){
		write ( colnames(dataObj$FACS)[ids], file="./filtered_genes.txt",ncolumn=1, append=T )
		dataObj$FACS <- dataObj$FACS[,-ids]
		if ( with(dataObj, exists('z'))) {
			dataObj$z$FACS <- dataObj$z$FACS[,-ids]
		}
	}
	else {
		print ( "No genes to filter out!" )
	}
	dataObj	
} )
#' @name remove.genes
#' @aliases remove.genes,Rscexv-method
#' @rdname remove.genes-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param ids  TEXT MISSING
#' @title description of function remove.genes
#' @export 
setGeneric('remove.genes', ## Name
	function ( dataObj, ids ) { ## Argumente der generischen Funktion
		standardGeneric('remove.genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('remove.genes', signature = c ('Tool_Pipe'),
	definition = function ( dataObj, ids ) {
	if ( length(ids) > 0 ){
		write ( colnames(dataObj$PCR)[ids], file="./filtered_genes.txt",ncolumn=1, append=T )
		dataObj$PCR <- dataObj$PCR[,-ids]
		if ( with(dataObj, exists('z'))) {
			dataObj$z$PCR <- dataObj$z$PCR[,-ids]
		}
	}
	else {
		print ( "No genes to filter out!" )
	}
	dataObj	
} )
#' @name filter.on.controls.no.inv
#' @aliases filter.on.controls.no.inv,Rscexv-method
#' @rdname filter.on.controls.no.inv-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param ref.nms  TEXT MISSING
#' @param thresh  TEXT MISSING
#' @param howm  TEXT MISSING
#' @title description of function filter.on.controls.no.inv
#' @export 
setGeneric('filter.on.controls.no.inv', ## Name
	function (dataObj,ref.nms,thresh,howm) { ## Argumente der generischen Funktion
		standardGeneric('filter.on.controls.no.inv') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('filter.on.controls.no.inv', signature = c ('Tool_Pipe'),
	definition = function (dataObj,ref.nms,thresh,howm) {
	
	conts <- dataObj$PCR[,ref.nms]
	conts[which(conts < thresh)] <- 0
	conts[conts>0] <- 1
	if ( is.null(dim(conts)) ){
		nums <- as.vector(conts)
	}
	else {
		nums <- apply(conts,1,sum)
	}
	
	dataObj <- remove.samples(dataObj, which(nums>howm) )
	dataObj
} )
#' @name sd.filter
#' @aliases sd.filter,Rscexv-method
#' @rdname sd.filter-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @title description of function sd.filter
#' @export 
setGeneric('sd.filter', ## Name
	function (dataObj) { ## Argumente der generischen Funktion
		standardGeneric('sd.filter') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('sd.filter', signature = c ('Tool_Pipe'),
	definition = function (dataObj) {
	
	sds1 <- NULL
	sds2 <- NULL
	if ( with(dataObj, exists('z'))) {
		sds1 <- apply(dataObj$z$PCR,1,sd)
		sds2 <- apply(dataObj$z$PCR,2,sd)
	}
	else {
		sds1 <- apply(dataObj$PCR,1,sd)
		sds2 <- apply(dataObj$PCR,2,sd)
	}
	cuto1 <- which(sds1==0)
	cuto2 <- which(sds2==0)
	
	write ( colnames(dataObj$PCR)[cuto2], file="./filtered_genes.txt",ncolumn=1, append=T )
	
	if(length(cuto1) >0){
		dataObj <- remove.samples(dataObj, cuto1 )
	}
	
	if(length(cuto2) >0){
		dataObj <- remove.genes(dataObj, cuto2 )
	}
	
	dataObj
	
} )
#' @name rank.normalize
#' @aliases rank.normalize,Rscexv-method
#' @rdname rank.normalize-methods
#' @docType methods
#' @description 
#' @param ap  TEXT MISSING
#' @title description of function rank.normalize
#' @export 
setGeneric('rank.normalize', ## Name
	function (ap) { ## Argumente der generischen Funktion
		standardGeneric('rank.normalize') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('rank.normalize', signature = c ('Tool_Pipe'),
	definition = function (ap) {
	o.ap <- apply(ap,2,order)
	ap.o <- NULLtab <- 
			for (i in 1:ncol(o.ap))
				ap.o <- cbind(ap.o,ap[o.ap[,i],i])
	m.ap <- apply(ap.o,1,median)
	ap.n <- array(0,dim(ap))
	for (i in 1:ncol(ap.n))
		ap.n[o.ap[,i],i] <- m.ap
	ap.n
} )
#' @name my.median
#' @aliases my.median,Rscexv-method
#' @rdname my.median-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @param max.cyc  TEXT MISSING
#' @title description of function my.median
#' @export 
setGeneric('my.median', ## Name
	function (x, max.cyc)  { ## Argumente der generischen Funktion
		standardGeneric('my.median') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('my.median', signature = c ('Tool_Pipe'),
	definition = function (x, max.cyc)  {
	median( x[which( x != max.cyc )] ) 
} )
#' @name norm.PCR
#' @aliases norm.PCR,Rscexv-method
#' @rdname norm.PCR-methods
#' @docType methods
#' @description 
#' @param tab  TEXT MISSING
#' @param meth  TEXT MISSING default=c("none"
#' @param "meancontrolgenes"  TEXT MISSING default=c("none"
#' @param "maxexpression"  TEXT MISSING
#' @param "medianexpression"  TEXT MISSING
#' @param "quantile")  TEXT MISSING
#' @param ctrl  TEXT MISSING default=NA
#' @param max.cyc  TEXT MISSING default=NA
#' @title description of function norm.PCR
#' @export 
setGeneric('norm.PCR', ## Name
	function (tab,meth=c("none","mean control genes","max expression","median expression","quantile"),ctrl=NA,max.cyc) { ## Argumente der generischen Funktion
		standardGeneric('norm.PCR') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('norm.PCR', signature = c ('Tool_Pipe'),
	definition = function (tab,meth=c("none","mean control genes","max expression","median expression","quantile"),ctrl=NA,max.cyc) {
	
	tab.new <- NULL
	tab.na <- which( tab == max.cyc )
	
	if(meth=="none"){
		tab.new <- tab
	}
	else if(meth=="mean control genes"){
		if(length(ctrl)>0){
			if ( length(ctrl)>1 ) {
				mean.ctrl <- apply(tab[,ctrl],1,mean)
			}else {
				mean.ctrl <- as.vector(tab[,ctrl])
			}
			for(i in 1:nrow(tab)){
				tab.new <- rbind(tab.new,(tab[i,]-mean.ctrl[i]))
			} 
		}
	}
	else if (meth== "max expression" ){
		max.expr <- apply( tab,1,min )
		for(i in 1:nrow(tab)){
			tab.new <- rbind(tab.new,(tab[i,]-max.expr[i]))
		} 
	} 
	else if (meth== "median expression" ){
		median.expr <- apply( tab,1,my.median, max.cyc  )
		for(i in 1:nrow(tab)){
			tab.new <- rbind(tab.new,(tab[i,]-median.expr[i]))
		} 
	} 
	else if(meth=="quantile"){
		tab.new <- rank.normalize(tab)
		colnames(tab.new) <- colnames(tab)
	}
	else {
		stop(paste ("norm method",meth, "is not implemented!" ) )
	}
	if ( is.null(rownames(tab.new))){
		rownames(tab.new) <- rownames(tab)
	}
	
	tab.ret <- max.cyc-tab.new
	tab.ret[tab.na] <- 0
	tab.ret
	
} )
#' @name z.score.PCR.mad
#' @aliases z.score.PCR.mad,Rscexv-method
#' @rdname z.score.PCR.mad-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @title description of function z.score.PCR.mad
#' @export 
setGeneric('z.score.PCR.mad', ## Name
	function (dataObj) { ## Argumente der generischen Funktion
		standardGeneric('z.score.PCR.mad') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('z.score.PCR.mad', signature = c ('Tool_Pipe'),
	definition = function (dataObj) {
	
	#tab.new <- NULL
	#dataObj$ArrayID
	arrays <- max(dataObj$ArrayID)
	
	rem.inds <- NULL
	
	tab.pre <-NULL
	
	for (j in 1:arrays ){
		
		mads <- NULL
		meds <- NULL
		
		tdat <- dataObj$PCR[which(dataObj$ArrayID==j),]
		
		tab.new <- NULL
		
		for(i in 1:ncol(tdat)){
			vec <- tdat[which(tdat[,i]!=0),i]
			mads <- c(mads,mad(vec))
			meds <- c(meds,median(vec))
		}
		
		for(i in 1:ncol(tdat)){
			#print ( paste ( "tdat[,i]-", meds[i],")/( 1.48* ",mads[i],")"))
			new.v <- (tdat[,i]-meds[i])/(1.48*mads[i])
			#print ( paste (i, colnames(tdat)[i], range(new.v) ) ) 
			if(all(is.na(range(new.v))==T)){
				rem.inds <- c(rem.inds,i)
				#print (paste( i, colnames(tdat)[i],'Kicked in array', j ) )
				plug.ind <- 1:length(tdat[,i])
				## NA's have been created for all values if (all(is.na(range(new.v))==T))
				## => we need to set even those to -20 that might have gotten a value
				## far from optimal, but hopefully not problematic
			}
			else {
				plug.ind <- which(tdat[,i] == 0)
				#plug.ind <- which(new.v==(0-meds[i])/(1.48*mads[i]))
			}
			
			new.v[plug.ind] <- -20
			tab.new <- cbind(tab.new,new.v)
			
		}
		
		tab.pre <- rbind(tab.pre,tab.new)
		
	}
	
	rem.ind.fin <- as.numeric( names(table(rem.inds))[which(table(rem.inds)==arrays)]) ### might fall
	#rem.ind.fin <- as.numeric( names(table(rem.inds)))
	#print ( table(rem.inds) )
	#print (rem.ind.fin)
	
	write ( colnames(dataObj$PCR)[rem.ind.fin],file="./filtered_genes.txt",ncolumn=1, append=T )
	
	dataObj$z$kicked.z.score <- rem.ind.fin
	if ( length(rem.ind.fin) > 0 ){
		dataObj$PCR <- dataObj$PCR[,-rem.ind.fin]
		tab.pre <- tab.pre[,-rem.ind.fin]
	}
	
	#print ( paste ( sum(is.na(tab.pre)), "NA's produced" ) ) 
	
	colnames(tab.pre) <- colnames(dataObj$PCR) 
	
	#print ( apply(tab.pre,2,sd))
	
	dataObj$z$PCR <- tab.pre
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
	function ( dataObj ) { ## Argumente der generischen Funktion
		standardGeneric('auto_problem') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
	function (tab) { ## Argumente der generischen Funktion
		standardGeneric('clus.PCR') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
	function (tab,k) { ## Argumente der generischen Funktion
		standardGeneric('clus.z.PCR') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
#' @name plot.histograms
#' @aliases plot.histograms,Rscexv-method
#' @rdname plot.histograms-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param cuts  TEXT MISSING default=vector('list'
#' @param 1)  TEXT MISSING
#' @title description of function plot.histograms
#' @export 
setGeneric('plot.histograms', ## Name
	function ( dataObj, cuts=vector('list',1) ) { ## Argumente der generischen Funktion
		standardGeneric('plot.histograms') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('plot.histograms', signature = c ('Tool_Pipe'),
	definition = function ( dataObj, cuts=vector('list',1) ) {
	ma <- t(dataObj$PCR)
	if ( ! is.null(dataObj$FACS)){
		ma <- rbind( ma, t( dataObj$FACS) )
	}
	n <- rownames(ma)
	arrays <- max(dataObj$ArrayID)
	cols <- rainbow(arrays)
	n.cuts <- names(cuts)
	for ( i in 1:nrow(ma) ) {
		png( file=paste(n[i],'.png',sep=''),width=800, height=800 )
		temp <- vector('list',arrays)
		m <- NULL
		for (a in 1:arrays ) {
			temp[[a]] <- density(ma[i,which(dataObj$ArrayID == a )])
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
	
} )
#' @name kick.expressed.negContr.samples
#' @aliases kick.expressed.negContr.samples,Rscexv-method
#' @rdname kick.expressed.negContr.samples-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param negContrGenes  TEXT MISSING default=NULL
#' @title description of function kick.expressed.negContr.samples
#' @export 
setGeneric('kick.expressed.negContr.samples', ## Name
	function ( dataObj, negContrGenes=NULL ) { ## Argumente der generischen Funktion
		standardGeneric('kick.expressed.negContr.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('kick.expressed.negContr.samples', signature = c ('Tool_Pipe'),
	definition = function ( dataObj, negContrGenes=NULL ) {
	if ( ! is.null(negContrGenes) ){
		rem.samples <- NULL
		for ( i in match (  negContrGenes, colnames( dataObj$PCR) ) ){
			#browser()
			rem.samples <- c( rem.samples, which(t(dataObj$PCR[,i]) < 999 ) )
		}
		rem.samples <- unique( rem.samples )
		dataObj <- remove.samples(dataObj, rem.samples )
	}
	dataObj
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
	function (dataObj, old.names, new.names ) { ## Argumente der generischen Funktion
		standardGeneric('update.FACS.rownames') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
	function ( x, keyword = 'Array' ) { ## Argumente der generischen Funktion
		standardGeneric('files.sorted') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('files.sorted', signature = c ('Tool_Pipe'),
	definition = function ( x, keyword = 'Array' ) {
	x[order( regmatches(x, regexpr( paste(keyword,"\\d+",sep=''),x ,perl=T)) )]
} )
#' @name createDataObj
#' @aliases createDataObj,Rscexv-method
#' @rdname createDataObj-methods
#' @docType methods
#' @description 
#' @param PCR  TEXT MISSING default=NULL
#' @param FACS  TEXT MISSING default=NULL
#' @param max.value  TEXT MISSING default=40
#' @param createDataObj  TEXT MISSING
#' @title description of function createDataObj
#' @export 
setGeneric('createDataObj', ## Name
	function ( PCR=NULL,  FACS=NULL, max.value=40, createDataObj { ## Argumente der generischen Funktion
		standardGeneric('createDataObj') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('createDataObj', signature = c ('Tool_Pipe'),
	definition = function ( PCR=NULL,  FACS=NULL, max.value=40, createDataObj {
		ref.genes=NULL, max.ct=25, max.control=0,  norm.function='none', negContrGenes=NULL, use_pass_fail = T, ... ){
	
	data <- vector ('list' ,2 )
	names(data) <- c( 'PCR', 'FACS')
	data$FACS <- NULL
	
	system ( 'rm filtered*' )
	
	try ( system ( 'rm R_file_read_error.txt' ) )
	
	if ( is.null(PCR) ){
		system ( 'echo "You need to upload at least one PCR data set in order to start the analysis" > R.error')
		stop("You need to upload at least one PCR data set in order to start the analysis")
	}
	
	data$PCR <- read.PCR.set(PCR, use_pass_fail)
	colnames(data$PCR$data) <- str_replace_all( colnames(data$PCR$data), '/', '_' )
	
	if ( ! is.null(FACS) ){
		data$FACS <- read.FACS.set(FACS)
		## bloody hack to not create problems!
		data$FACS$data[which(is.na(data$FACS$data))] <- 0
		## black magick with the FACS gene names
		colnames(data$FACS$data) <- str_replace_all( colnames(data$FACS$data), '^P[0-9]+.', '' )
		colnames(data$FACS$data) <- str_replace_all( colnames(data$FACS$data), '.Min$', '' )
		if ( is.null(data$FACS$data ) ){
			system ( 'echo "FACS sample names do not match with the PCR sample names - please fix that!" >> R.error')
		}
	}
	data <- check.dataObj(data)
	data$ArrayID <- data$PCR$order
	data$PCR <- data$PCR$data
	if ( ! is.null(FACS) ){
		data$FACS <- data$FACS$data
	}
	
	data <- kick.expressed.negContr.samples(data, negContrGenes )
	
	data$PCR <- plug.999(data$PCR, max.value ) ## does nothing for pre-processed data
	
	if ( all ( data$PCR == 40 ) ){
		system ( 'echo "Please check your filter settings - all samples have been removed from the analysis!" >> R_file_read_error.txt' )
	}
	
	if ( ! is.null(ref.genes)){
		data <- filter.on.controls.no.inv(data,ref.genes,max.ct,max.control)	
	}
	
	## export the unfiltered_not_modified PCR data for publication
	write.table( data$PCR, file="../PCR_data_4_publication.xls", sep='\t' )
	
	data.filtered <- sd.filter( data )
	
	#plot.heatmap( list( data = t(data.filtered$PCR), genes=colnames(data.filtered$PCR)), 'SD_filtered_not_norm', title='SD filtered RAW data', width=12,height=6,Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},distfun = function(x){ 1- cor(t(x),method='spearman')} )
	
	plot.histograms( data.filtered )
	
	data.filtered$PCR <- norm.PCR(data.filtered$PCR,norm.function,max.cyc=40, ctrl=ref.genes )
	#plot.heatmap( list( data = t(data.filtered$PCR), genes=colnames(data.filtered$PCR)), 'Contr_filtered_inverted_norm', title='SD filtered inverted data', width=12,height=6,Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},distfun = function(x){ 1- cor(t(x),method='spearman')} )
	write.table( data.filtered$z$PCR, file="../PCR_data_normalized_4_publication.xls", sep='\t' )

	data.filtered <- z.score.PCR.mad(data.filtered)
	#data.filtered$z$PCR <- data.filtered$PCR
	
	arrays <- arrays <- max(data.filtered$ArrayID)
	cols <- rainbow(arrays)
	cmethod <- 'ward'
	try(PCR.heatmap( 
					list( data = t(data.filtered$PCR), genes=colnames(data.filtered$PCR)),
					'Contr_filtered_inverted_norm', 
					title='SD filtered and normlizied data', 
					width=12,height=6,
					Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},
					distfun = function(x){ as.dist(1- cor(t(x),method='spearman'))}, 
					margins = c(0,15), 
					lwid = c( 0.05, 0.45),
					ColSideColors=cols[data.filtered$ArrayID] )
	)
	
	try(PCR.heatmap( 
					list( data = t(data.filtered$z$PCR), genes=colnames(data.filtered$z$PCR)),
					'Contr_filtered_inverted_norm_Zscore', 
					title='SD filtered normalized data and Z scored', 
					width=12,height=6,
					Colv=F,hclustfun = function(c){hclust( c, method=cmethod)},
					distfun = function(x){as.dist( 1- cor(t(x),method='spearman'))}, 
					margins = c(0,15), 
					lwid = c( 0.05, 0.45),
					ColSideColors=cols[data.filtered$ArrayID] )
	)
	
	devSVG( file="boxplot_filtered_samples.svg", height =8, width =18 )
	par(mar=c(12,5,2,2))
	tmp <- t(data.filtered$PCR)
	tmp [ which(tmp == 0 ) ] <- NA
	boxplot( tmp , las=2,cex.lab=0.5, main ="Normalized expression in all used Samples")
	for ( i in 1:arrays){
		da <- tmp[,which(data.filtered$ArrayID == i ) ]
		abline( h=median(da[which(! is.na(da) ) ]), col = cols[i], lwd=3 )
	}
	dev.off()
	
	devSVG( file="boxplot_filtered_zscored_samples.svg", height =8, width =18 )
	par(mar=c(12,5,2,2))
	tmp <- t(data.filtered$z$PCR)
	tmp [ which(tmp == -20 ) ] <- NA
	boxplot( tmp , las=2,cex.lab=0.5, main="Normalized and z-scored")
	for ( i in 1:arrays){
		da <- tmp[,which(data.filtered$ArrayID == i ) ]
		abline( h=median(da[which(! is.na(da) ) ]), col = cols[i], lwd=3 )
	}
	dev.off()
	
	data.filtered
} )
