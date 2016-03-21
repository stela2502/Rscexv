#' @name force.unique.sample
#' @aliases force.unique.sample,Rscexv-method
#' @rdname force.unique.sample-methods
#' @docType methods
#' @description This function will just make sure, that all samples have uniaue sample names, even if they come from different arrays.
#' @param x  TEXT MISSING
#' @title description of function force.unique.sample
#' @export 
setGeneric('force.unique.sample', ## Name
	function ( x ) { ## Argumente der generischen Funktion
		standardGeneric('force.unique.sample') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('force.unique.sample', signature = c ('character'),
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

#' @name read.PCR
#' @aliases read.PCR,Rscexv-method
#' @rdname read.PCR-methods
#' @docType methods
#' @description reads one fluidigm data file
#' @param fname  the fliename to read
#' @param use_pass_fail exclude data from failed wells default=T
#' @title description of function read.PCR
#' @export 
setGeneric('read.PCR', ## Name
	function (fname,use_pass_fail=T) { ## Argumente der generischen Funktion
		standardGeneric('read.PCR') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('read.PCR', signature = c ('character'),
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
#' @description reads one fluidigm data file in heatmap format
#' @param fname  TEXT MISSING
#' @param use_pass_fail  TEXT MISSING default=T
#' @title description of function read.PCR.heatmap
#' @export 
setGeneric('read.PCR.heatmap', ## Name
		function (fname,use_pass_fail=T) { ## Argumente der generischen Funktion
			standardGeneric('read.PCR.heatmap') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('read.PCR.heatmap', signature = c ('character'),
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
		} 
)

		
#' @name read.FACS
#' @aliases read.FACS,Rscexv-method
#' @rdname read.FACS-methods
#' @docType methods
#' @description THis function reads one INDEX data file.
#' @param fname the file containg the INDEX data
#' @title description of function read.FACS
#' @export 
setGeneric('read.FACS', ## Name
		function (fname) { ## Argumente der generischen Funktion
			standardGeneric('read.FACS') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('read.FACS', signature = c ('character'),
		definition = function (fname) {
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
			data.frame(ftab)
		} 
)


#' @name read.FACS.set
#' @aliases read.FACS.set,Rscexv-method
#' @rdname read.FACS.set-methods
#' @docType methods
#' @description reads a whole set of FACS data files
#' @param fnames  TEXT MISSING
#' @title description of function read.FACS.set
#' @export 
setGeneric('read.FACS.set', ## Name
		function (fnames) { ## Argumente der generischen Funktion
			standardGeneric('read.FACS.set') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('read.FACS.set', signature = c ('character'),
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
			etab[which(is.na(etab))] <- 0
			list ( data = scale.FACS.data(etab), order = order )
		}
)


#' @name read.PCR.set
#' @aliases read.PCR.set,Rscexv-method
#' @rdname read.PCR.set-methods
#' @docType methods
#' @description reads a whole set of PCR data files
#' @param fnames  a list of file names
#' @param use_pass_fail  whether to remove values marked as failed in the files (default=T)
#' @title description of function read.PCR.set
#' @export 
setGeneric('read.PCR.set', ## Name
		function (fnames, use_pass_fail=T) { ## Argumente der generischen Funktion
			standardGeneric('read.PCR.set') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('read.PCR.set', signature = c ('character'),
		definition = function (fnames, use_pass_fail=T) {
			
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
		} 
)



