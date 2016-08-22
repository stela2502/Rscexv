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


