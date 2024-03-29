#' @name read.FACS
#' @aliases read.FACS,Rscexv-method
#' @rdname read.FACS-methods
#' @docType methods
#' @description THis function reads one INDEX data file.
#' @param fname the file containg the INDEX data
#' @title description of function read.FACS
#' @export 
setGeneric('read.FACS', ## Name
		function (fname) { 
			standardGeneric('read.FACS')
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
				if (length(grep ('^\\.', colnames( ftab ))) > 0 ) {
					if ( length(length(grep ('^\\.', colnames( ftab )))) == ncol(ftab) ) {
						message( "FACS Data is not as expected - keeping all columns to be sure!" )
					}else {
						ftab <- ftab[,-grep( '^\\.' , colnames(ftab))]
					}
				}
			}
			data.frame(ftab)
		} 
)


