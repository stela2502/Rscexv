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
		function (fnames, use_pass_fail=T) { 
			standardGeneric('read.PCR.set')
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



