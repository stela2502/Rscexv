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
			tab <- as.matrix(etab)
			tab[which(is.na(tab))] <- 0
			etab <- as.data.frame( tab )
			list ( data = scale.FACS.data(etab), order = order )
		}
)


