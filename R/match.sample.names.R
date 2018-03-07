#' @name match.sample.names
#' @aliases match.sample.names,Rscexv-method
#' @rdname match.sample.names-methods
#' @docType methods
#' @description This function makes sure, that both the FACS and the PCR data samples are matched in the object
#' @param PCR.s  the PCR sample names
#' @param FACS.s the FACS sample names
#' @title description of function match.sample.names
#' @export 
setGeneric('match.sample.names', ## Name
		function ( PCR.s, FACS.s ) { ## Argumente der generischen Funktion
			standardGeneric('match.sample.names') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('match.sample.names', signature = c ('character'),
		definition = function ( PCR.s, FACS.s ) {		
			FACS.s <- str_replace_all( FACS.s, '.P[0-9]+', '' )			
			pos.in.FACS <- vector( length=length(PCR.s) )
			for ( i in 1:length(PCR.s) ) {
				pos <- which ((is.na(match (PCR.s, FACS.s[i] )) == FALSE) == T)
				if ( length(pos) == 1 ) {
					pos.in.FACS[pos] = i
					next
				}
				pos <- grep( paste(FACS.s[i], '$',sep=""), PCR.s, perl=T )
				if ( length(pos) == 1 ) {
					pos.in.FACS[pos] = i
					next
				}
				pos <- grep( paste(FACS.s[i], '[\\s_\\-\\.]',sep=""), PCR.s, perl=T )
				if ( length(pos) == 1 ) {
					pos.in.FACS[pos] = i
					next
				}
				pos <- grep( FACS.s[i], PCR.s, perl=T )
				if ( length(pos) == 1 ) {
					pos.in.FACS[pos] = i
					next
				}
			}
			system ( paste ('echo "unmatched PCR in the FACS data',paste( PCR.s[which(pos.in.FACS==0)], collapse=', '),'" >> R_file_read_error.txt', collapse=' ' ) )
			pos.in.FACS
		} 
)


