#' @name check.dataObj
#' @aliases check.dataObj,Rscexv-method
#' @rdname check.dataObj-methods
#' @docType methods
#' @description checks the two data frames and creates a sample and annotation table
#' @param PCR the object from read.PCR.set
#' @param FACS the object from read.FACS.set
#' @title description of function check.dataObj
#' @export 
setGeneric('check.dataObj', ## Name
		function ( PCR, FACS ) { ## Argumente der generischen Funktion
			standardGeneric('check.dataObj') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('check.dataObj', signature = c ('list'),
		definition = function ( PCR, FACS ) {
			colnames( PCR$data ) = make.names( colnames( PCR$data ),unique=T)
			new.facs = NULL
			if ( ! is.null(FACS ) ){
				colnames( FACS$data ) = make.names( colnames( FACS$data ),unique=T)
				reject = 0
				if ( length(unique(PCR$order)) != length(unique(FACS$order))){
					reject = 1
				}else {
					new.facs <- matrix( nrow = nrow(PCR$data), ncol=ncol(FACS$data), 0 )
					colnames( new.facs ) = colnames ( FACS$data )
					for ( i in 1:length(unique(PCR$order))){
						## all samples existant in PCR but not in FACS (putative control wells) get a random negative expression
						f1 <- which( PCR$order == i )
						F <- which( FACS$order == i )
						map <- match.sample.names ( rownames(PCR$data)[f1], rownames(FACS$data)[F] )
						new.facs[f1[which(map > 0) ],] <- FACS$data[F[ map[which(map > 0) ]], ]
						missing = which( map == 0 )
						if ( length ( missing) > 0) {
							all.f1 <- rownames(PCR$data)[f1]
							system ( paste ('echo "missing cell in FACS data',paste(all.f1[missing],collapse=', '),'" >> R_file_read_warn.txt' ) )
							for ( a in 1:length(missing) ){
								#print ( paste("Problematic a =", a, "?", missing[a]) )
								if ( length(missing[a]) == 1 ){
									new.facs[f1[ missing[a] ], ] = log10(abs(rnorm ( ncol(FACS$data), mean = 5 , sd = 1 )))
								}
							}
						}
					}
					rownames(new.facs) <- rownames(PCR$data)
					## there might be empty lines in the FACS data!
					missing <- which (apply( new.facs,1,sd) == 0)
					if ( length ( missing ) > 0) {
						for ( a in 1:length(missing) ){
							new.facs[missing[a],] = log10(abs(rnorm ( ncol(FACS$data), mean = 5 , sd = 1 )))
						}
					}
					FACS$data <- new.facs
				}
			}
			samples <- data.frame( SampleName = rownames(PCR$data), ArrayID = PCR$order )
			list( PCR= PCR$data, FACS = new.facs, samples= samples, annotation=data.frame('Gene Name' = colnames(PCR$data) ) )
		}
)

