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
			list( PCR= PCR$data, FACS = new.facs, samples= samples, annotation=data.frame('GeneID' = colnames(PCR$data) ) )
		}
)

#' @name filter.on.controls.no.inv
#' @aliases filter.on.controls.no.inv,Rscexv-method
#' @rdname filter.on.controls.no.inv-methods
#' @docType methods
#' @description This function filteres out samples that have an expression of less than thresh in there ref.nms(gene names)
#' @param dataObj the Rscexv object
#' @param ref.nms the gene names used as positive control
#' @param thresh the minimum expression in each ref gene
#' @param howm drop if one sample has howm genes that fail the test
#' @title description of function filter.on.controls.no.inv
#' @export 
setGeneric('filter.on.controls.no.inv', ## Name
		function (dataObj,ref.nms,thresh,howm) { ## Argumente der generischen Funktion
			standardGeneric('filter.on.controls.no.inv') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('filter.on.controls.no.inv', signature = c ('Rscexv'),
		definition = function (dataObj,ref.nms,thresh,howm) {
			
			conts <- dataObj@data[,ref.nms]
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
		} 
)
