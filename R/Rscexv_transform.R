#' @name kick.expressed.negContr.samples
#' @aliases kick.expressed.negContr.samples,Rscexv-method
#' @rdname kick.expressed.negContr.samples-methods
#' @docType methods
#' @description IF you are sure, that a list of genes must not be expressed in cells of interest you can fore to exclude these ceells using this function.
#' @description But this is a quite bad absolutely not recommended workflow!
#' @description This function drops samples with expression different from 999 - therefore you must call this function before removing the 999 values.
#' @param dataObj the Rscexv object
#' @param negContrGenes which genes should be checked default=NULL
#' @title description of function kick.expressed.negContr.samples
#' @export 
setGeneric('kick.expressed.negContr.samples', ## Name
		function ( dataObj, negContrGenes=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('kick.expressed.negContr.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('kick.expressed.negContr.samples', signature = c ('Rscexv'),
		definition = function ( dataObj, negContrGenes=NULL ) {
			if ( ! is.null(negContrGenes) ){
				rem.samples <- NULL
				rem.samples <- unique( 
						apply( dataObj@data[,negContrGenes], 2, function(x) { 
									which(t(dataObj@data[,i]) < 999 )
								} ) 
				)
				dataObj <- remove.samples(dataObj, rem.samples )
			}
			dataObj }
)

#' @name remove.samples
#' @aliases remove.samples,Rscexv-method
#' @rdname remove.samples-methods
#' @docType methods
#' @description Remove samples by id. 
#' @param dataObj the Rscexv object
#' @param ids which samples (ids!) to remove
#' @title description of function remove.samples
#' @export 
setGeneric('remove.samples', ## Name
		function ( dataObj, ids ) { ## Argumente der generischen Funktion
			standardGeneric('remove.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('remove.samples', signature = c ('Rscexv'),
		definition = function ( dataObj, ids ) {
			
			if ( length(ids) > 0 ){
				write ( rownames(dataObj@data)[ids], file="./filtered_samples.txt",ncolumn=1, append=T )
				dataObj@data <- dataObj@data[-ids,]
				if ( dataObj@wFACS ){
					dataObj@facs <- dataObj@facs[-ids,]
				}
				dataObj@samples <- dataObj@samples[-ids,]
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
#' @description this function removes a list of FACS genes from the analysis
#' @param dataObj the Rscexv object
#' @param ids the gene ids to remove
#' @title description of function remove.FACS.genes
#' @export 
setGeneric('remove.FACS.genes', ## Name
		function ( dataObj, ids ) { ## Argumente der generischen Funktion
			standardGeneric('remove.FACS.genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('remove.FACS.genes', signature = c ('Rscexv'),
		definition = function ( dataObj, ids ) {
			if ( length(ids) > 0 ){
				write ( colnames(dataObj@facs)[ids], file="./filtered_genes.txt",ncolumn=1, append=T )
				dataObj@facs <- dataObj@facs[,-ids]
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
#' @description this function removes a list of PCR genes from the analysis
#' @param dataObj the Rscexv object
#' @param ids the gene ids to remove
#' @title description of function remove.genes
#' @export 
setGeneric('remove.genes', ## Name
		function ( dataObj, ids ) { ## Argumente der generischen Funktion
			standardGeneric('remove.genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('remove.genes', signature = c ('Rscexv'),
		definition = function ( dataObj, ids ) {
			if ( length(ids) > 0 ){
				write ( colnames(dataObj@data)[ids], file="./filtered_genes.txt",ncolumn=1, append=T )
				dataObj@data <- dataObj@data[,-ids]
				dataObj@annotation <- dataObj@annotation[-ids,]
			}
			else {
				print ( "No genes to filter out!" )
			}
			dataObj	
		} )

