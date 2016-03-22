
#' @name plug.999
#' @aliases plug.999,Rscexv-method
#' @rdname plug.999-methods
#' @docType methods
#' @description this function replaces the 999 (failed) values by the plug value.
#' @param x the Rscexv object
#' @param plug the replacement value (default=45)
#' @title description of function plug.999
#' @export 
setGeneric('plug.999', ## Name
		function (x,plug=45) { ## Argumente der generischen Funktion
			standardGeneric('plug.999') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plug.999', signature = c ('Rscexv'),
		definition = function (x,plug=45) {
			ma <- as.matrix(x@data)
			ind <- which(ma==999)
			ma[ind] <- plug
			x@data <- data.frame(ma)
			x
		} 
)

#' @name scale.FACS.data
#' @aliases scale.FACS.data,Rscexv-method
#' @rdname scale.FACS.data-methods
#' @docType methods
#' @description This function applies a simple log10 onto all expression data.
#' @param dataObj  TEXT MISSING
#' @title description of function scale.FACS.data
#' @export 
setGeneric('scale.FACS.data', ## Name
		function (dataObj) { ## Argumente der generischen Funktion
			standardGeneric('scale.FACS.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('scale.FACS.data', signature = c ('data.frame'),
		definition = function (dataObj) {
			rown <- rownames( dataObj )
			dataObj <-  apply( dataObj,2, as.numeric)			
			dataObj[which(dataObj < 1)] <- 1
			dataObj <- log10(dataObj)
			rownames(dataObj) <- rown
			dataObj
		} 
)

#' @name sd.filter
#' @aliases sd.filter,Rscexv-method
#' @rdname sd.filter-methods
#' @docType methods
#' @description This function is removing cells and genes that do not show any variation in the data.
#' @param dataObj the Rscexv object
#' @title description of function sd.filter
#' @export 
setGeneric('sd.filter', ## Name
		function (dataObj) { ## Argumente der generischen Funktion
			standardGeneric('sd.filter') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('sd.filter', signature = c ('Rscexv'),
		definition = function (dataObj) {
			
			sds1 <- NULL
			sds2 <- NULL

			sds1 <- apply(dataObj@data,1,sd)
			sds2 <- apply(dataObj@data,2,sd)
			
			cuto1 <- which(sds1==0)
			cuto2 <- which(sds2==0)
			
			if(length(cuto1) >0){
				write ( colnames(dataObj@data)[cuto2], file="./filtered_samples.txt",ncolumn=1, append=T )
				dataObj <- remove.samples(dataObj, cuto1 )
			}
			
			if(length(cuto2) >0){
				write ( colnames(dataObj@data)[cuto2], file="./filtered_genes.txt",ncolumn=1, append=T )
				dataObj <- remove.genes(dataObj, cuto2 )
			}
			
			dataObj			
		} 
)

