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
			dataObj <-  apply( dataObj,2, function(x) {x = as.vector(x); x <- str_replace_all(x,',','.');as.numeric(x) })
			t <- as.matrix( dataObj )
			if ( all(floor(t) == t, na.rm = TRUE)) {
				## this data is not in log space - log it
				dataObj[which(dataObj < 1)] <- 1
				dataObj <- log10(dataObj)
			}else {
				warn ("FACS data was detected as logged - no log applied")
			}
			rownames(dataObj) <- rown
			dataObj
		} 
)

