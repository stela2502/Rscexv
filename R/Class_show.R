#' @name show
#' @aliases show,Rscexv-method
#' @rdname show-methods
#' @docType methods
#' @description  print the Rscexv
#' @param object the Rscexv object
#' @return nothing
#' @title description of function show
#' @export 
setMethod('show', signature(object='Rscexv') ,
		definition = function (object) {
			cat (paste("An object of class", class(object)),"\n" )
			cat("named ",object@name,"\n")
			cat (paste( 'with',nrow(object@data),'samples and', ncol(object@data),' genes.'),"\n")
			if ( object@wFACS ) {
				cat ( paste("The object also contains FACS for",ncol(object@facs),'fluorocromes',sep=' '),"\n")
			}
			cat (paste("Annotation datasets (",paste(dim(object@annotation),collapse=','),"): '",paste( colnames(object@annotation ), collapse="', '"),"'  ",sep='' ),"\n")
			cat (paste("Sample annotation (",paste(dim(object@samples),collapse=','),"): '",paste( colnames(object@samples ), collapse="', '"),"'  ",sep='' ),"\n")
			cat ( paste ( "A total of",length(object@usedObj),"Other objects have been added to this Rscexv object:\n", paste( collapse=', ', names(object@usedObj))),"\n")
		}
)


