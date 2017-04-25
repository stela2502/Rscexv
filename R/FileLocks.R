#' @name set.lock
#' @aliases set.lock,Rscexv-method
#' @rdname set.lock-methods
#' @docType methods
#' @description set a lock for a file (threading)
#' @param filename  the locked file
#' @title description of function set_lock
#' @export 
setGeneric('set.lock', ## Name
		function ( filename ) { ## Argumente der generischen Funktion
			standardGeneric('set.lock') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('set.lock', signature = c ('character'),
		definition = function ( filename ) {
			system ( paste('touch ',filename,'.lock', sep='') )
		} )
#' @name release.lock
#' @aliases release.lock,Rscexv-method
#' @rdname release.lock-methods
#' @docType methods
#' @description releases the lock of a file (threading)
#' @param filename  the locked file
#' @title description of function release_lock
#' @export 
setGeneric('release.lock', ## Name
		function ( filename ) { ## Argumente der generischen Funktion
			standardGeneric('release.lock') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('release.lock', signature = c ('character'),
		definition = function ( filename ) {
			system ( paste('rm -f ',filename,'.lock', sep='') )
		} )
#' @name locked
#' @aliases locked,Rscexv-method
#' @rdname locked-methods
#' @docType methods
#' @description simple check for a file lock (threading)
#' @param filename  lock this file
#' @title description of function locked
#' @export 
setGeneric('locked', ## Name
		function ( filename ) { ## Argumente der generischen Funktion
			standardGeneric('locked') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('locked', signature = c ('character'),
		definition = function ( filename ) {
			ret = TRUE
			if (file.exists( filename )){
				ret <- file.exists( paste(filename,'.lock', sep='') )
			}
			ret
		} )