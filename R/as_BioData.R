#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @description create a BioData object from a Rscexv object
#' @param dat the Rscexv object
#' @title description of function as
#' @export as_BioData
if ( ! isGeneric('as_BioData') ){ setGeneric('as_BioData', ## Name
			function ( dat ) { 
				standardGeneric('as_BioData')
			}
	)
}else {
	print ("Onload warn generic function 'as_BioData' already defined - no overloading here!")
}

setMethod('as_BioData', signature = c ('Rscexv'),
		definition = function ( dat ) {
			#cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = ""
			if ( ncol(dat@facs) != 0 ) {
				print ("The Rscexv FACS data is lost in the conversion process" )
			}
			## the BioData sores samples and genes in a transposed form!
			
			ok = which(lapply(colnames(dat@samples) , function(x) { all.equal( as.character(as.vector(dat@samples[,x])),make.names(rownames(dat@data))) == T } )==T)
			if ( length(ok) == 0) {
				dat@samples = cbind( cell.name = make.names(rownames(dat@data)), dat@samples)
				rownames(dat@data) = make.names(rownames(dat@data))
				namecol = 'cell.name'
			}
			else {
				namecol = colnames(dat@samples)[ok]
				namecol = namecol[1]
			}
			
			namerow = 'Gene.Name'
			if ( is.null(dat@annotation[,'Gene.Name']) ) {
				dat@annotation[,'Gene.Name'] = colnames(dat@data)
			}
			
			d <- data.frame(cbind( dat@annotation,t(dat@data) ))
			
			ret <- BioData$new( d, Samples=dat@samples, name= 'from.Rscexv', namecol= namecol, namerow=namerow, outpath='./' )
			ret$usedObj <- dat@usedObj
			ret
		}
)