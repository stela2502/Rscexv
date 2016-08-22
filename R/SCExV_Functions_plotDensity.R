#' @name plotDensity
#' @aliases plotDensity,Rscexv-method
#' @rdname plotDensity-methods
#' @docType methods
#' @description This function creates the density plot for the analysis page
#' @param obj the Rscexv object
#' @title description of function analyse.data
#' @export 
setGeneric('plotDensity', ## Name
		function ( obj ){	## Argumente der generischen Funktion
			standardGeneric('plotDensity') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plotDensity', signature = c ('Rscexv'),
		definition = function ( obj ){
			usable <- is.na(match( obj@usedObj[['clusters']], which(table(as.factor(obj@usedObj[['clusters']])) < 4 ) )) == T
			
			use <- obj
			use@usedObj[['clusters']] <- obj@usedObj[['clusters']][usable]
			use@usedObj[['mds.proj']] <- obj@usedObj[['mds.proj']][usable,]
			cols <- rainbow(max(as.numeric(obj@usedObj[['clusters']])))
			H <- Hkda( use@usedObj[['mds.proj']], use@usedObj[['clusters']], bw='plugin')
			kda.fhat <- kda( use@usedObj[['mds.proj']], use@usedObj[['clusters']],Hs=H, compute.cont=TRUE)
			try(plot(kda.fhat, size=0.001, colors = cols[as.numeric(names(table(use@usedObj[['clusters']])))] ),silent=F)
			try( writeWebGL(
							dir = file.path(obj@outpath,'densityWebGL'), 
							width=470, height=470, prefix='K', 
							template= system.file(file.path("densityWebGL.html"), 
									package = "Rscexv") 
							),
					silent=F 
			)
		}
)
