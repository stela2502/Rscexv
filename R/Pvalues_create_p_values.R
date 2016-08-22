#' @name create_p_values
#' @aliases create_p_values,Rscexv-method
#' @rdname create_p_values-methods
#' @docType methods
#' @description This is the main function called by the SCExV server.
#' @param obj the Rscexv object
#' @param boot for the boot strap approach - how many runs default= 1000
#' @param lin_lang_file the own stat outfile default='lin_lang_stats.xls'
#' @param sca_ofile the SingleCellsAssay p values outfile default="Significant_genes.csv"
#' @title description of function create_p_values
#' @export 
setGeneric('create_p_values', ## Name
		function ( obj, boot = 1000, lin_lang_file='lin_lang_stats.xls', sca_ofile="Significant_genes.csv" ) { ## Argumente der generischen Funktion
			standardGeneric('create_p_values') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('create_p_values', signature = c ('Rscexv'),
		definition = function ( obj, boot = 1000, lin_lang_file='lin_lang_stats.xls', sca_ofile="Significant_genes.csv" ) {
			stat_res = try ( SingleCellAssay_Pvalues ( obj, sca_ofile ))
			if ( obj@wFACS ){
				ma <- cbind( obj@data,  obj@facs) 
			}else {
				ma <- obj@data
			}
			
			groups.n = max(as.vector(obj@usedObj[['clusters']]))
			ma <- as.matrix(t(ma))
			n <- rownames(ma)
			cols = rainbow( groups.n )
			ma[which( ma == -20)] <- NA
			obj@usedObj[['stats']] <- vector('list', length=nrow( ma ))
			names(obj@usedObj[['stats']]) = rownames(ma)
			for ( i in 1:nrow( ma ) ) {
				obj@usedObj[['stats']][[i]] = p.lin.lang ( ma[i,], groups.n, obj@usedObj[['clusters']], n=boot )
			}
			obj@usedObj[['lin_lang']] <- write.stats( obj@usedObj[['stats']], file=file.path(obj@outpath,lin_lang_file) )
			write.table( cbind( stat_res, 'linear_model' = unlist( lapply( obj@usedObj[['stats']], function(x) { x$p_value } ))) , file=file.path(obj@outpath,'Summary_Stat_Outfile.xls') ,  sep='\t',quote=F )
			obj
		} 
)

