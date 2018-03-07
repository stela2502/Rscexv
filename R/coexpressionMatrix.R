#' @name coexpressionMatrix
#' @aliases coexpressionMatrix,Rscexv-method
#' @rdname coexpressionMatrix-methods
#' @docType methods
#' @description this function calculates the coexpression matrix for the expression data.
#' @param dataObj the Rscexv object
#' @title description of function coexpressionMatrix
#' @export 
setGeneric('coexpressionMatrix', ## Name
		function ( dataObj ) { ## Argumente der generischen Funktion
			standardGeneric('coexpressionMatrix') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('coexpressionMatrix', signature = c ('Rscexv'),
		definition = function ( dataObj ) {
			if ( is.null( dataObj@usedObj[['clusters']] ) ) {
				print ( 'Please calculate the clusters information first');
				return ;
			}
			# data$PCR cols == genes; rows == cells
			m <- max(dataObj@usedObj[['clusters']])
			coma <- matrix ( rep(0, m*m), ncol=m, nrow=m)
			mm <-  matrix ( rep(0,m * ncol(dataObj@data)), nrow=m)
			for ( i in 1:m ){
				## calculate the mean expression for each gene over all cells of the group
				if ( length(which(dataObj@usedObj[['clusters']] == i )) == 1 ) {
					mm[i,] <- dataObj@data[which(dataObj@usedObj[['clusters']] == i ),]
				}else{
					mm[i,] <- apply( dataObj@data[which(dataObj@usedObj[['clusters']] == i ),],2,mean)
				}
			}
			rownames(mm) <- paste('Group',1:m)
			coma <- cor(t(mm))
			dataObj@usedObj[['coma']] <- coma
			write.table (cbind(groups = rownames(coma), coma) , 
					file=file.path( dataObj@outpath,'correlation_matrix_groups.xls'), sep='\t',  row.names=F,quote=F )
			colnames(mm) <- colnames(dataObj@data)
			write.table (cbind(groups = rownames(mm), mm) , 
					file=file.path( dataObj@outpath,'mean_expression_per_groups.xls'), sep='\t',  row.names=F,quote=F )
			dataObj
		} 
)

