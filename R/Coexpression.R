#' @name plotcoma
#' @aliases plotcoma,Rscexv-method
#' @rdname plotcoma-methods
#' @docType methods
#' @description This function plots the correlation plot for the grouping. 
#' @description It will create an svg file if the global variable plotsvg ==1.
#' @param dataObj  TEXT MISSING
#' @param fname  TEXT MISSING default='CorrelationPlot'
#' @title description of function plotcoma
#' @export 
setGeneric('plotcoma', ## Name
		function ( dataObj, fname='/CorrelationPlot' ) { ## Argumente der generischen Funktion
			standardGeneric('plotcoma') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plotcoma', signature = c ('Rscexv'),
		definition = function ( dataObj, fname='/CorrelationPlot' ) {
			if ( is.null(dataObj@usedObj[['coma']] ) ) {
				if ( max(dataObj@usedObj[['clusters']]) > 1){
					dataObj <- coexpressionMatrix ( dataObj )
				}else {
					print ("I can not calculate a co-expression dataset on 1 groups!" )
					return (0)
				}
			}
			png ( file=paste(dataObj@outpath,fname,'.png',sep=''), width=800, height=800 )
			corrplot(dataObj@usedObj[['coma']], order = "hclust", method = "square", hclust.method ='single' )
			dev.off()
			if ( plotsvg == 1 ) {
				devSVG ( file=paste(fname,'.svg',sep=''), width=6, height=6 )
				corrplot(dataObj@usedObj[['coma']], order = "hclust", method = "square", hclust.method ='single' )
				dev.off()
			}
		}  
)
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

#' @name coexpressGenes
#' @aliases coexpressGenes,Rscexv-method
#' @rdname coexpressGenes-methods
#' @docType methods
#' @description  calculates the coexpression for all genes in all groups in the data set
#' @param dataObj the Rscexv object
#' @param grouping the column in the samples table that describes the grouping to use
#' @title description of function coexpressGenes
#' @export 
setGeneric('coexpressGenes', ## Name
		function ( dataObj, grouping=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('coexpressGenes') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('coexpressGenes', signature = c ('Rscexv'),
		definition = function ( dataObj, grouping=NULL ) {
			
			cor.funct <- function ( ma ){
				ma <- ma[, which( apply( ma, 2, function ( x) { length(which( x != 0)) }) > 9 )]
				if ( ncol(ma) < 2 ) {
					NULL
				}
				else {
					cor.t <- cor( ma , method='spearman')
					cor.p <- cor.t
					diag(cor.p) <- 1
					for ( i in 1:(ncol(ma)-1) ) {
						for (a in (i+1):ncol(ma) ) {
							if ( length( as.vector(ma[,i]) ) != length(as.vector(ma[,a]))){
								browser()
							}
							cor.p[a,i] <- cor.test( as.vector(ma[,i]), as.vector(ma[,a]),method='spearman')$p.value
						}
					}
					cor.t.m <- melt(cor.t)
					cor.p.m <- melt(cor.p)
					cor.t.m <- cbind(cor.t.m, cor.p.m[,3])
					cor.t.m <- cor.t.m[which(cor.t.m[,4] < 0.05), ]
					cor.t.m
				}
			}
			ret <- NULL
			for (i in unique(as.vector(dataObj@samples[,grouping] ))){
				t <- cor.funct ( dataObj@snorm[which(dataObj@samples[,grouping]== i),] )
				if ( ! is.null(t)){
					if ( nrow(t) > 0  ){
						t[,5] <- i
						ret <- rbind(ret, t)
					}
				}
			}
			i <- 0
			t <- cor.funct ( dataObj@snorm )
			if ( ! is.null(t)){
				if ( nrow(t) > 0  ){
					t[,5] <- i
					ret <- rbind(ret, t)
				}
			}
			colnames(ret) <- c('Source.Node','Target.Node', 'rho', 'p.value','Group' )
			ret
		} 
)

