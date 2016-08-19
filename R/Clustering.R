#' @name mds.and.clus
#' @aliases mds.and.clus,Rscexv-method
#' @rdname mds.and.clus-methods
#' @docType methods
#' @description this function calculates the MDS and creates the grouping if necessary.
#' @param dataObj  TEXT MISSING
#' @param clusterby cluster by raw or MDS default="raw"
#' @param mds.type which MDS function to use default="PCA"
#' @param groups.n the numer of groups to create
#' @param LLEK the LLEK option for the LLE and ISOMAP mds functions default=2
#' @param cmethod the hyclust method to use default='ward.D'
#' @param ctype the clustering method ('hierarchical clust' or kmeans) default='hierarchical clust'
#' @param onwhat cluster on FACS or Expression data default="Expression"
#' @param useGrouping do not calculate a new grouping - use this column in the samples table (default=NULL)
#' @title description of function mds.and.clus
#' @export 
setGeneric('mds.and.clus', ## Name
		function (dataObj, ..., clusterby="raw", mds.type="PCA", groups.n, LLEK=2, cmethod='ward.D', 
				ctype='hierarchical clust',onwhat="Expression",useGrouping=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('mds.and.clus') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('mds.and.clus', signature = c ('Rscexv'),
		definition = function (dataObj, ..., clusterby="raw", mds.type="PCA", groups.n, LLEK=2, 
				cmethod='ward.D', ctype='hierarchical clust',onwhat="Expression", useGrouping=NULL ) {
			if(onwhat=="Expression"){
				tab <- as.matrix(dataObj@data)
			} 
			else {
				print ( paste ( "I work on the FACS data!" ) )
				tab <- as.matrix(dataObj@facs)
			}
			mds.proj <- NULL
			pr <- NULL
			system ( 'rm  loadings.png' )
			if(mds.type == "PCA"){
				pr <- prcomp(tab)
				mds.proj <- pr$x[,1:3]
				png ( file=file.path( dataObj@outpath,'loadings.png'), width=1000, height=1000 )
				plot (  pr$rotation[,1:2] , col='white' );
				text( pr$rotation[,1:2], labels= rownames(pr$rotation), cex=1.5 )
				dev.off()
				write.table( cbind( Genes = rownames(pr$rotation), pr$rotation[,1:2] ), 
						file=file.path( dataObj@outpath,'gene_loadings.xls') , row.names=F, sep='\t',quote=F )
				#	mds.trans <- prcomp(t(tab))$x[,1:3]
			} else if ( mds.type == "LLE"){
				mds.proj <- LLE( tab, dim = 3, k = as.numeric(LLEK) )
				#	mds.trans <- LLE( t(tab), dim = 3, k = as.numeric(LLEK) )
			}else if ( mds.type == "ISOMAP"){
				mds.proj <- Isomap( tab, dim = 3, k = as.numeric(LLEK) )$dim3
				#	mds.trans <- Isomap( t(tab), dim = 3, k = as.numeric(LLEK) )$dim3
			}
			else {
				print( paste("Sory I can not work on the option",mds.type) )
			}
			dataObj@usedObj[['mds.proj']] <- mds.proj
			
			dataObj <- clusters ( dataObj, onwhat=onwhat, clusterby=clusterby, groups.n = groups.n,
					ctype = ctype, cmethod=cmethod, useGrouping=useGrouping )
			
			dataObj
		} 
)
#' @name clusters
#' @aliases clusters,Rscexv-method
#' @rdname clusters-methods
#' @docType methods
#' @description create the grouping based on the MDS or raw data.
#' @param dataObj  TEXT MISSING
#' @param clusterby  TEXT MISSING default="raw"
#' @param groups.n  TEXT MISSING default= 3
#' @param ctype  TEXT MISSING default='hierarchical clust'
#' @param onwhat  TEXT MISSING default="Expression"
#' @param cmethod  TEXT MISSING default='ward.D'
#' @param useGrouping do not calculate a new grouping - use this column in the samples table (default=NULL)
#' @title description of function clusters
#' @export 
setGeneric('clusters', ## Name
		function (dataObj,clusterby="raw", useGrouping=NULL, mds.proj=NULL,groups.n = 3,
				ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D') { ## Argumente der generischen Funktion
			standardGeneric('clusters') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('clusters', signature = c ('Rscexv'),
		definition = function (dataObj,clusterby="raw", useGrouping=NULL, groups.n = 3,
				ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D' ) {
			## custering	
			clusters <- NULL
			hc <- NULL
			if(onwhat=="Expression"){
				tab <- dataObj@data
			}
			else {
				print ( paste ( "I work on the FACS data!" ) )
				tab <- dataObj@facs
			}
			if ( ! is.null(useGrouping) ) {
				clusters <- dataObj@samples[,useGrouping]
			}else if(clusterby=="MDS"){
				if ( ctype=='hierarchical clust'){
					hc <- hclust(dist( dataObj@usedObj[['mds.proj']] ),method = cmethod)
					clusters <- cutree(hc,k=groups.n)
				}else if (  ctype=='kmeans' ) {
					hc <- hclust(dist( dataObj@usedObj[['mds.proj']] ),method = cmethod)
					clusters <- kmeans( dataObj@usedObj[['mds.proj']] ,centers=groups.n)$cluster
				}
				else { stop( paste('ctype',ctype, 'unknown!' ) )}
			}
			else{#...do mds on tab
				if ( ctype=='hierarchical clust'){
					hc <- hclust(as.dist( 1- cor(t(tab), method='pearson') ),method = cmethod)
					clusters <- cutree(hc,k=groups.n)
				}else if (  ctype=='kmeans' ) {
					hc <- hclust(as.dist( 1- cor(t(tab), method='pearson') ),method = cmethod)
					clusters <- kmeans( dataObj@usedObj[['mds.proj']] ,centers=groups.n)$cluster
				}
				else { stop( paste('ctype',ctype, 'unknown!' ) )}
			}
			if ( is.null(useGrouping) ){
				png ( file=file.path( dataObj@outpath,'hc_checkup_main_clustering_function.png'), width=1600, height=800 )
				plot ( hc);
				dev.off()
				if(is.null(dataObj@usedObj[['auto_clusters']])){
					dataObj@usedObj[['auto_clusters']] = 0
				}
				dataObj@usedObj[['auto_clusters']] <- dataObj@usedObj[['auto_clusters']] +1
				dataObj@samples <- cbind ( dataObj@samples, clusters )
				n <- paste( 'auto_clusters', 
						dataObj@usedObj[['auto_clusters']] ,sep='.')
				colnames(dataObj@samples)[ncol(dataObj@samples)] = n
				dataObj <- colors_4(dataObj, n )
			}
			dataObj@usedObj[['clusters']] <- clusters
			dataObj@usedObj[['hc']] <- hc
			dataObj
		}
)

#' @name quality_of_fit
#' @aliases quality_of_fit,Rscexv-method
#' @rdname quality_of_fit-methods
#' @docType methods
#' @description Calculates a quality of fit
#' @param obj  TEXT MISSING
#' @title description of function quality_of_fit
#' @export 
setGeneric('quality_of_fit', ## Name
		function ( obj ) { ## Argumente der generischen Funktion
			standardGeneric('quality_of_fit') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('quality_of_fit', signature = c ('Rscexv'),
		definition = function ( obj ) {
			test <- as.matrix(obj@data)
			test[which(test ==  -20 ) ] = NA
			ret <- list ( 'per_expression' = apply(test,2, difference, obj ) )
			ret$Expression = round(sum(ret$per_expression))
			if ( obj@wFACS ) {
				test <- obj@facs
				ret$per_FACS = apply(test,2, difference, obj ) 
				ret$FACS = round(sum(ret$per_FACS))
			}
			else {
				ret$per_FACS <- NA
				ret$FACS <- NA
			}
			ret
		} 
)


#' @name difference
#' @aliases difference,Rscexv-method
#' @rdname difference-methods
#' @docType methods
#' @description This function calculates the 
#' @param x  TEXT MISSING
#' @param obj  TEXT MISSING
#' @title description of function difference
#' @export 
setGeneric('difference', ## Name
		function ( x, obj ) { ## Argumente der generischen Funktion
			standardGeneric('difference') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('difference', signature = c ('numeric'),
		definition = function ( x, obj ) {
			ret = 0 
			for ( i in 1:groups.n  ) {
				a <- x[which( obj@usedObj[['clusters']] == i)]
				a <- a[- (is.na(a))==F]
				if ( length(a) > 1 ) {  ret = ret + sum( (a- mean(a) )^2 ) }
			}
			ret
		} 
)