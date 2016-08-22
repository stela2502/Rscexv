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

