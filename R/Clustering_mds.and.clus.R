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
			this.k <- paste(onwhat,mds.type)
			if ( is.null( dataObj@usedObj$MDSkey) ) {
				dataObj@usedObj$MDSkey = "none"
			}
			if ( (dataObj@usedObj$MDSkey != this.k) ||  all.equal( rownames(dataObj@usedObj$mds.proj), rownames(dataObj@data) )==F ) {
				dataObj@usedObj$MDSkey = mds.type
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
			}else if ( mds.type == "ZIFA" ) {
				print ( "Running external python script to apply ZIFA dimensional reduction (PCR data only)" )
				ZIFA <- 40.000001 - dataObj@raw
				write.table( ZIFA, file="ZIFA_input.dat", sep=" ", col.names=F, row.names=F , quote=F)
				write( c("from ZIFA import ZIFA","from ZIFA import block_ZIFA", "import numpy as np",
						"Y = np.loadtxt('ZIFA_input.dat')", "Z, model_params = ZIFA.fitModel( Y, 3 )", 
						"np.savetxt('TheMDS_ZIFA.xls', Z )" ), 
				       file= 'ZIFA_calc.py' )
			    system( "python ZIFA_calc.py" )
				Sys.sleep(5)
				mds.proj <- read.delim( "TheMDS_ZIFA.xls", sep=' ', header=F)
				rownames(mds.proj) <- rownames(ZIFA)
				colnames(mds.proj) <- c( 'x','y','z')
				
			}
			else {
				print( paste("Sory I can not work on the option",mds.type) )
			}
			dataObj@usedObj[['mds.proj']] <- mds.proj
			}
			dataObj <- clusters ( dataObj, onwhat=onwhat, clusterby=clusterby, groups.n = groups.n,
					ctype = ctype, cmethod=cmethod, useGrouping=useGrouping )
			
			dataObj
		} 
)
