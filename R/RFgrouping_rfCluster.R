#' @name rfCluster
#' @aliases 'rfCluster,Rscexv-method
#' @title rfCluster
#' @name rfCluster-methods
#' @docType methods
#' @description This fucntion uses the RFclust.SGE to create fandomForest based unsupervised clusters on a subset of the data.
#' @description Default is on 200 cells using all (provided) genes with 500 forests and 500 trees per forest for 5 repetitions.
#' @description You are asked to give a k numer of expected clusters (better too many than too little), classifies the total 
#' @description data using the 5 different unsupervised runs and all cluster ids from these runs are merged into the final cluster id.
#' @description This <summaryCol> will be part of the return objects samples table, together with a <usefulCol> where
#' @description all clusters with less than 10 cells have been merged into the 'gr. 0'.
#' @description The final results will be reported as new columns in the samples table containing the 'name'
#' @param x the single cells ngs object
#' @param email your email to use together with the SGE option
#' @param SGE whether to use the sun grid engine to calculate the rf grouping
#' @param rep how many repetitions for the random forest grouping should be run (default = 5)
#' @param slice how many processes should be started for each random forest clustering (default = 30)
#' @param bestColname the column name to store the results in
#' @param k the numer of expected clusters (better more than to view)
#' @param subset how many cells should be randomly selected for the unsupervised clustering (default = 200)
#' @param name if you want to run multiple RFclusterings on e.g. using different input genes you need to specify a name (default ='RFclust')
#' @param pics create a heatmap for each grouping that has been accessed (in the outpath folder; default = FALSE)
#' @param nforest the numer of forests to grow for each rep (defualt = 500)
#' @param ntree the numer of trees per forest (default = 500)
#' @param recover in case some problem has occured you could recover the objects from the saved files (default =F)
#' @return a Rscexv object including the results and storing the RF object in the usedObj list (bestColname)
#' @export 
setGeneric('rfCluster',
		function ( x, rep=5, SGE=F, email, k=16, slice=30, subset=200, pics=F ,nforest=500, ntree=500, name='RFclust', recover=F){
			standardGeneric('rfCluster')
		}
)
setMethod('rfCluster', signature = c ('Rscexv'),
		definition = function ( x, rep=5, SGE=F, email, k=16, slice=30, subset=200, pics=F ,nforest=500, ntree=1000, name='RFclust', recover=F) {
			summaryCol=paste( 'All_groups', name,sep='_')
			usefulCol=paste ('Usefull_groups',name, sep='_')
			n= paste(x@name, name,sep='_')
			m <- max(k)
			OPATH <- paste( x@outpath,"/",str_replace( x@name, '\\s', '_'), sep='')
			opath = paste( OPATH,"/RFclust.mp",sep='' )
			
			if ( ! dir.exists(OPATH)){
				dir.create( OPATH )
			}
			processed = FALSE
			single_res_col <- paste('RFgrouping',name)
			if ( is.null(x@usedObj[['rfExpressionSets']])){
				x@usedObj[['rfExpressionSets']] <- list()
				x@usedObj[['rfObj']] <- list()
			}
			for ( i in 1:rep) {
				tname = paste(n,i,sep='_')
				
				if ( is.null(x@usedObj[['rfExpressionSets']][[tname]]) ){
					## start the calculations!
					if ( dir.exists(opath)){
						if ( opath == '' ) {
							stop( "Are you mad? Not giving me an tmp path to delete?")
						}
						if ( ! recover ){
							system( paste('rm -f ',opath,"/*",tname,'*', sep='') )
						}
					}else {
						dir.create( opath )
					}
					if ( recover ){
						l <- function ( fn ) {
							load( fn )
							x
						}
						if ( file.exists(file.path( opath, paste(tname,'.RData',sep=''))) ){
							loaded <- l( file.path( opath, paste(tname,'.RData',sep='')) )
							x@usedObj[['rfObj']][[tname]] <- loaded
							x@usedObj[['rfExpressionSets']][[ tname ]] <-  remove.samples( x, 
									c(1:nrow(x@data))[is.na(match( make.names(rownames(x@data)), colnames(loaded@dat) )) ==T ]
							)
							## now add all analysis files into the object!
							v <- NULL
							x@usedObj[['rfObj']][[ tname ]]@RFfiles[[tname]] = NULL
							for ( i in 1:slice ){
								v =c(v,  file.path(opath,paste( paste('runRFclust',tname,i,sep='_'),'RData',sep='.') ) )
							}
							x@usedObj[['rfObj']][[ tname ]]@RFfiles[[tname]] = v
						}
						print ( "Please re-run to load the results into the object." )
					}
					else {
					total <- nrow(x@data)
					if ( total-subset <= 20 ) {
						stop( paste( 'You have only', total, 'samples in this dataset and request to draw random',subset, "samples, which leaves less than 20 cells to draw on random!") )
					}
					
					if ( length( x@usedObj[['rfExpressionSets']] ) < i  ) {
						x@usedObj[['rfExpressionSets']][[ tname ]] <- remove.samples( x, sample(c(1:total),total-subset) )
						x@usedObj[['rfObj']][[ tname ]] <- RFclust.SGE ( dat=data.frame(t(x@usedObj[['rfExpressionSets']][[ tname ]]@data)), SGE=SGE, slice=slice, email=email, tmp.path=opath, name= tname )
					}
					#names(x@usedObj[['rfExpressionSets']]) [i] <- tname
					#names(x@usedObj[['rfObj']]) [i] <- tname
					x@usedObj[['rfObj']][[ tname ]] <- runRFclust ( x@usedObj[['rfObj']][[ tname ]] , nforest=nforest, ntree=ntree, name=tname )
					if ( SGE){
						print ( "You should wait some time now to let the calculation finish! check: system('qstat -f') -> re-run the function")
					}
					else {
						print ( "You should wait some time now to let the calculation finish! -> re-run the function")
						print ( "check: system( 'ps -Af | grep Rcmd | grep -v grep')")
					}
					}
				}
				else {
					
					## read in the results
					try ( x@usedObj[['rfObj']][[ tname ]] <- runRFclust (
									x@usedObj[['rfObj']][[tname]] , 
									nforest=nforest, 
									ntree=ntree, 
									name=tname 
					) )
					if ( ! is.null(x@usedObj[['rfObj']][[ tname ]]@RFfiles[[tname]]) ){
						stop( "please re-run this function later - the clustring process has not finished!")
					}
					for ( a in k ){
						x@usedObj[["rfExpressionSets"]][[tname]]@samples <- 
								x@usedObj[["rfExpressionSets"]][[tname]]@samples[ ,
										is.na(match ( colnames(x@usedObj[["rfExpressionSets"]][[tname]]@samples), paste('group n=',a) ))==T 
								]
					}
					groups <- createGroups( x@usedObj[['rfObj']][[tname]], k=k, name=tname )
					x@usedObj[['rfExpressionSets']][[tname]]@samples <- cbind ( x@usedObj[['rfExpressionSets']][[tname]]@samples, groups[,3:(2+length(k))] )
										
					le <- ncol(x@usedObj[['rfExpressionSets']][[tname]]@samples)
					colnames(x@usedObj[['rfExpressionSets']][[tname]]@samples)[(le-length(k)+1):le] <- paste('group n=',k)
					
					## create the required RF object
					m <- max(k)
					x@usedObj[['rfExpressionSets']][[tname]] <- bestGrouping( x@usedObj[['rfExpressionSets']][[tname]], group=paste('group n=', m), bestColname = paste('OptimalGrouping',m ,name) )
					## the 'predictive RFobj group n=' object is created by the bestGrouping call
					x@samples[, paste( single_res_col, i) ] <-
							predict( x@usedObj[['rfExpressionSets']][[tname]]@usedObj[[paste( 'predictive RFobj group n=',m) ]], as.matrix(x@data) )
					x@usedObj[['colorRange']][[paste( single_res_col, i)]] <- rainbow( length(levels( x@samples[, paste( single_res_col, i) ])))
					if ( pics ){
						fn <- paste(OPATH,'/heatmap_rfExpressionSets_',i,'.png', sep='')
						png ( file=fn, width=800, height=1600 )
						gg.heatmap.list( x, groupCol=paste( single_res_col , i) )
						dev.off()
						print ( paste('heatmap stored in', fn) )
					}
					print ( paste("Done with cluster",i))
					processed = TRUE
				}
			}
			if ( processed ) {
				print( 'If you want to re-run with more trees/forests you need to use a new "name" option' ) 			
			}
			x		
		}
)
