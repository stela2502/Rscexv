



#' @name createGroups_randomForest
#' @aliases createGroups_randomForest,Rscexv-method
#' @rdname createGroups_randomForest-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param fname  TEXT MISSING default='RandomForest_groupings.txt'
#' @title description of function createGroups_randomForest
#' @export 
setGeneric('createGroups_randomForest', ## Name
	function (dataObj, fname='RandomForest_groupings.txt' ) { 
		standardGeneric('createGroups_randomForest')
	}
)

setMethod('createGroups_randomForest', signature = c ('Tool_grouping'),
	definition = function (dataObj, fname='RandomForest_groupings.txt' ) {
	## load('RandomForestdistRFobject.RData') <- this has to be done before calling this function!!
	persistingCells <- rownames( dataObj$PCR )
	if ( exists('distRF') ) {
		expected_groupings <- unique(scan ( fname ))
		for ( i in 1:length(expected_groupings) ) {
			res = pamNew(distRF$cl1, expected_groupings[i] )
			N <- names( res )
			## probably some cells have been kicked in the meantime - I need to kick them too
			N <- intersect( persistingCells, N )
			userGroups <- matrix(ncol=3, nrow=0)
			for ( a in 1:length(N) ){
				userGroups <- rbind (userGroups, c( N[a], 'no info', as.numeric(res[[N[a]]]) ) )
			}
			colnames(userGroups) <- c('cellName', 'userInput',  'groupID' )
			## write this information into a file that can be used as group
			userGroups = data.frame( userGroups)
			save ( userGroups , file= paste("forest_group_n", expected_groupings[i],'.RData', sep=''))
			
			fileConn<-file(paste("Grouping.randomForest.n",expected_groupings[i],".txt", sep="") )
			writeLines(c(paste("load('forest_group_n",expected_groupings[i],".RData')",sep=""),
							"userGroups <- checkGrouping ( userGroups[is.na(match(userGroups$cellName, rownames(data.filtered$PCR) ))==F, ], data.filtered )" 
						), fileConn)
			close(fileConn)
		}
	}
} 
)

#' @name createGeneGroups_randomForest
#' @aliases createGeneGroups_randomForest,Rscexv-method
#' @rdname createGeneGroups_randomForest-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param expected_grouping  TEXT MISSING default=10
#' @title description of function createGeneGroups_randomForest
#' @export 
setGeneric('createGeneGroups_randomForest', ## Name
	function (dataObj, expected_grouping=10 ) { 
		standardGeneric('createGeneGroups_randomForest')
	}
)

setMethod('createGeneGroups_randomForest', signature = c ('Tool_grouping'),
	definition = function (dataObj, expected_grouping=10 ) {
	## load('RandomForestdistRFobject_genes.RData') <- this has to be done before calling this function!!
	persistingGenes <- colnames( dataObj$PCR )
	if ( round(length(persistingGenes)/4) < expected_grouping ){
		expected_grouping <- round(length(persistingGenes)/4)
	}
	if (expected_grouping < 2 ){
		expected_grouping <- 2
	}
	if ( exists('distRF') ) {
			res = pamNew(distRF$cl1, expected_grouping )
			N <- names( res )
			## probably some cells have been kicked in the meantime - I need to kick them too
			N <- intersect( persistingGenes , N )
			geneGroups <- matrix(ncol=3, nrow=0)
			for ( a in 1:length(N) ){
				geneGroups <- rbind (geneGroups, c( N[a], 'no info', as.numeric(res[[N[a]]]) ) )
			}
			colnames(geneGroups) <- c('geneName', 'userInput',  'groupID' )
			## write this information into a file that can be used as group
			geneGroups = data.frame( geneGroups)
			save ( geneGroups , file= paste("forest_gene_group_n", expected_grouping,'.RData', sep=''))
			
			fileConn<-file(paste("Gene_grouping.randomForest.txt", sep="") )
			writeLines(c(paste("load('forest_gene_group_n",expected_grouping,".RData')",sep=""),
							"geneGroups <- checkGrouping ( geneGroups[is.na(match(geneGroups$geneName, colnames(data.filtered$PCR) ))==F, ] )",
							"write.table( geneGroups[order(geneGroups[,3]),], file='GeneClusters.xls' , row.names=F, sep='\t',quote=F )"
					), fileConn)
			close(fileConn)
		
	}
}
)



