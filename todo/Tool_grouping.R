#' @name group_1D
#' @aliases group_1D,Rscexv-method
#' @rdname group_1D-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param gene  TEXT MISSING
#' @param ranges  TEXT MISSING
#' @title description of function group_1D
#' @export 
setGeneric('group_1D', ## Name
	function (dataObj, gene, ranges) { ## Argumente der generischen Funktion
		standardGeneric('group_1D') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('group_1D', signature = c ('Tool_grouping'),
	definition = function (dataObj, gene, ranges) {
	userGroups <- group_1D_worker ( dataObj$PCR, gene, ranges)
	if ( max(userGroups$groupID) == 0 ){
		userGroups <- group_1D_worker ( dataObj$FACS, gene, ranges)
	}
	userGroups <- checkGrouping ( userGroups, dataObj )
	userGroups
} )
#' @name group_1D_worker
#' @aliases group_1D_worker,Rscexv-method
#' @rdname group_1D_worker-methods
#' @docType methods
#' @description 
#' @param ma  TEXT MISSING
#' @param gene  TEXT MISSING
#' @param ranges  TEXT MISSING
#' @title description of function group_1D_worker
#' @export 
setGeneric('group_1D_worker', ## Name
	function (ma, gene, ranges ) { ## Argumente der generischen Funktion
		standardGeneric('group_1D_worker') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('group_1D_worker', signature = c ('Tool_grouping'),
	definition = function (ma, gene, ranges ) {
	
	position <- which ( colnames(ma) == gene )
	
	userGroups <- data.frame( cellName = rownames(ma), userInput = rep.int(0, nrow(ma)), groupID = rep.int(0, nrow(ma)) )
	
	if ( length(position) > 0 ){
		min <- min(ma[,position])
		max <- max(ma[,position])+1
		
		ranges = ranges[order(ranges)]
		minor = 0
		now <- as.vector( which( ma[,position] >= min & ma[,position] < ranges[1] ))
		userGroups$userInput[now] = paste ('min <= x <',ranges[1] )
		userGroups$groupID[now] = 1
		for ( i in 2:length(ranges) ) {
			now <- as.vector( which( ma[,position] >= ranges[i-1] & ma[,position] < ranges[i] ))
			userGroups$userInput[now] = paste(ranges[i-1],'<= x <',ranges[i])
			if ( length(now) > 0 ){
				userGroups$groupID[now] = i
			}
			else {
				minor = minor + 1
			}
		}
		now <- as.vector( which( ma[,position] >= ranges[length(ranges)] & ma[,position] < max ))
		userGroups$userInput[now] = paste(ranges[length(ranges)],'<= x < max')
		userGroups$groupID[now] = length(ranges) +1
		userGroups <- checkGrouping ( userGroups )
	}
	
	userGroups
} )

#' @name regroup
#' @aliases regroup,Rscexv-method
#' @rdname regroup-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param group2sample  TEXT MISSING default= list ( '1' = c( 'Sample1'
#' @param 'Sample2'))  TEXT MISSING
#' @title description of function regroup
#' @export 
setGeneric('regroup', ## Name
	function ( dataObj, group2sample = list ( '1' = c( 'Sample1', 'Sample2' ) ) ) { ## Argumente der generischen Funktion
		standardGeneric('regroup') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('regroup', signature = c ('Tool_grouping'),
	definition = function ( dataObj, group2sample = list ( '1' = c( 'Sample1', 'Sample2' ) ) ) {
	userGroups <-  data.frame( cellName = rownames(dataObj$PCR), userInput = rep.int(0, nrow(dataObj$PCR)), groupID = rep.int(0, nrow(dataObj$PCR)) )
	n <- names(group2sample)
	n <- n[order( n )]
	minor = 0
	for ( i in 1:length(n) ){
		if ( sum(is.na(match(group2sample[[i]], userGroups$cellName))==F) == 0 ){
			minor = minor +1
		}
		else if (  sum(is.na(match(group2sample[[i]], userGroups$cellName))==F) < length( group2sample[[i]] ) ) {
			rows <- match(group2sample[[i]], userGroups$cellName)
			rows <- rows[-which(is.na(rows))]
			userGroups[ rows ,3] <- i
		}
		else { ## the whole group has failed??
			userGroups[ match(group2sample[[i]], userGroups$cellName) ,3] = i - minor
		}
	}
	if ( length(which(userGroups[,3] == 0)) > 0 ){
		userGroups[,3]
		system (paste('echo "', length(which(userGroups[,3] == 0)),"cells were not grouped using the updated grouping' > Grouping_R_Error.txt", collaps=" ") )
	}
	checkGrouping ( userGroups, dataObj )
} )
#' @name group_on_strings
#' @aliases group_on_strings,Rscexv-method
#' @rdname group_on_strings-methods
#' @docType methods
#' @description 
#' @param dataObj  TEXT MISSING
#' @param strings  TEXT MISSING default= c()
#' @title description of function group_on_strings
#' @export 
setGeneric('group_on_strings', ## Name
	function (dataObj, strings = c() ) { ## Argumente der generischen Funktion
		standardGeneric('group_on_strings') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('group_on_strings', signature = c ('Tool_grouping'),
	definition = function (dataObj, strings = c() ) {
	userGroups <-  data.frame( cellName = rownames(dataObj$PCR), userInput = rep.int(0, nrow(dataObj$PCR)), groupID = rep.int(0, nrow(dataObj$PCR)) )
	minor = 0	
	for ( i in 1:length(strings) ) {
		g <- grep(strings[i], userGroups$cellName)
		if ( length(g) == 0 ){
			system (paste ('echo "The group name',strings[i] ,'did not match to any sample" > Grouping_R_Error.txt', collaps=" ") )
			minor = minor +1
		}
		else {
			userGroups[g ,3] = i - minor
			userGroups[g ,2] = strings[i]
		}
	}
	checkGrouping ( userGroups, dataObj )
} )
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
	function (dataObj, fname='RandomForest_groupings.txt' ) { ## Argumente der generischen Funktion
		standardGeneric('createGroups_randomForest') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
} )
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
	function (dataObj, expected_grouping=10 ) { ## Argumente der generischen Funktion
		standardGeneric('createGeneGroups_randomForest') ## der Aufruf von standardGeneric sorgt für das Dispatching
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
} )
