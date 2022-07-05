#' @name read.PCR.heatmap
#' @aliases read.PCR.heatmap,Rscexv-method
#' @rdname read.PCR.heatmap-methods
#' @docType methods
#' @description reads one fluidigm data file in heatmap format
#' @param fname  TEXT MISSING
#' @param use_pass_fail  TEXT MISSING default=T
#' @title description of function read.PCR.heatmap
#' @export 
setGeneric('read.PCR.heatmap', ## Name
		function (fname,use_pass_fail=T) { 
			standardGeneric('read.PCR.heatmap')
		}
)

setMethod('read.PCR.heatmap', signature = c ('character'),
		definition = function (fname,use_pass_fail=T) {
			
			top20 <- readLines(fname)
			
			ftab <- NULL
			
			line.start <- grep(",,1,2,3,4,5,6,7",top20)
			
			if ( length(line.start) == 0){
				system( paste('echo "File format error: file',fname, 'not readable\n" >> R_file_read_error.txt' ) )
			}
			else{
				tmp <- read.delim(fname,sep=',')
				line.end <- grep("Quality Results,,,,,,,",top20)-1[1]
				while (tmp[line.end,2] == '' ){
					line.end = line.end -1 
				}
				ftab <- tmp[(line.start[1]+1):(line.end),3:ncol(tmp)]
				rownames(ftab) <- tmp[(line.start[1]+1):(line.end),2]
				colnames(ftab) <- as.vector(t(tmp[line.start[1],3:ncol(tmp)]))
				for ( i in 1:ncol(ftab) ) {
					ftab[,i] <- as.numeric(as.vector(ftab[,i]))
				}
				if ( use_pass_fail ){
					## now we read the quality control data
					line.end <-nrow(tmp)
					while ( length(grep( 'Pass',t(tmp[line.end,3:ncol(tmp)]) ) ) == 0 ){ 
						line.end = line.end -1 
					}
					qval <- tmp[(line.start[2]+1):line.end,3:ncol(tmp)]
					rownames(qval) <- tmp[(line.start[2]+1):line.end,2]
					colnames(qval) <- as.vector(t(tmp[line.start[2],3:ncol(tmp)]))
					for ( i in 1:ncol(ftab) ) {
						ftab[which(qval[,i]=='Fail'),i] <- 999
					}
				}
				as.matrix(ftab)
			}
		} 
)

		
