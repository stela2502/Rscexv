#' @name z.score.PCR.mad
#' @aliases z.score.PCR.mad,Rscexv-method
#' @rdname z.score.PCR.mad-methods
#' @docType methods
#' @description This normlization method is re-implemented from the MAST package.
#' @param dataObj the Rscexv object
#' @title description of function z.score.PCR.mad
#' @export 
setGeneric('z.score.PCR.mad', ## Name
		function (dataObj) { 
			standardGeneric('z.score.PCR.mad')
		}
)

setMethod('z.score.PCR.mad', signature = c ('Rscexv'),
		definition = function (dataObj) {
			if ( ! dataObj@zscored ){

			arrays <- unique(dataObj@samples$ArrayID)
			
			rem.inds <- NULL
			
			tab.pre <-NULL
			
			for (j in arrays ){
				
				mads <- NULL
				meds <- NULL
				
				tdat <- as.matrix(dataObj@data[which(dataObj@samples$ArrayID==j),])
				
				tab.new <- NULL
				
				for(i in 1:ncol(tdat)){
					vec <- tdat[which(tdat[,i]!=0),i]
					mads <- c(mads,mad(vec))
					meds <- c(meds,median(vec))
				}
				
				for(i in 1:ncol(tdat)){
					new.v <- (tdat[,i]-meds[i])/(1.48*mads[i])
					if(all(is.na(range(new.v))==T)){
						rem.inds <- c(rem.inds,i)
						plug.ind <- 1:length(tdat[,i])
					}
					else {
						plug.ind <- which(tdat[,i] == 0)
					}
					
					new.v[plug.ind] <- -20
					tab.new <- cbind(tab.new,new.v)
					
				}
				tab.pre <- rbind(tab.pre,tab.new)
				
			}
			
			rem.ind.fin <- as.numeric( names(table(rem.inds))[which(table(rem.inds)==length(arrays))]) ### might fail

						
			if ( length(rem.ind.fin) > 0 ){
				dataObj <- remove.genes (dataObj, rem.ind.fin)
				tab.pre <- tab.pre[,-rem.ind.fin]
			}
						
			colnames(tab.pre) <- colnames(dataObj@data) 
			dataObj@snorm <- dataObj@data 
			dataObj@data <- data.frame(tab.pre)
			dataObj@zscored <- T
			}
			else {
				print ( "Uncanged as data was already zscored!")
			}
			dataObj
			
		} 
)
