#' @name norm.PCR
#' @aliases norm.PCR,Rscexv-method
#' @rdname norm.PCR-methods
#' @docType methods
#' @description This function normalizes the expression data
#' @param x the Rscexv object 
#' @param meth The normalization method; one of "none","mean control genes","max expression","median expression","quantile" default="none"
#' @param ctrl the control genes default=NA
#' @param max.cyc this function does also invert the Ct values this is the value that is used as maximum default=NA
#' @title description of function norm.PCR
#' @export 
setGeneric('norm.PCR', ## Name
		function (x,meth=c("none","mean control genes","max expression","median expression","quantile"),ctrl=NA,max.cyc) { ## Argumente der generischen Funktion
			standardGeneric('norm.PCR') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('norm.PCR', signature = c ('Rscexv'),
		definition = function (x, meth=c("none","mean control genes","max expression","median expression","quantile"),ctrl=NA,max.cyc) {
			if ( !x@norm ){
			tab.new <- NULL
			tab.na <- which(x@data == max.cyc )
			if(meth=="none"){
				tab.new <- as.matrix(x@data)
			}
			else if(meth=="mean control genes"){
				if(length(ctrl)>0){
					if ( length(ctrl)>1 ) {
						mean.ctrl <- apply(x@data[,ctrl],1,mean)
					}else {
						mean.ctrl <- as.vector(x@data[,ctrl])
					}
					for(i in 1:nrow(x@data)){
						tab.new <- rbind(tab.new,(x@data[i,]-mean.ctrl[i]))
					}
					x@samples$norm_to <- mean.ctrl
				}
			}
			else if (meth== "max expression" ){
				max.expr <- apply( x@data,1,min )
				for(i in 1:nrow(x@data)){
					tab.new <- rbind(tab.new,(x@data[i,]-max.expr[i]))
				} 
				x@samples$norm_to <- max.expr
			} 
			else if (meth== "median expression" ){
				my.median <- function (x, max.cyc)  {median( x[which( x != max.cyc )] ) } 
				median.expr <- apply( x@data,1,my.median, max.cyc  )
				for(i in 1:nrow(x@data)){
					tab.new <- rbind(tab.new,(x@data[i,]-median.expr[i]))
				}
				x@samples$norm_to <- median.expr
			} 
			else if(meth=="quantile"){
				rank.normalize <- function(ap) {
					o.ap <- apply(ap,2,order)
					ap.o <- NULL
					for (i in 1:ncol(o.ap))
						ap.o <- cbind(ap.o,ap[o.ap[,i],i])
					m.ap <- apply(ap.o,1,median)
					ap.n <- array(0,dim(ap))
					for (i in 1:ncol(ap.n))
						ap.n[o.ap[,i],i] <- m.ap
					ap.n
				}
				tab.new <- rank.normalize(x@data)
				colnames(tab.new) <- colnames(x@data)
			}
			else {
				stop(paste ("norm method",meth, "is not implemented!" ) )
			}
			if ( is.null(rownames(tab.new))){
				rownames(tab.new) <- rownames(x@data)
			}
			
			#tab.ret <- as.matrix(max.cyc-tab.new)
			tab.ret <- as.matrix(max.cyc-tab.new)
			tab.ret[tab.na] <- 0
			x@norm = T
			x@raw <- x@data
			x@data <- data.frame(tab.ret)
			x@normFunct <- meth
			}
			else {
				print ("Unchanged as data was already normalized")
			}
			x
		} 
)

