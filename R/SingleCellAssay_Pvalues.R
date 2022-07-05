#' @name SingleCellAssay_Pvalues
#' @aliases SingleCellAssay_Pvalues,Rscexv-method
#' @rdname SingleCellAssay_Pvalues-methods
#' @docType methods
#' @description USes the SingleCellAssay package to calculate p values.
#' @param obj the Rscexv object
#' @param ofile the outfile default="Significant_genes.csv"
#' @title description of function SingleCellAssay_Pvalues
#' @export 
setGeneric('SingleCellAssay_Pvalues', ## Name
		function ( obj, ofile="Significant_genes.csv" ) { 
			standardGeneric('SingleCellAssay_Pvalues')
		}
)

setMethod('SingleCellAssay_Pvalues', signature = c ('Rscexv'),
		definition = function ( obj, ofile="Significant_genes.csv" ) {
			d <- as.matrix(obj@data)
			d[which(d==-20)] <- NA
			x <- as.matrix(d)
			d[is.na(d)] <- 0
			sca <- MAST::FromMatrix(class='SingleCellAssay', 
					#exprsArray= t(as.matrix(d)), 
					exprsArray= as.matrix(d), 
					cData= data.frame(
							wellKey=rownames(d),
							GroupName = obj@usedObj[['clusters']] 
					), 
					fData=data.frame(primerid=colnames(d)) 
			)
			#sca <- FromMatrix('SingleCellAssay', as.matrix(d), data.frame(wellKey=rownames(d)), data.frame(primerid=colnames(d)) )
			zlm.output <- MAST::zlm.SingleCellAssay(~ GroupName, sca, method='glm', ebayes=T)
			zlm.lr <- MAST::lrTest(zlm.output,'GroupName')
			pvalue <- ggplot(melt(zlm.lr[,,'Pr(>Chisq)']), aes(x=primerid, y=-log10(value)))+ geom_bar(stat='identity')+facet_wrap(~test.type) + coord_flip()
			png ( file.path( obj@outpath,'Analysis1.png'), width=800, height=800)
			print(pvalue)
			dev.off()
			write.table( zlm.lr[,,'Pr(>Chisq)'], file=file.path(obj@outpath, ofile), sep='\t')
			zlm.lr[,,'Pr(>Chisq)']
		} 
)

