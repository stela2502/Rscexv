#' @name create_p_values
#' @aliases create_p_values,Rscexv-method
#' @rdname create_p_values-methods
#' @docType methods
#' @description This is the main function called by the SCExV server.
#' @param obj the Rscexv object
#' @param boot for the boot strap approach - how many runs default= 1000
#' @param lin_lang_file the own stat outfile default='lin_lang_stats.xls'
#' @param sca_ofile the SingleCellsAssay p values outfile default="Significant_genes.csv"
#' @title description of function create_p_values
#' @export 
setGeneric('create_p_values', ## Name
		function ( obj, boot = 1000, lin_lang_file='lin_lang_stats.xls', sca_ofile="Significant_genes.csv" ) { ## Argumente der generischen Funktion
			standardGeneric('create_p_values') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('create_p_values', signature = c ('Rscexv'),
		definition = function ( obj, boot = 1000, lin_lang_file='lin_lang_stats.xls', sca_ofile="Significant_genes.csv" ) {
			stat_res = try ( SingleCellAssay_Pvalues ( obj, sca_ofile ))
			if ( obj@wFACS ){
				ma <- cbind( obj@data,  obj@facs) 
			}else {
				ma <- obj@data
			}
			
			groups.n = max(as.vector(obj@usedObj[['clusters']]))
			ma <- as.matrix(t(ma))
			n <- rownames(ma)
			cols = rainbow( groups.n )
			ma[which( ma == -20)] <- NA
			obj@usedObj[['stats']] <- vector('list', length=nrow( ma ))
			names(obj@usedObj[['stats']]) = rownames(ma)
			for ( i in 1:nrow( ma ) ) {
				obj@usedObj[['stats']][[i]] = p.lin.lang ( ma[i,], groups.n, obj@usedObj[['clusters']], n=boot )
			}
			obj@usedObj[['lin_lang']] <- write.stats( obj@usedObj[['stats']], file=file.path(obj@outpath,lin_lang_file) )
			write.table( cbind( stat_res, 'linear_model' = unlist( lapply( obj@usedObj[['stats']], function(x) { x$p_value } ))) , file=file.path(obj@outpath,'Summary_Stat_Outfile.xls') ,  sep='\t',quote=F )
			obj
		} 
)

#' @name p.lin.lang
#' @aliases p.lin.lang,Rscexv-method
#' @rdname p.lin.lang-methods
#' @docType methods
#' @description Calculates the linear hypothesis p value.
#' @param x data matrix
#' @param groups.n the amount of groups
#' @param clus the grouping vector
#' @param n the amount of boot strap operations default=1000
#' @title description of function p.lin.lang
#' @export 
setGeneric('p.lin.lang', ## Name
		function ( x, groups.n, clus, n=1000 ) { ## Argumente der generischen Funktion
			standardGeneric('p.lin.lang') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('p.lin.lang', signature = c ('numeric'),
		definition = function ( x, groups.n, clus, n=1000 ) {
			x[which( x == -20)] <- NA
			real_list <- vector ( 'list', groups.n)
			random_length <- vector ( 'numeric', groups.n)
			for( a in 1:groups.n){
				real_list[[a]]=x[which(clus == a)]
				random_length[a] = length(real_list[[a]])
				#real_list[[a]] = real_list[[a]][is.na(real_list[[a]]) == F]
			}
			## get one obj with all data
			stat_real = calc.lin.lang.4_list ( real_list )
			p_value = 1
			#print (stat_real$cor)
			if ( is.null(stat_real$cor) ) {
				stat_real$cor = 0.001
			}
			if ( is.na(stat_real$cor) ) {
				stat_real$cor = 0.001
			}
			if ( length(stat_real$cor) == 0){
				stop ( "Some crap happened!" )
			}
			
			if ( length(stat_real$weight) > 2 &&  stat_real$cor > 0.001 ){
				## bootstrap
				cmp = vector( 'numeric', n )
				medians <- vector ( 'numeric', groups.n)
				or <- vector ( 'numeric', groups.n)
				expo = 2
				for ( i in 1:n ) {
					for( gid in stat_real$groupIDs){
						tmp =sample(x,random_length[gid])
						tmp = tmp[is.na(tmp) ==F]
						or[gid] = length(tmp) / random_length[gid]
						if ( length(tmp) == 0 ){
							medians[gid] = 0
						}
						else if ( length(tmp) == 1 ){
							medians[gid] = tmp[1]
						}else {
							medians[gid] = median(tmp)
						}
					}
					cmp[i] = corr( cbind(or[stat_real$groupIDs], medians[stat_real$groupIDs]), w = stat_real$weight / sum(stat_real$weight) )
					if ( i %% 10^expo == 0 ){
						#	expo = expo + 1
						t <- boot_p_value ( cmp, stat_real$cor, i )
						#print (paste( "I check whether the p value (",t,") is already saturated",i))
						if ( t > 20/i ) {
							
							break
						}
					}
					if (is.na(cmp[i]) ){
						cmp[i]=1/n
					}
				} 
				#hist( cmp )
				#abline(v=stat_real$cor, col='red')
				stat_real$bootstrapedcor = cmp
				p_value <- boot_p_value ( cmp, stat_real$cor, sum( cmp != 0 ) )
				#print ( paste('Final p value = ',p_value, 'Using stat real cor =',stat_real$cor,"and n =",sum( cmp != 0 ),"values"))
			}
			else {
				p_value = 1 
			}
			stat_real$p_value = p_value
			stat_real
		} 
)

#' @name calc.lin.lang.4_list
#' @aliases calc.lin.lang.4_list,Rscexv-method
#' @rdname calc.lin.lang.4_list-methods
#' @docType methods
#' @description This function calculates the correlation between the median expression value of a group and the fraction of expression cells in the group.
#' @param l the list of values
#' @title description of function calc.lin.lang.4_list
#' @export 
setGeneric('calc.lin.lang.4_list', ## Name
		function ( l =list() ) { ## Argumente der generischen Funktion
			standardGeneric('calc.lin.lang.4_list') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('calc.lin.lang.4_list', signature = c ('list'),
		definition = function ( l =list() ) {
			n = length(l)
			or <- NULL
			medians <- NULL
			groupIDs<- NULL
			weight <- NULL
			for (i in 1:n) {
				data <- l[[i]][ is.na(l[[i]]) ==F ]
				if ( length(data) == 0) {
					data <- c(0)
				}
				data.min <- min(data)
				data.max <- max(data)
				if ( data.min != data.max ) {
					or <- c( or, length(data) / length(l[[i]] ))
					weight <- c( weight,  length(l[[i]] ) )
					medians <- c( medians,  median(data) )
					groupIDs <- c( groupIDs, i )
				}
			}
			ret <- list (  'weight' = weight, medians = medians, groupIDs = groupIDs, or = or , cor = NULL)
			if( length(weight) > 2 ){
				ret$cor <- corr( cbind(or, medians), w = weight / sum(weight) )
			}
			
			ret
		} 
)


#' @name boot_p_value
#' @aliases boot_p_value,Rscexv-method
#' @rdname boot_p_value-methods
#' @docType methods
#' @description calculates the p value from the boot strap distribution for the lin_lang analysis
#' @param cmp the boot strap list of correlation values
#' @param real_val the correlation value
#' @param i the numer of bootstrap runs
#' @title description of function boot_p_value
#' @export 
setGeneric('boot_p_value', ## Name
		function (cmp, real_val, i ) { ## Argumente der generischen Funktion
			standardGeneric('boot_p_value') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('boot_p_value', signature = c ('numeric'),
		definition = function (cmp, real_val, i ) {
			a <- length(which(cmp > real_val))
			if ( a == 0 ) { a<-1}
			p_value = a / i
			p_value
		} 
)

#' @name write.stats
#' @aliases write.stats,Rscexv-method
#' @rdname write.stats-methods
#' @docType methods
#' @description write a statistics table from the lin lang list
#' @param stats the lin lang list default= NULL
#' @param file the outfile default='lin_lang_stats.xls'
#' @title description of function write.stats
setGeneric('write.stats', ## Name
		function ( stats = NULL, file='lin_lang_stats.xls' ) { ## Argumente der generischen Funktion
			standardGeneric('write.stats') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('write.stats', signature = c ('list'),
		definition = function ( stats = NULL, file='lin_lang_stats.xls' ) {
			
			groupL <- function(x) {
				if ( ! is.vector(x$medians)){ x$medians = c(-1,-2) }
				if ( ! is.vector(x$groupIDs)){ x$groupIDs = c(-1,-2) }
				if ( ! is.vector(x$weight)){ x$weight = c(-1,-2) }
				c( x$cor, x$p_value, 
						paste(x$groupIDs[order(x$medians)], collapse =', '), 
						paste(x$medians[order(x$medians)], collapse =', '), 
						paste(x$weight[order(x$medians)], collapse =', ') 
				) }
			ma <- NULL
			if ( ! is.null(stats) ) {
				ma <- t(as.data.frame(lapply(stats, groupL )))
				rownames(ma) <- names(stats)
				colnames(ma)<- c('Correlation', 'p value', 'groups in order', 'median expression in group', 'weight of group' )
				write.table( ma, file=file ,  sep='\t',quote=F ) 
			}
			else {
				print ( "No starts to print!" )
			}
			ma
		} 
)


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
		function ( obj, ofile="Significant_genes.csv" ) { ## Argumente der generischen Funktion
			standardGeneric('SingleCellAssay_Pvalues') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('SingleCellAssay_Pvalues', signature = c ('Rscexv'),
		definition = function ( obj, ofile="Significant_genes.csv" ) {
			d <- as.matrix(obj@data)
			d[which(d==-20)] <- NA
			x <- as.matrix(d)
			d[is.na(d)] <- 0
			sca <- FromMatrix('SingleCellAssay', as.matrix(d), data.frame(wellKey=rownames(d)), data.frame(primerid=colnames(d)) )
			groups <- cData(sca)$GroupName <- obj@usedObj[['clusters']]
			zlm.output <- zlm.SingleCellAssay(~ GroupName, sca, method='glm', ebayes=T)
			zlm.lr <- lrTest(zlm.output,'GroupName')
			pvalue <- ggplot(melt(zlm.lr[,,'Pr(>Chisq)']), aes(x=primerid, y=-log10(value)))+ geom_bar(stat='identity')+facet_wrap(~test.type) + coord_flip()
			png ( file.path( obj@outpath,'Analysis1.png'), width=800, height=800)
			print(pvalue)
			dev.off()
			write.table( zlm.lr[,,'Pr(>Chisq)'], file=file.path(obj@outpath, ofile), sep='\t')
			zlm.lr[,,'Pr(>Chisq)']
		} 
)

