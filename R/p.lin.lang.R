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
		function ( x, groups.n, clus, n=1000 ) { 
			standardGeneric('p.lin.lang')
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

