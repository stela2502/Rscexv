data(analyzed)

context('re-use.groups with FACS')

onwhat='Expression'
clusterby='MDS'
mds.type='PCA'
move.neg <- TRUE
plot.neg <- TRUE
beanplots = TRUE
plotsvg = 0
zscoredVioplot = 1
cmethod='ward.D'
LLEK='2'
ctype= 'hierarchical clust'
groups.n <- 5
data@outpath = tempdir()


data2 <- analyse.data ( 
		data,
		groups.n=4,  
		onwhat='Expression', 
		clusterby='MDS', 
		mds.type='PCA', 
		cmethod='ward.D', 
		LLEK='2',  
		ctype= 'hierarchical clust',  
		zscoredVioplot = zscoredVioplot, 
		move.neg = move.neg, 
		plot.neg=plot.neg, 
		beanplots=beanplots,
		useGrouping='ArrayID',
		plotsvg = 0
)

expect_that( colnames(data2@samples) , equals(colnames(data@samples)) )

data2 <- group_1D (data, 'CD41.Alexa.Fluor.488.A', c(0.5, 1, 1.5, 2, 2.5, 3 ) )


expect_that( 
		colnames(data2@samples) , 
		equals(c(colnames(data@samples) , 'CD41.Alexa.Fluor.488.A 1D Group' ) ) 
)

## check that the values are also OK
expect_that ( all( data2@facs[ which(data2@samples[,'CD41.Alexa.Fluor.488.A 1D Group'] == 1), 'CD41.Alexa.Fluor.488.A' ] <= 0.5), is_true() )
expect_that ( all( data2@facs[ which(data2@samples[,'CD41.Alexa.Fluor.488.A 1D Group'] == 2), 'CD41.Alexa.Fluor.488.A' ] <= 1), is_true() )
expect_that ( all( data2@facs[ which(data2@samples[,'CD41.Alexa.Fluor.488.A 1D Group'] == 3), 'CD41.Alexa.Fluor.488.A' ] <= 1.5), is_true() )
expect_that ( all( data2@facs[ which(data2@samples[,'CD41.Alexa.Fluor.488.A 1D Group'] == 4), 'CD41.Alexa.Fluor.488.A' ] <= 2), is_true() )
expect_that ( all( data2@facs[ which(data2@samples[,'CD41.Alexa.Fluor.488.A 1D Group'] == 5), 'CD41.Alexa.Fluor.488.A' ] <= 2.5), is_true() )
expect_that ( all( data2@facs[ which(data2@samples[,'CD41.Alexa.Fluor.488.A 1D Group'] == 6), 'CD41.Alexa.Fluor.488.A' ] <= 3), is_true() )
expect_that ( all( data2@facs[ which(data2@samples[,'CD41.Alexa.Fluor.488.A 1D Group'] == 7), 'CD41.Alexa.Fluor.488.A' ] >= 3), is_true() )

data2 <- group_1D (data2, 'Actb', c( 26,30,32))

expect_that ( all( data2@data[ which(data2@samples[,'Actb 1D Group'] == 1), 'Actb' ] <= 26), is_true() )
expect_that ( all( data2@data[ which(data2@samples[,'Actb 1D Group'] == 2), 'Actb' ] <= 30), is_true() )
expect_that ( all( data2@data[ which(data2@samples[,'Actb 1D Group'] == 3), 'Actb' ] <= 32), is_true() )
expect_that ( all( data2@data[ which(data2@samples[,'Actb 1D Group'] == 4), 'Actb' ] >= 32), is_true() )

gr_old <- split( data2@samples[,1], data2@samples[,4])
l <- list()
expect_that(names(gr_old), equals( paste(1:7)) )
for ( i in  c( '2','5','4','1','3','6','7') ) { l[[length(l)+1]] = gr_old[[i]]}
names(l) <- 1:7

data2@samples$reorderedGroup <- data2@samples[,4]
data2 <- regroup(data2, group2sample=l, name =  'reorderedGroup' )
t <- table (data2@samples[,c(4,5)])
expect_that(t[c(2,7+5,7*2+4,7*3+1,7*4+3,7*5+6,7*6+7)],equals( c(12,62,161,9,33,4,1)) )

data2 <- group_on_strings( data2, strings=c( ' P1', ' P2', ' P3' ) )
t <- table (data2@samples[,c(2,6)])
expect_that(t[c(1,2,3,4,8,12)],equals( c(5,5,3,89,89,91)) )

colColors = list(
		ArrayID=rainbow(max(data2@samples$ArrayID)), 
		'auto_clusters.1' = rainbow(max(data2@samples$'auto_clusters.1')), 
		'CD41.Alexa.Fluor.488.A 1D Group' = gray.colors( max(data2@samples$'CD41.Alexa.Fluor.488.A 1D Group')) 
) 

complexHeatmap( data2, ofile= 'notUsed', colGroups= c('ArrayID', 'auto_clusters.1', 'CD41.Alexa.Fluor.488.A 1D Group' ),colColors= colColors)

