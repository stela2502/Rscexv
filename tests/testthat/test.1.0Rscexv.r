PCR <- NULL
for ( f in c( 'PCR_Array1.csv','PCR_Array2.csv','PCR_Array3.csv') ) {
	PCR <- c(PCR, system.file(file.path(f), 
					package = "Rscexv") )
}

FACS = NULL
for ( f in c( 'Array1_Index_sort.csv', 'Index_sort_Array2.csv','Index_sort_Array3.csv') ) {
	FACS <- c(FACS, system.file(file.path(f), 
					package = "Rscexv") )
}

## test read and normalize data

negContrGenes <- NULL
max.value=40
ref.genes=c("Actb", "Gapdh")
max.ct=25
max.control=0
norm.function='none'
negContrGenes=NULL
use_pass_fail = T
methods <- c("mean control genes","max expression","median expression","quantile")

InbuiltData <- Rscexv( PCR, FACS, use_pass_fail)
save(InbuiltData, file=file.path("..","..","data","InbuiltData.RData") )

## otherwise I should be able to use the inbuilt data - right?

context('createDataObj with FACS')
InbuiltData@outpath= tempdir()

data <- kick.expressed.negContr.samples(InbuiltData, negContrGenes )

data <- plug.999(data, max.value ) ## does nothing for pre-processed data

data <- filter.on.controls.no.inv(data,ref.genes,max.ct,max.control)	

expect_that( c("NTC2.P0","NTC1.P0","NTC2.P1","NTC1.P1","NTC4.P2"," NTC3.P2"), equals(setdiff( InbuiltData@samples[,1], data@samples[,1]) ))

data.filtered <- sd.filter( data )

plot.histograms( data.filtered ) ## this is needed for the web tool

opath = file.path(data.filtered@outpath, 'preprocess')
for ( n in c( colnames(data.filtered@data), colnames(data.filtered@facs) ) ) {
	expect_that( file.exists( file.path( opath, paste(n,'png',sep='.'))), is_true())
}

for ( i in methods ) { 
	t <- norm.PCR(data.filtered,i,max.cyc=max.value, ctrl=ref.genes ) 
	## the original data is stored in the raw slot
	expect_that(t@raw, equals(data.filtered@data))
	expect_that( dim(t@data), equals(dim(data.filtered@data)))
}

data.filtered <- norm.PCR(data.filtered,norm.function,max.cyc=max.value, ctrl=ref.genes )
expect_that( dim(data.filtered@data), equals( c(282,95) ))


## test analysis no grouping

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

context('analyse.data with FACS')

data <- analyse.data ( 
		data.filtered,
		groups.n=groups.n,  
		onwhat='Expression', 
		clusterby='MDS', 
		mds.type='PCA', 
		cmethod='ward.D', 
		LLEK='2',  
		ctype= 'hierarchical clust',  
		zscoredVioplot = zscoredVioplot, 
		move.neg = move.neg, 
		plot.neg=plot.neg, 
		beanplots=beanplots
)

expect_that( 
		names(data@usedObj),
		equals( c("mds.proj","auto_clusters","clusters","hc","colors","quality_of_fit","for.plot"))
)

## check that all files have been created....

opath = data.filtered@outpath
for ( n in c( colnames(data@data), colnames(data@facs) ) ) {
	try(expect_that( file.exists( file.path( opath, paste(n,'png',sep='.')) ), is_true()))
}

xls.files <- c( "2D_data_color.xls","2D_data.xls","correlation_matrix_groups.xls" , "gene_loadings.xls","mean_expression_per_groups.xls", "merged_data_Table.xls" ) 

for ( n in xls.files ) {
	try(expect_that( file.exists( file.path( opath, n ) ), is_true()))
}

xls.files <- c("facs_color_groups_Heatmap.png","PCR_color_groups_Heatmap.png","facs_Heatmap.png","PCR_Heatmap.png")
for ( n in xls.files ) {
	try(expect_that( file.exists( file.path( opath, n ) ), is_true()))
}

context('plotDensity')

plotDensity( data )
try(expect_that( file.exists( file.path( opath, 'densityWebGL', 'index.html' ) ), is_true()))

if ( ! file.exists(file.path("..","..","data","analyzed.RData") ) ){
	save(data, file=file.path("..","..","data","analyzed.RData") )
	
}





