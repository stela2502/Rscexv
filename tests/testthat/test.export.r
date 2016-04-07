data(InbuiltData)
InbuiltData@outpath <- tempdir()

saveObj(InbuiltData)

expect_that( file.exists( file.path(InbuiltData@outpath, 'SCExV_Grps.txt')), is_true())
expect_that( file.exists( file.path(InbuiltData@outpath, 'norm_data.RData')), is_true())
expect_that( file.exists( file.path(InbuiltData@outpath, 'analysis.RData')), is_false())


expect_that( read.delim( file=file.path(InbuiltData@outpath, 'SCExV_Grps.txt'), header=F), equals ( data.frame( V1=c('none', 'Group by plateID') ) ) )

data <- group_1D (InbuiltData, 'Gata1', c( 3,5,10,15,20,30 ) )

saveObj(data )
