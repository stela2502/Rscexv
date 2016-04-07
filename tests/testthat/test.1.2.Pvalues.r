data(analyzed)

context('Pvalues with FACS')

data@outpath = tempdir()

stat_obj <- create_p_values( 
		data, 
		boot = 1000,
		lin_lang_file= 'lin_lang_stats.xls' , 
		sca_ofile ="Significant_genes.csv"
)


try(expect_that( file.exists( file.path(stat_obj@outpath , 'lin_lang_stats.xls' ) ), is_true()))
try(expect_that( file.exists( file.path(stat_obj@outpath , 'Significant_genes.csv' ) ), is_true()))

