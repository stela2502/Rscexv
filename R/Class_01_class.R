#' @name Rscexv
#' @title Rscexv
#' @description  This S4 class implements the SCExV backend logics
#' @slot data a data.frame containing the expression values for each gene x sample (gene = row)
#' @slot facs a data.frame containing all the FACS data 
#' @slot samples a data.frame describing the columnanmes in the data column
#' @slot annotation a data.frame describing the rownames in the data rows
#' @slot outpath the default outpath for the plots and tables from this package
#' @slot name the name for this package (all filesnames contain that)
#' @slot zscored genes are normalized?
#' @slot snorm samples normalized?
#' @slot normFunct which PCR norm function has been called
#' @slot rownamescol the column name in the annotation table that represents the rownames of the data table
#' @slot sampleNamesCol the column name in the samples table that represents the colnames of the data table
#' @slot usedObj here a set of used and probably lateron important objects can be saved. Be very carful using any of them!
#' @exportClass Rscexv
setClass(
		Class='Rscexv', 
		representation = representation ( 
				data='data.frame',
				facs='data.frame',
				wFACS='logical',
				samples='data.frame',
				raw="data.frame",
				snorm="data.frame",
				annotation='data.frame',
				norm='logical',
				normFunct='character',
				outpath='character',
				name='character',
				rownamescol='character',
				sampleNamesCol='character',
				zscored = 'logical',
				simple = 'character',
				baseSamplesCol='numeric',
				usedObj = 'list'
		),
		prototype(outpath ='', name = 'Rscexv',
				sampleNamesCol=NA_character_,
				facs=data.frame(),
				snorm=data.frame(),
				raw=data.frame(),
				wFACS=F,
				norm=F,
				normFunct='none',
				zscored=F,
				usedObj= list(),
				simple= c( 'outpath', 'rownamescol', 'sampleNamesCol', 'simple', 'zscored') )
)


#' @name InbuiltData
#' @title This is the test dataset also available at the SCExV server page
#' @description Please read PMID:26437766 for more details. No analysis performed.
#' @docType data
#' @usage data(InbuiltData)
#' @format Rscexv object
'InbuiltData'

#' @name analyzed
#' @title This is the analyzed data used for downstream testing
#' @description This file can be re-created running the test scripts.
#' @docType data
#' @usage data(analyzed)
#' @format Rscexv object
'data'



