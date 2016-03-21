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
				annotation='data.frame',
				snorm='logical',
				normFunct='character',
				outpath='character',
				name='character',
				rownamescol='character',
				sampleNamesCol='character',
				zscored = 'logical',
				simple = 'character',
				usedObj = 'list'
		),
		prototype(outpath ='', name = 'Rscexv',
				sampleNamesCol=NA_character_,
				facs=data.frame(),
				wFACS=F,
				snorm=F,
				normFunct='none',
				zscored=F,
				usedObj= list(),
				simple= c( 'outpath', 'rownamescol', 'sampleNamesCol', 'simple', 'snorm', 'zscored') )
)


#' @name show
#' @aliases show,Rscexv-method
#' @rdname show-methods
#' @docType methods
#' @description  print the Rscexv
#' @param object the Rscexv object
#' @return nothing
#' @title description of function show
#' @export 
setMethod('show', signature(object='Rscexv') ,
		definition = function (object) {
			cat (paste("An object of class", class(object)),"\n" )
			cat("named ",object@name,"\n")
			cat (paste( 'with',nrow(object@data),'samples and', ncol(object@data),' genes.'),"\n")
			if ( object@wFACS ) {
				cat ( paste("The object also contains FACS for",ncol(object@facs),'fluorocromes',sep=' '),"\n")
			}
			cat (paste("Annotation datasets (",paste(dim(object@annotation),collapse=','),"): '",paste( colnames(object@annotation ), collapse="', '"),"'  ",sep='' ),"\n")
			cat (paste("Sample annotation (",paste(dim(object@samples),collapse=','),"): '",paste( colnames(object@samples ), collapse="', '"),"'  ",sep='' ),"\n")
		})