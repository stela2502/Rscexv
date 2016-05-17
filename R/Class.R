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
			cat ( paste ( "A total of",length(object@usedObj),"Other objects have been added to this Rscexv object:\n", paste( collapse=', ', names(object@usedObj))),"\n")
		}
)


#' @name Rscexv
#' @aliases Rscexv,Rscexv-method
#' @rdname Rscexv-methods
#' @docType methods
#' @description create a new Rscexv object from the data files.
#' @param PCR  the pcr data file names default=NULL
#' @param FACS the FACS data file names MISSING default=NULL
#' @param use_pass_fail whether or not to use the pass_fail, information in the PCR data files.
#' @title description of function createDataObj
#' @export 
setGeneric('Rscexv', ## Name
		function  ( PCR=NULL,  FACS=NULL, use_pass_fail = T ){ ## Argumente der generischen Funktion
			standardGeneric('Rscexv') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('Rscexv', signature = c ('character'),
		definition = function ( PCR=NULL,  FACS=NULL, use_pass_fail = T ){
			
			fnamesPCR <- PCR
			fnamesFACS <- FACS
			PCR <- read.PCR.set ( PCR, use_pass_fail )
			
			system ( 'rm filtered*' )
			
			try ( system ( 'rm R_file_read_error.txt' ) )
			
			if ( is.null(PCR) ){
				system ( 'echo "You need to upload at least one PCR data set in order to start the analysis" > R.error')
				stop("You need to upload at least one PCR data set in order to start the analysis")
			}
			
			colnames(PCR$data) <- str_replace_all( colnames(PCR$data), '/', '_' )
			wFACS=F
			if ( ! is.null(FACS)){
				FACS <- read.FACS.set ( FACS)
				wFACS=T
				## black magick with the FACS gene names
				colnames(FACS$data) <- str_replace_all( colnames(FACS$data), '^P[0-9]+.', '' )
				colnames(FACS$data) <- str_replace_all( colnames(FACS$data), '.Min$', '' )
			}
			
			## check if the samples do overlapp
			
			
			data <- check.dataObj(PCR, FACS)
			
			data$samples$fnamesPCR <- fnamesPCR[data$samples$ArrayID]
			if ( wFACS ) {
				data$samples$fnamesFACS <- fnamesFACS[data$samples$ArrayID]
			}
			missing <- setdiff( 1:length(fnamesPCR), unique(data$samples$ArrayID) )
			if ( length(missing) > 0 ){
				stop( paste( length(missing),"files failed to load into the analysis:", paste( collapse= ', ', fnamesPCR[missing] ),". Most likely due to gene name missmatches." ))
			}
			## now create the object and done
			res <- new('Rscexv', data=data.frame(data$PCR), 
					facs=data.frame(data$FACS), samples=data$samples, 
					annotation=data$annotation, wFACS=wFACS, outpath=pwd(), baseSamplesCol=ncol(data$samples) )
			
			res
		}
)

#' @name pwd
#' @aliases pwd
#' @rdname pwd-methods
#' @docType methods
#' @description  uses the linux pwd command to determin the working directory 
#' @return A string containing the working directory 
#' @title description of function pwd
#' @export 
setGeneric('pwd', ## Name
		function ( a ) { ## Argumente der generischen Funktion
			standardGeneric('pwd') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('pwd', signature = c () ,
		definition = function ( a ) {
			rm(a)
			system( 'pwd > __pwd' )
			t <- read.delim( file = '__pwd', header=F)
			t <- as.vector(t[1,1])
			t <- paste(t,"/",sep='')
			unlink( '__pwd')
			t
		}
)
