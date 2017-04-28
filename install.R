if (!library("devtools", quietly = TRUE,logical.return=TRUE )) {
	install.packages(c('devtools'),  repos='https://ftp.acc.umu.se/mirror/CRAN/')
	library(devtools)
}
source("https://bioconductor.org/biocLite.R")
biocLite()
if (!library("MAST", quietly = TRUE,logical.return=TRUE )) {
	install_github('RGLab/MAST')
}
if (!library("RFclust", quietly = TRUE,logical.return=TRUE )) {
	install_github('stela2502/RFclust.SGE')
}
install_github( 'RGLab/MAST', ref="MASTClassic" )
install()
