# Rscexv

This package is meant to replace the R code in the SCExV server in the next update.

During the next update of the SCExV server an upgrade to the grouping is planned, which most likely will break the existing code. Therefore the rigid testing available with an R::S4 class might become important.

This package is under heavy development at the moment and should not be used!

# Requirements

R and several R packages, that should all be installed when installing this package.
Only the two dependencies hosted on github need to be installed 'manualy':

library(devtools)

install_github('RGLab/MAST')
install_github('stela2502/RFclust.SGE')


Recently I added another MDS option (ZIFA https://github.com/epierson9/ZIFA).
It is a Python library which has to be installed in order to use the ZIFA mds option.

# Installation

Installation is extremely simple:

library(devtools)

install_github('stela2502/Rscexv')

## installation on a fresh Ubuntu 16.04


sudo apt-get install r-base r-base-html r-base-core libcurl4-openssl-dev libssl-dev libssh2-1-dev libx11-dev libglu1-mesa-dev libfreetype6-dev


### Start a superuser R

install.packages(c('httr','git2r', 'devtools','Rcpp') )

source("http://bioconductor.org/biocLite.R")
biocLite("RDRToolbox")
biocLite("Biobase", 'BiocGenerics')

library(devtools)
install_github('stela2502/RFclust.SGE')
install_github('RGLab/MAST')
install_github('stela2502/Rscexv')
 


# Usage

Please take a look at the test scripts and especially the SCExV_Functions.R file.
