# Rscexv

This package is meant to replace the R code in the SCExV server in the next update.

During the next update of the SCExV server an upgrade to the grouping is planned, which most likely will break the existing code. Therefore the rigid testing available with an R::S4 class might become important.

This package is under heavy development at the moment and should not be used!

# Requirements

R and several R packages, that should all be installed when installing this package.
All, but the MAST dependency, which is also hosted on GitHub:

library(devtools)

install_github('RGLab/MAST')

# Installation

Installation is extremely simple:

library(devtools)

install_github('stela2502/Rscexv')


# Usage

Please take a look at the test scripts and especially the SCExV_Functions.R file.