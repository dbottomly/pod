The `pod` package
=========

The `pod` package contains implementations of the methods compared in 'Comparison of methods to identify aberrant expression patterns in individual patients: augmenting our toolkit for precision medicine'.  Currently, this is a bare-bones package however, I will be updating it periodically with documentation, new methods etc. 

How to install
--------

Using a recent version of R (version >= 3.1) use the following commands:

source("http://bioconductor.org/biocLite.R")

biocLite(c("Biobase", "devtools"))

library(devtools)

install_github(rep="dbottomly/pod", ref="master") 


Basic usage
--------

See the vignette for typical usage:

library(pod)
vignette("pod")

How to contribute
---------

Please send me an email if you are interested in contributing.

Contact
---------

Dan Bottomly
bottomly@ohsu.edu
