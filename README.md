# MSnID

<!-- badges: start -->
[![R-CMD-check](https://github.com/PNNL-Comp-Mass-Spec/MSnID/workflows/R-CMD-check/badge.svg)](https://github.com/PNNL-Comp-Mass-Spec/MSnID/actions)
<!-- badges: end -->

`MSnID` is an R/Bioconductor package that provides convenient tools for handling and filtering of MS/MS identification results. The official page is the Bioconductor landing page ([release](http://www.bioconductor.org/packages/release/bioc/html/MSnID.html) or [devel](http://www.bioconductor.org/packages/devel/bioc/html/MSnID.html) versions). The [github page](https://github.com/PNNL-Comp-Mass-Spec/MSnID) is for sharing, testing, issue tracking and forking/pulling purposes.


## R Installation and Usage

```R
# Option 1: Installing from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MSnID")

# Option 2: Installing from GitHub
if (!require("remotes", quietly = TRUE)) 
    install.packages("remotes")
remotes::install_github("PNNL-Comp-Mass-Spec/MSnID@pnnl-master")

library(MSnID)
```

By default, `install_github` does not build vignettes. Specify `build_vignettes = TRUE`.


## Original Location

The original location is on the [vladpetyuk](https://github.com/vladpetyuk) account, repo [MSnID](https://github.com/vladpetyuk/MSnID).

