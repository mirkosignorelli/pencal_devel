# `R` package `pencal` (development version)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/pencal)](https://cran.r-project.org/package=pencal)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/pencal?color=orange)](http://cranlogs.r-pkg.org/badges/grand-total/pencal)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/pencal)](http://cranlogs.r-pkg.org/badges/pencal)
[![R-CMD-check](https://github.com/mirkosignorelli/pencal_devel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mirkosignorelli/pencal_devel/actions/workflows/R-CMD-check.yaml)

<img src="https://user-images.githubusercontent.com/20061736/162180793-072613f0-a93e-4ef6-b0c4-b8d8a45d770a.png" align="right" alt="" width="250" />

## `pencal`: what's that?

`pencal` is the `R` package that implements Penalized Regression Calibration (PRC), a statistical method that we proposed in Signorelli *et al.* (2021). [Penalized regression calibration: A method for the prediction of survival outcomes using complex longitudinal and high-dimensional data](https://onlinelibrary.wiley.com/doi/10.1002/sim.9178). *Statistics in Medicine*, 40 (27), 6178-6196. A detailed description of how to use PRC to do dynamic prediction of survival can be found in Signorelli (in review). [`pencal`: an `R` Package for the Dynamic Prediction of Survival with Many Longitudinal Predictors](https://cran.r-project.org/web/packages/pencal/vignettes/vignette.pdf).

## About this repository

This repository contains the **development** version of the `R` package. The latest CRAN version can be found [here](https://cran.r-project.org/web/packages/pencal/index.html).

## Installing `pencal`

### CRAN version (recommended)

To install the latest `CRAN` version of the package you may use:

```
install.packages('pencal')
```

If you encounter problems with packages on which `pencal` depends, you may alternatively install the package using `BiocManager`:

```
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install('pencal')
```

### Development version

You can install the development version of pencal from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mirkosignorelli/pencal_devel")
```

## Further information

More information on PRC can be found at the following pages:
* [CRAN package page](https://cran.r-project.org/web/packages/pencal/index.html);
* [my personal website](https://mirkosignorelli.github.io/r.html);
* [vignette that illustrates with examples how to use `pencal`](https://cran.r-project.org/web/packages/pencal/vignettes/vignette.pdf).
