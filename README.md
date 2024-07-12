# fastTest
Fast Rcpp/C++ translations of popular statistical hypothesis tests in R

## General motivation

Undoubtedly, R is a great analysis platform, but memory handling and speed is of great conern, especially in a 'big' data setting such as transcriptomics, proteomics, or public health data - just to name few domains of my interest. While the speed limitations may be partially addressed by [parallelization](https://furrr.futureverse.org/), the R interface to C++, [Rcpp](https://www.rcpp.org/) provides the most elegant solution. The `fastTest` package bundles a growing set of statistical hypothesis testing tools, which are, except for ANOVA and correlation tests, almost 1:1 C++ translation of genuine R code and offer approximately 3- to 10-fold faster computation with minimal memory burden. The secondary motivation was to generate statistical hypothesis testing tools which return [effect size statistics](https://en.wikipedia.org/wiki/Effect_size) as well, at the moment the testing and effect size tools live in base R, stats and various packages with quite often incompatible interfaces. 'Standard' R functions and their `fastTest` counterparts are shown in the scheme below:

![package_usage](https://github.com/user-attachments/assets/a9877df3-5e2b-4669-9feb-c096800b3232)


## Installation

You may fetch the package with `devtools`: 

```r

devtools::install_github('PiotrTymoszuk/fastTest')

```
## Acknowledgements

Credits to creators, contributors and maintainers of [Rcpp](https://www.rcpp.org/), [tidyverse](https://www.tidyverse.org/), and [furrr](https://furrr.futureverse.org/).


## Terms of use

The package is available under a [GPL-3 license](https://github.com/PiotrTymoszuk/fastTest/blob/main/LICENSE).

## Contact

The package maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).
