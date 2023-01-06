
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BipartiteModularityMaximization

The goal of this R package is to partition a bipartite network into
non-overlapping biclusters by maximizing bipartite modularity defined in
[Barber (2007)](https://doi.org/10.1103/PhysRevE.76.066102) using the
bipartite version of the algorithm described in [Treviño
(2015)](https://doi.org/10.1088/1742-5468/2015/02/P02003).

## Installation

### Environment

The package contains C/C++ code that needs compilation.

-   For Windows users, please download and install
    [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html),
    and add it to `PATH` by the follwing command:

``` r
write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
```

-   Mac users can install Xcode from the Mac AppStore.

-   Linux users can install a compiler and various development libraries
    (details vary across different flavors of Linux).

More details can be found in documents of
[RStudio](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites)
or [R](https://cran.r-project.org/doc/manuals/r-devel/R-admin.pdf).

### Installation from CRAN

Use the following command to install the released version of the package
from
[CRAN](https://cran.r-project.org/package=BipartiteModularityMaximization),
where pre-compiled “binary packages” might be available.

``` r
install.packages("BipartiteModularityMaximization")
```

### Installation from GitHub

Once you set up the compiler, you can install the latest version of the
package from
[GitHub](https://github.com/DIVA-Lab-UTMB/BipartiteModularityMaximization)
with:

``` r
install.packages("remotes")
remotes::install_github("DIVA-Lab-UTMB/BipartiteModularityMaximization")
```

## Example

This is a basic example which shows you how to use the main function
`bipmod` (short for bipartite modularity) to partition the example
bipartite network (represented as an incidence matrix of 798 rows and 8
columns):

``` r
## basic example code
library(BipartiteModularityMaximization)
data(example_data)
str(example_data)
#> 'data.frame':    798 obs. of  8 variables:
#>  $ Symptom_1: int  1 1 0 1 1 1 0 1 0 0 ...
#>  $ Symptom_2: int  1 1 0 1 1 1 0 1 1 0 ...
#>  $ Symptom_3: int  1 1 0 1 0 0 0 0 1 0 ...
#>  $ Symptom_4: int  1 1 0 1 1 1 1 1 1 0 ...
#>  $ Symptom_5: int  1 1 1 1 0 1 0 1 1 1 ...
#>  $ Symptom_6: int  0 1 0 1 1 1 0 1 1 0 ...
#>  $ Symptom_7: int  1 1 0 0 1 1 1 1 1 0 ...
#>  $ Symptom_8: int  1 0 1 1 1 1 0 1 1 0 ...
Q_part=bipmod(example_data)
str(Q_part)
#> List of 2
#>  $ MODULARITY: num 0.262
#>  $ ASSIGN    : int [1:806] 2 2 2 4 4 2 3 2 4 2 ...
```

## Documentation

Please read the documentation using `?bipmod` or `?example_data` for
more details.

## Related library

The biclusters can be visualized using
[ExplodeLayout](https://github.com/DIVA-Lab-UTMB/ExplodeLayout) or
[epl](https://github.com/UTMB-DIVA-Lab/epl) described in [Bhavnani
(2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5543384/pdf/2613038.pdf).
A more convenient wrapper of
[BipartiteModularityMaximization](https://github.com/DIVA-Lab-UTMB/BipartiteModularityMaximization)
and [ExplodeLayout](https://github.com/DIVA-Lab-UTMB/ExplodeLayout) can
be found in
[UtilitiesDIVA](https://github.com/DIVA-Lab-UTMB/UtilitiesDIVA).
