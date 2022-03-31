
<!-- README.md is generated from README.Rmd. Please edit that file -->
BipartiteModularityMaximization
===============================

The goal of this R package is to partition a bipartite network into non-overlapping biclusters by maximizing bipartite modularity defined in Barber (2007) <doi:10.1103/PhysRevE.76.066102> using the bipartite version of the algorithm described in Treviño (2015) <doi:10.1088/1742-5468/2015/02/p02003>.

Installation
------------

You can install the released version of BipartiteModularityMaximization from [GitHub](https://github.com/DIVA-Lab-UTMB/BipartiteModularityMaximization) with:

``` r
install.packages("devtools")
devtools::install_github("UTMB-DIVA-Lab/BipartiteModularityMaximization")
```

Example
-------

This is a basic example which shows you how to use the main function `bipmod` (short for bipartite modularity) to partition the example bipartite network (represented as an incidence matrix of 798 rows and 8 columns):

``` r
## basic example code
library(BipartiteModularityMaximization)
data(example_data)
str(example_data)
#> 'data.frame':    798 obs. of  8 variables:
#>  $ not.happy          : int  0 1 0 0 0 1 1 1 0 0 ...
#>  $ didnot.enjoyed.life: int  0 1 0 0 0 0 0 1 0 0 ...
#>  $ felt.depressed     : int  0 1 1 0 0 1 0 1 0 1 ...
#>  $ everything.effort  : int  1 0 1 1 0 1 1 1 0 1 ...
#>  $ felt.lonely        : int  0 1 0 0 0 0 1 1 0 1 ...
#>  $ felt.sad           : int  0 1 0 0 0 1 0 1 0 1 ...
#>  $ couldnt.get.going  : int  1 0 0 1 0 1 1 1 0 0 ...
#>  $ sleep.restless     : int  0 1 0 0 1 0 0 1 1 1 ...
Q_part=bipmod(example_data)
str(Q_part)
#> List of 2
#>  $ MODULARITY: num 0.262
#>  $ ASSIGN    : int [1:806] 3 4 4 3 1 4 3 4 1 2 ...
```

Documentation
-------------

Please read the documentation using `?bipmod` or `?example_data` for more details.

Related library
---------------

The biclusters can be visualized using ExplodeLayout [epl](https://github.com/UTMB-DIVA-Lab/epl) described in <http://www.skbhavnani.com/DIVA/papers/Dang-et-al-AMIA-Demo-2016.pdf>.
