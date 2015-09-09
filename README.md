# netrep
##### Fast permutation procedure for testing network module replication

An R package containing functions for assessing the replication/preservation of 
network topology for weighted gene coexpression network modules in one or more
independent datasets through permutation testing.

The main function for this package is `modulePreservation`. Other
useful functions include `networkProperties` for calculating the
topological properties of a module, and `plotModule` and other
plotting functions, e.g. `plotCoexpression` for visualising a
gene coexpression network module.

## Installation

The latest stable version of this package as well as its package dependencies are maintained 
in a lightweight CRAN-like repository at http://InouyeLab.github.io/miniCRAN and can
be installed using `install.packages`:

```{r}
install.packages("netrep", repos="http://InouyeLab.github.io/miniCRAN")
```

The latest developmental version of this package can be installed directly from github using
the devtools package:

```{r}
library(devtools)
install_github("InouyeLab/netrep")
```

## Software dependencies: Installation troubleshooting

`netrep` and its dependencies require several third party libraries to be
installed. If not found, installation of the package will fail.

The `g++` compiler is required for the `bigmemory` package to install, 
`gfortran` is required for the `statmod` package to install, and a `BLAS` 
library is required for `RcppArmadillo` to install.

BLAS libraries must be installed prior to the installation of R, otherwise R 
won't link to them correctly and `RcppArmadillo` will fail to install. You've 
encountered this error if `RcppArmadillo` partially compiles, but then throws 
an error about failing to link to `-llapack`. LAPACK libraries come bundled with
most BLAS libraries.

Operating specific instructions to follow. 

## Testing
To ensure the package has installed correctly and will run on your system, run the following:

```{r}
library(testthat)
test_package("netrep")
```

