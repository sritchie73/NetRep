# NetRep
##### Fast permutation procedure for testing network module replication

An R package containing functions for assessing the replication/preservation of 
network topology for weighted gene coexpression network modules in one or more
independent datasets through permutation testing.

A preprint is available on BioRxiv: [A scalable permutation approach reveals replication and preservation patterns of gene coexpression modules](http://biorxiv.org/content/early/2015/10/21/029553)

The main function for this package is `modulePreservation`. Other
useful functions include `networkProperties` for calculating the
topological properties of a module, and `plotModule` and other
plotting functions, e.g. `plotCoexpression` for visualising a
gene coexpression network module.

## Installation

The latest stable version of this package can be installed directly from this
Github repository:

```{r}
library(devtools)
install_github("InouyeLab/NetRep")
```

The latest developmental version of this package can be installed directly from github using
the devtools package:

```{r}
library(devtools)
install_github("InouyeLab/NetRep", branch="dev")
```

## Software dependencies: Installation troubleshooting

`NetRep` and its dependencies require several third party libraries to be
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
test_package("NetRep")
```
