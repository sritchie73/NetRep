# Fast permutation procedure for testing network module replication

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
in a lightweight CRAN-like repository at https://sritchie73.github.io/miniCRAN and can
be installed using `install.packages`:

```{r}
install.packages("fastModPres", repos="http://sritchie73.github.io/miniCRAN")
```

The most up-to-date version of this package can be installed directly from github using
the devtools package:

```{r}
library(devtools)
install_github("sritchie73/fastModPres")
```

## Testing
To ensure the package has installed correctly and will run on your system, run the following:

```{r}
library(testthat)
test_package("fastModPres")
```

