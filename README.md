# NetRep
##### Fast permutation procedure for testing network module replication

An R package containing functions for assessing the replication/preservation of 
network topology for weighted gene coexpression network modules in one or more
independent datasets through permutation testing.

A preprint is available on BioRxiv: [A scalable permutation approach reveals replication and preservation patterns of gene coexpression modules](http://biorxiv.org/content/early/2015/10/21/029553)

The main function for this package is `modulePreservation`. Other
useful functions include `networkProperties` for calculating the
topological properties of a module, and `plotModule` for visualising a
gene coexpression network module.

## Installation

The latest stable version of this package can be installed directly from this
Github repository:

```{r}
library(devtools)
install_github("InouyeLab/NetRep")
```

If you have `rmarkdown` installed, a package vignette (tutorial) is also 
available to install:

```{r}
install_github("InouyeLab/NetRep", build_vignettes=TRUE)
vignette("NetRep")
```
But can otherwise be read online at [vignettes/NetRep.md](vignettes/NetRep.md)

Older versions of NetRep can be installed by specifying the version number in the `ref` argument:

```{r}
install_github("InouyeLab/NetRep", ref="v0.61.0")
```
### Testing the package installation

To ensure the package has installed correctly and will run on your system, run the following:

```{r}
library(testthat)
test_package("NetRep")
```

## Installation troubleshooting

**NetRep** and its dependencies require several third party libraries to be
installed. If not found, installation of the package will fail.

 1. A compiler with `C++11` support for the `<thread>` libary.
 2. A `fortran` compiler.
 3. `BLAS` and `LAPACK` libraries.
 
### OSX

The necessary `fortran` and `C++11` compilers are provided with the `Xcode` 
application and subsequent installation of `Command line tools`. The most
recent version of OSX should prompt you to install these tools when 
installing the `devtools` package from RStudio. Those with older versions of 
OSX should be able to install these tools by typing the following command into 
their Terminal application: `xcode-select --install`.

Some users on OSX Mavericks have reported that even after this step they 
receive errors relating to `-lgfortran` or `-lquadmath`. This is reportedly 
solved by installing the version of `gfortran` used to compile the R binary for
your system: `gfortran-4.8.2`, using the following commands in your `Terminal` 
application

```{r, engine="bash", eval=FALSE}
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

### Windows

**NetRep** can be installed on Windows in R version 3.3.0 or later. The 
necessary `fortran` and `C++11` compilers are provided with the `Rtools` 
program. We recommend installation of `NetRep` through `RStudio`, which should
prompt the user and install these tools when running 
`devtools::install_github("InouyeLab/NetRep")`. You may need to run this 
command again after Rtools finishes installing.

### Linux

If installation fails when compiling **NetRep** at `permutations.cpp` with an 
error about `namespace thread`, you will need to install a newer version of 
your compiler that supports this `C++11` feature. We have found that this works
on versions of `gcc` as old as `gcc-4.6.3`.

If installation fails prior to this it is likely that you will need to install
the necessary compilers and libraries, then reinstall R. For `C++` and 
`fortran` compilers we recommend installing `g++` and `gfortran` from the
appropriate package manager for your operating system (e.g. `apt-get` for 
Ubuntu). `BLAS` and `LAPACK` libraries can be installed by installing 
`libblas-dev` and `liblapack-dev`. Note that these libraries **must** be
installed prior to installation of R.