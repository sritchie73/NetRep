# NetRep
##### Fast permutation procedure for testing network module replication

An R package containing functions for assessing the replication and 
preservation of a network module's topology across datasets through 
permutation testing. This is suitable for networks that can be meaningfully 
inferred from multiple datasets. These include gene coexpression networks,
protein-protein interaction networks, and microbial interaction
networks. Modules within these networks consist of groups of nodes
that are particularly interesting: for example a group of tightly
connected genes associated with a disease, groups of genes
annotated with the same term in the Gene Ontology database, or
groups of interacting microbial species, i.e. communities.
Application of this method can answer questions such as; (1) do
the relationships between genes in a module replicate in an
independent cohort? (2) are these gene coexpression modules
preserved across tissues or tissue specific? (3) are these modules
conserved across species? (4) are microbial communities preseved
across multiple spatial locations?

The main function for this package is `modulePreservation`, which 
performs the permutation test procedure. Other useful functions include 
`networkProperties` for calculating the topological properties of a 
module, and `plotModule` for visualising a network module.

For more information see the associated publication in Cell Systems,
[A Scalable Permutation Approach Reveals Replication and Preservation Patterns of Network Modules in Large Datasets](http://dx.doi.org/10.1016/j.cels.2016.06.012). 

## Installation

The latest stable version of NetRep can be installed either directly from
CRAN or from this GitHub repository:

```{r}
# From CRAN
install.packages("NetRep")

# Alternatively From GitHub
library(devtools)
install_github("InouyeLab/NetRep")
```

Developmental / pre-release versions of NetRep can be installed from this repository:

```{r}
library(devtools)
install_github("InouyeLab/NetRep", ref="devel")
```

Older versions of NetRep can be installed by specifying the version number in the `ref` argument:

```{r}
install_github("InouyeLab/NetRep", ref="v0.61.0")
```

## Package Tutorial

A vignette (tutorial) is available online at [vignettes/NetRep.md](vignettes/NetRep.md),
or can be loaded directly from R by running `vignette("NetRep")` if you have installed
the package from CRAN. 

If you are installing NetRep from GitHub and wish to make the vignette available on your 
local machine, you will need to install `rmarkdown` and specify `build_vignettes=TRUE` 
when running `install_github`.

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
