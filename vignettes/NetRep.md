# NetRep
Scott Ritchie  
`r Sys.Date()`  



## Introduction

The *NetRep* package provides functions for assessing the reproducibility and
preservation of network modules across datasets. The preservation of network 
modules in a second dataset is quantified through a permutation procedure that 
tests whether topological properties are more similar in the test dataset than
expected by chance: \emph{i.e.} when calculated between a module in the 
*discovery* dataset and random groups of nodes in the *test* dataset. 

This type of analysis is suitable for networks that can be meaningfully inferred
from multiple datasets. These include gene coexpression networks,
protein-protein interaction networks, and microbial interaction networks.
Modules within these networks consist of groups of nodes that are particularly
interesting: for example a group of tightly connected genes associated with a
disease, groups of genes annotated with the same term in the Gene Ontology
database, or groups of interacting microbial species, i.e. communities.
Application of this method can answer questions such as; (1) do the
relationships between genes in a module replicate in an independent cohort? (2)
are these gene coexpression modules preserved across tissues or tissue specific?
(3) are these modules conserved across species? (4) are microbial communities
preseved across multiple spatial locations?

In this tutorial, we will learn how to:

 1. Assess module preservation using the `modulePreservation` function.
 2. Visualise a module's network topology in both the discovery and test 
    datasets using the `plotModule` function.
 3. Calculate the topological properties of a network module for other 
    downstream analyses using the `networkProperties` function.
    
## Tutorial data

For this tutorial, we will use gene expression data simulated for two independent
cohorts. This data is provided with the package to demonstrate function usage: 


```r
library("NetRep")
data("NetRep")
ls()
```

```
## [1] "discovery_correlation" "discovery_data"        "discovery_network"     "module_labels"        
## [5] "oldLC"                 "test_correlation"      "test_data"             "test_network"
```

This loads seven objects into the R session:

  - `discovery_data`: a matrix with 150 columns (genes) and 30 rows (samples) 
     whose entries correspond to the expression level of each gene in each 
     sample in the discovery dataset.
  - `discovery_correlation`: a matrix with 150 columns and 150 rows containing 
     the correlation-coefficients between each pair of genes calculated from the 
    `discovery_data` matrix.
  - `discovery_network`: a matrix with 150 columns and 150 rows containing the 
     network edge weights encoding the interaction strength between each pair of 
     genes in the discovery dataset.
  - `module_labels`: a named vector with 150 entries containing the module 
     assignment for each gene as identified in the discovery dataset.
  - `test_data`: a matrix with 150 columns (genes) and 30 rows (samples) whose 
     entries correspond to the expression level of each gene in each sample in 
     the test dataset.
  - `test_correlation`: a matrix with 150 columns and 150 rows containing the 
     correlation-coefficients between each pair of genes calculated from the 
    `test_data` matrix.
  - `test_network`: a matrix with 150 columns and 150 rows containing the 
     network edge weights encoding the interaction strength between each pair of
     genes in the test dataset.
     
The *discovery* dataset was simulated to contain four modules of varying size,
two of which (Modules 1 and 4) replicate in the *test* dataset. Details of the
simulation are provided in the documentation for the package data 
(see `help("NetRep-data")`).

## Setting up the input data 

The topological properties used to quantify module preservation are calculated 
not only from the interaction networks inferred in each dataset, but also from 
the data used to infer those networks (e.g. gene expression data) as well as 
the correlation structure between variables/nodes. 

All functions in the *NetRep* package have the following arguments:

 - `network`: a list of interaction networks, one for each dataset.
 - `data`: a list of data matrices used to infer those networks, one for each 
    dataset.
 - `correlation`: a list of matrices containing the pairwise correlation 
    coefficients between variables/nodes in each dataset.
 - `moduleAssignments`: a list of vectors, one for each *discovery* dataset, 
    containing the module assignments for each node in that dataset.
 - `modules`: a list of vectors, one vector for each *discovery* dataset, 
    containing the names of the modules from that dataset to analyse.
 - `discovery`: a vector indicating the names or indices to use as the 
   *discovery* datasets in the `network`, `data`, `correlation`, 
   `moduleAssignments`, and `modules` arguments.
 - `test`: a list of vectors, one vector for each *discovery* dataset, 
    containing the names or indices of the `network`, `data`, and `correlation` 
    argument lists to use as the *test* dataset(s) for the analysis of each 
   *discovery* dataset.

*NetRep* requires the matrix data for each of the `network`, `data`, and 
`correlation` arguments to be stored in the `bigMatrix` format. This is a class
provided by *NetRep* that stores the data in shared memory, which allows the 
data to be accessed from multiple parallel R sessions (see `help("bigMatrix")` 
for more details).

First, we will convert the tutorial data into the `bigMatrix` format:


```r
# Convert the data to the 'bigMatrix' format:
discovery_data <- as.bigMatrix(discovery_data)
discovery_correlation <- as.bigMatrix(discovery_correlation)
discovery_network <- as.bigMatrix(discovery_network)
test_data <- as.bigMatrix(test_data)
test_correlation <- as.bigMatrix(test_correlation)
test_network <- as.bigMatrix(test_network)
```

If you are using your own data, you can convert it either using `as.bigMatrix`,
or by reading it in directly using `read.bigMatrix`. 

Next, we will set up the input lists:


```r
data_list <- list(discovery=discovery_data, test=test_data)
correlation_list <- list(discovery=discovery_correlation, test=test_correlation)
network_list <- list(discovery=discovery_network, test=test_network)
```

Since we only have one discovery and one test dataset we do not have to set up
list structures for the other common arguments.

## Assessing module preservation

Now we can use the `modulePreservation` function to assess the reproducibility 
of the four simulated modules in the test dataset. We recommend running the 
procedure with at least 10,000 permutations to ensure the permutation procedure
generates representative null distributions for each statistic.

*NetRep* will automatically detect the number of cores and use all but one. 
Alternatively, the number of cores can be set by the user using the `nCores`
argument. For a full list of arguments see `help("modulePreservation")`. On a 
single core machine the following should take an hour to an hour and a half 
to run. On a laptop with 8 cores it took just under 10 minutes when parallelised
over seven cores.

**Note:** scalable memory usage depends on clean R sessions. All objects in the 
R session are copied to each parallel core, with the exception of the
`bigMatrix` objects which are stored in shared memory. Provided there are no
additional objects in the R session only the memory required to store each 
`bigMatrix` in memory once will be used, along with approximately additional 200
MB per core used by each vanilla R session.


```r
# NetRep will assess module preservation for *all* modules by default
preservation <- modulePreservation(
 data=data_list, correlation=correlation_list, network=network_list,
 moduleAssignments=module_labels, nPerm=10000, discovery="discovery", 
 test="test"
)
```

```
##  Validating user input...
##    Running on 7 cores.
##    Checking matrices for non-finite values...
##  User input ok!
##  Calculating preservation of network subsets from dataset "discovery" in dataset "test".
##    Calculating observed test statistics...
##    Calculating null distributions with 10000 permutations...
##    Calculating P-values...
##    Collating results...
##  Cleaning up temporary objects...
##  Done!
```

Now we can look at the observed values and permutation P-values for each module 
preservation statistic:


```r
# Each row corresponds to a module
preservation$observed
```

```
##    avg.weight coherence    cor.cor  cor.degree cor.contrib      avg.cor avg.contrib
## 1 0.161069393 0.6187688 0.78448573  0.90843993   0.8795006  0.550004272  0.76084777
## 2 0.001872928 0.1359063 0.17270312 -0.03542772   0.5390504  0.034040922  0.23124826
## 3 0.001957475 0.1263280 0.01121223 -0.17179855  -0.1074944 -0.007631867  0.05412794
## 4 0.046291489 0.4871179 0.32610667  0.68122446   0.5251965  0.442614173  0.68239136
```

```r
preservation$p.values
```

```
##   avg.weight  coherence    cor.cor cor.degree cor.contrib    avg.cor avg.contrib
## 1 0.00009999 0.00009999 0.00009999 0.00009999  0.00009999 0.00009999  0.00009999
## 2 0.97710229 0.97100290 0.01019898 0.56204380  0.00359964 0.01799820  0.00659934
## 3 0.98890111 0.98550145 0.41805819 0.81081892  0.71152885 0.99420058  0.88131187
## 4 0.00009999 0.00009999 0.00009999 0.00009999  0.00059994 0.00009999  0.00009999
```

Details for each statistic are provided in the documentation (see 
`help("modulePreservation")`).

Here, we will consider a module to be reprocible if all statistics have P < 0.01:


```r
max_pval <- apply(preservation$p.value, 1, max)
max_pval
```

```
##          1          2          3          4 
## 0.00009999 0.97710229 0.99420058 0.00059994
```

Only modules 1 and 4 are reproducible at this significance threshold.

## Visualising network modules

The topological properties that contribute to each module preservation statistic
can be visualised using `plotModule`. First, let's look at the four module in
the *discovery* dataset:


```r
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,2,3,4),
  discovery="discovery", test="discovery"
)
```

```
##  Validating user input...
##    Running on 7 cores.
##    Checking matrices for non-finite values...
##  User input ok!
##  Ordering nodes...
##  Ordering samples...
##  rendering plot components...
```

<img src="NetRep_files/figure-html/unnamed-chunk-8-1.png" title="" alt="" style="display: block; margin: auto;" />

```
##  Cleaning up temporary objects...
##  Done!
```

The `discovery` argument tells `plotModule` which dataset the modules were 
identified in, while the `test` argument tells `plotModule` which dataset to
calculate those modules' properties in. By default, nodes are ordered from left
to right in decreasing order of *weighted degree* and modules are ordered by the
similarity of their summary profiles. For visualisation, the *weighted degree* 
is normalised within each module by the maximum value since the *weighted degree*
of nodes can be dramatically different for modules of different sizes.

If we set `test = "test"`, then `plotModule` will show the properties in the 
*test* dataset while keeping the ordering of nodes and modules the same as if 
calculated in the `discovery` dataset:


```r
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,2,3,4),
  discovery="discovery", test="test"
)
```

```
##  Validating user input...
##    Running on 7 cores.
##    Checking matrices for non-finite values...
##  User input ok!
##  Ordering nodes...
##  Ordering samples...
##  rendering plot components...
```

<img src="NetRep_files/figure-html/unnamed-chunk-9-1.png" title="" alt="" style="display: block; margin: auto;" />

```
##  Cleaning up temporary objects...
##  Done!
```

There are many arguments to adjust the way modules are drawn (see 
`help("plotModule"))`). For example, we can change the previous plot so that 
nodes are ordered within the `test` dataset and modules are drawn in the order 
provided rather than by the similarity of their summary profiles:


```r
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,2,3,4),
  discovery="discovery", test="test", orderNodesBy="test", orderModules=FALSE
)
```

```
##  Validating user input...
##    Running on 7 cores.
##    Checking matrices for non-finite values...
##  User input ok!
##  Ordering nodes...
##  Ordering samples...
##  rendering plot components...
```

<img src="NetRep_files/figure-html/unnamed-chunk-10-1.png" title="" alt="" style="display: block; margin: auto;" />

```
##  Cleaning up temporary objects...
##  Done!
```

We can also plot individual components of the plot separately. For example, 
a heatmap of the correlation structure:


```r
par(mar=c(6,6,6,4)+0.1)
plotCorrelation(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=0:4, discovery="discovery",
  test="discovery", symmetric=TRUE, orderModules=FALSE
)
```

```
##  Validating user input...
##    Running on 7 cores.
##    Checking matrices for non-finite values...
##  User input ok!
##  Ordering nodes...
##  rendering plot components...
```

<img src="NetRep_files/figure-html/unnamed-chunk-11-1.png" title="" alt="" style="display: block; margin: auto;" />

```
##  Cleaning up temporary objects...
##  Done!
```

A full list of function and arguments for these individual plots can be found
at `help("plotTopology")`.

## Calculating the network properties

Finally, the network properties for each module can be calculated in both 
datasets through the `networkProperties` function:


```r
properties <- networkProperties(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, 
  # Only calculate for the reproducible modules
  modules=c(1,4),
  # what dataset were the modules identified in?
  discovery="discovery", 
  # which datasets do we want to calculate their properties in?
  test=c("discovery", "test")
)
```

```
##  Validating user input...
##    Running on 7 cores.
##    Checking matrices for non-finite values...
##  User input ok!
##  Calculating properties for:
##  Cleaning up temporary objects...
##  Done!
```

These properties can be useful for downstream analysis. For example, the 
module summary profiles can be used to assess associations between modules and
external information, e.g. case or control status:


```r
# The summary profile of module 1 in the discovery dataset:
properties[["discovery"]][["1"]][["summary"]]
```

```
##  Discovery_1  Discovery_2  Discovery_3  Discovery_4  Discovery_5  Discovery_6  Discovery_7 
##  -0.15173019  -0.09817810  -0.10356266  -0.21351111  -0.06424053  -0.25787365  -0.06191222 
##  Discovery_8  Discovery_9 Discovery_10 Discovery_11 Discovery_12 Discovery_13 Discovery_14 
##  -0.05886898   0.04544493   0.16790065  -0.16163254  -0.07158769  -0.16775343   0.39457572 
## Discovery_15 Discovery_16 Discovery_17 Discovery_18 Discovery_19 Discovery_20 Discovery_21 
##   0.10762551   0.25872801   0.01187731   0.57266243   0.15737963   0.02368060  -0.07088476 
## Discovery_22 Discovery_23 Discovery_24 Discovery_25 Discovery_26 Discovery_27 Discovery_28 
##   0.03726126  -0.13770047  -0.01978039  -0.06336512  -0.06360727  -0.30044215   0.14682841 
## Discovery_29 Discovery_30 
##   0.07036710   0.07229971
```

```r
# Along with the proportion of variance in the module data explained by the 
# summary profile:
properties[["discovery"]][["1"]][["coherence"]]
```

```
## [1] 0.585781
```

```r
# The same information in the test dataset:
properties[["test"]][["1"]][["summary"]]
```

```
##       Test_1       Test_2       Test_3       Test_4       Test_5       Test_6       Test_7 
## -0.099957918  0.061501299  0.043541623  0.051055323  0.056572949  0.136605203  0.116491092 
##       Test_8       Test_9      Test_10      Test_11      Test_12      Test_13      Test_14 
## -0.395294200 -0.099564626  0.092715774 -0.005526985  0.256963062  0.028746029 -0.076793357 
##      Test_15      Test_16      Test_17      Test_18      Test_19      Test_20      Test_21 
## -0.435677499  0.100475978 -0.339161521 -0.195830382 -0.104643904  0.050046780  0.238180614 
##      Test_22      Test_23      Test_24      Test_25      Test_26      Test_27      Test_28 
##  0.144114251  0.211841029  0.228291634 -0.171340087 -0.188991911 -0.093239829  0.063972325 
##      Test_29      Test_30 
##  0.278339356  0.046567899
```

```r
properties[["test"]][["1"]][["coherence"]]
```

```
## [1] 0.6187688
```
