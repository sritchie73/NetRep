## ---- echo=FALSE, cache=FALSE---------------------------------------------------------------------
options(width = 100)

## -------------------------------------------------------------------------------------------------
library("NetRep")
data("NetRep")

## ---- message=FALSE, echo=FALSE-------------------------------------------------------------------
# This is hidden so the script will run. The other code blocks are demonstrative:
# they show the user how to use the 'backingfile' argument. This is primarily for
# people who are going to copy code from the vignette and modify it to their 
# data (e.g. most people), however the vignette should not save things to the
# user's current directory when they build the vignette 
discovery_data <- as.bigMatrix(discovery_data)
discovery_correlation <- as.bigMatrix(discovery_correlation)
discovery_network <- as.bigMatrix(discovery_network)
test_data <- as.bigMatrix(test_data)
test_correlation <- as.bigMatrix(test_correlation)
test_network <- as.bigMatrix(test_network)

## ---- eval=FALSE----------------------------------------------------------------------------------
#  # Save data in the 'bigMatrix' format in your current directory.
#  save.as.bigMatrix(discovery_data, backingfile="discovery_data.bm")
#  save.as.bigMatrix(discovery_correlation, backingfile="discovery_correlation.bm")
#  save.as.bigMatrix(discovery_network, backingfile="discovery_network.bm")
#  save.as.bigMatrix(test_data, backingfile="test_data.bm")
#  save.as.bigMatrix(test_correlation, backingfile="test_correlation.bm")
#  save.as.bigMatrix(test_network, backingfile="test_network.bm")
#  
#  # Write out the module assignments vector to file to be read in by the module
#  # preservation script we write later.
#  write.csv(module_labels, file="discovery_modules.csv")

## ---- eval=FALSE----------------------------------------------------------------------------------
#  discovery_network <- load.bigMatrix(backingfile="discovery_network.bm")

## -------------------------------------------------------------------------------------------------
discovery_network[1:5, 1:5]

## ---- eval=FALSE----------------------------------------------------------------------------------
#  # This converts a `bigMatrix` to a `matrix` in R, but leaves the
#  # backingfile on disk so you can still instantly load it in other
#  # R sessions using 'load.bigMatrix'
#  as.matrix(discovery_network)

## ---- eval=FALSE----------------------------------------------------------------------------------
#  # First, we need to load in the data we previously saved in
#  # the `bigMatrix` format:
#  discovery_data <- load.bigMatrix("discovery_data.bm")
#  discovery_correlation <- load.bigMatrix("discovery_correlation.bm")
#  discovery_network <- load.bigMatrix("discovery_network.bm")
#  test_data <- load.bigMatrix("test_data.bm")
#  test_correlation <- load.bigMatrix("test_correlation.bm")
#  test_network <- load.bigMatrix("test_network.bm")
#  
#  # As well as read in the module labels:
#  module_labels <- read.csv("module_labels.csv", stringsAsFactors=FALSE)
#  # Convert the 'data.frame' to a 'vector'
#  module_labels <- structure(module_labels[,2], names=module_labels[,1])

## ---- eval=FALSE----------------------------------------------------------------------------------
#  # Set up the input data structures for NetRep.
#  data_list <- list(cohort1=discovery_data, cohort2=test_data)
#  correlation_list <- list(cohort1=discovery_correlation, cohort2=test_correlation)
#  network_list <- list(cohort1=discovery_network, cohort2=test_network)
#  
#  # We do not need to set up a list for containing the 'module_labels' because
#  # there is only one "discovery" dataset.

## ---- eval=FALSE----------------------------------------------------------------------------------
#  # Assess the preservation of modules in the test dataset
#  preservation <- modulePreservation(
#   data=data_list, correlation=correlation_list, network=network_list,
#   moduleAssignments=module_labels, nPerm=10000, discovery="cohort1",
#   test="cohort2"
#  )
#  
#  # Write out the results object:
#  saveRDS(preservation, "preservation-analysis-results.rds")

## ---- echo=FALSE, hold=TRUE-----------------------------------------------------------------------
# This is the code that actually gets run in the Rmarkdown document
data_list <- list(cohort1=discovery_data, cohort2=test_data)
correlation_list <- list(cohort1=discovery_correlation, cohort2=test_correlation)
network_list <- list(cohort1=discovery_network, cohort2=test_network)

preservation <- modulePreservation(
 data=data_list, correlation=correlation_list, network=network_list,
 moduleAssignments=module_labels, nPerm=100, discovery="cohort1", 
 test="cohort2", verbose=TRUE
)

## ---- eval=FALSE----------------------------------------------------------------------------------
#  preservation <- readRDS("preservation-analysis-results.rds")
#  
#  # The results are stored as a list. The table of permutation test p-values is
#  # stored in the element named "p.value".
#  preservation$p.value

## ---- echo=FALSE----------------------------------------------------------------------------------
# This is the code that actually gets run in the Rmarkdown document
preservation$p.value

## -------------------------------------------------------------------------------------------------
# Get the maximum permutation test p-value
max_pval <- apply(preservation$p.value, 1, max)
max_pval

## ---- eval=FALSE----------------------------------------------------------------------------------
#  # First, we need to load in the data we previously saved in
#  # the `bigMatrix` format:
#  discovery_data <- load.bigMatrix("discovery_data.bm")
#  discovery_correlation <- load.bigMatrix("discovery_correlation.bm")
#  discovery_network <- load.bigMatrix("discovery_network.bm")
#  test_data <- load.bigMatrix("test_data.bm")
#  test_correlation <- load.bigMatrix("test_correlation.bm")
#  test_network <- load.bigMatrix("test_network.bm")
#  
#  # As well as read in the module labels:
#  module_labels <- read.csv("module_labels.csv", stringsAsFactors=FALSE)
#  # Convert the 'data.frame' to a 'vector'
#  module_labels <- structure(module_labels[,2], names=module_labels[,1])
#  
#  # Set up the input data structures for NetRep.
#  data_list <- list(cohort1=discovery_data, cohort2=test_data)
#  correlation_list <- list(cohort1=discovery_correlation, cohort2=test_correlation)
#  network_list <- list(cohort1=discovery_network, cohort2=test_network)

## ----modules_in_discovery, dev="png", fig.height=6, fig.width=6, fig.align="center", results="hold", fig.keep="last", fig.show="hold"----
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,2,3,4),
  discovery="cohort1", test="cohort1"
)

## ----modules_in_test, dev="png", fig.height=6, fig.width=6, fig.align="center", results="hold", fig.keep="last", fig.show="hold"----
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,2,3,4),
  discovery="cohort1", test="cohort2"
)

## ----mean_degree, dev="png", fig.height=6, fig.width=6, fig.align="center", results="hold", fig.keep="last", fig.show="hold"----
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,4), # only the preserved modules
  discovery="cohort1", test="cohort2",
  orderNodesBy=c("cohort1", "cohort2") # this can be any number of datasets
)

## ----dry_run, dev="png", fig.height=6, fig.width=6, fig.align="center", results="hold", fig.keep="last", fig.show="hold"----
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,4),
  discovery="cohort1", test="cohort2",
  orderNodesBy=c("cohort1", "cohort2"),
  dryRun=TRUE
)

## ----dry_run_customised, dev="png", fig.height=6, fig.width=6, fig.align="center", results="hold", fig.keep="last", fig.show="hold"----
# Change the margins so the plot is more compressed. Alternatively we could 
# change the device window.
par(mar=c(3,10,3,10)) # bottom, left, top, right margin sizes
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,4),
  discovery="cohort1", test="cohort2",
  orderNodesBy=c("cohort1", "cohort2"),
  dryRun=TRUE, 
  # Title of the plot
  main = "Preserved modules", 
  # Use the maximum edge weight as the highest value instead of 1 in the
  # network heatmap
  netRange=NA,
  # Turn off the node and sample labels:
  plotNodeNames=FALSE, plotSampleNames=FALSE,
  # The distance from the bottom axis should the module labels be drawn:
  maxt.line=0,
  # The distance from the legend the legend titles should be drawn:
  legend.main.line=2
)

## ----mean_degree_customised, dev="png", fig.height=6, fig.width=6, fig.align="center", results="hold", fig.keep="last", fig.show="hold"----
par(mar=c(3,10,3,10)) 
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,4),
  discovery="cohort1", test="cohort2",
  orderNodesBy=c("cohort1", "cohort2"), main = "Preserved modules", 
  netRange=NA, plotNodeNames=FALSE, plotSampleNames=FALSE,
  maxt.line=0, legend.main.line=2
)

## ----correlation_heatmap, dev="png", fig.height=6, fig.width=6, fig.align="center", results="hold", fig.keep="last", fig.show="hold"----
par(mar=c(5,5,3,3)) 
plotCorrelation(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=0:4, discovery="cohort1",
  test="cohort1", symmetric=TRUE, orderModules=FALSE
)

## ---- eval=FALSE----------------------------------------------------------------------------------
#  # First, we need to load in the data we previously saved in
#  # the `bigMatrix` format:
#  discovery_data <- load.bigMatrix("discovery_data.bm")
#  discovery_correlation <- load.bigMatrix("discovery_correlation.bm")
#  discovery_network <- load.bigMatrix("discovery_network.bm")
#  test_data <- load.bigMatrix("test_data.bm")
#  test_correlation <- load.bigMatrix("test_correlation.bm")
#  test_network <- load.bigMatrix("test_network.bm")
#  
#  # As well as read in the module labels:
#  module_labels <- read.csv("module_labels.csv", stringsAsFactors=FALSE)
#  # Convert the 'data.frame' to a 'vector'
#  module_labels <- structure(module_labels[,2], names=module_labels[,1])
#  
#  # Set up the input data structures for NetRep.
#  data_list <- list(cohort1=discovery_data, cohort2=test_data)
#  correlation_list <- list(cohort1=discovery_correlation, cohort2=test_correlation)
#  network_list <- list(cohort1=discovery_network, cohort2=test_network)

## -------------------------------------------------------------------------------------------------
properties <- networkProperties(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, 
  # Only calculate for the reproducible modules
  modules=c(1,4),
  # what dataset were the modules identified in?
  discovery="cohort1", 
  # which datasets do we want to calculate their properties in?
  test=c("cohort1", "cohort2")
)

# The summary profile of module 1 in the discovery dataset:
properties[["cohort1"]][["1"]][["summary"]]
# Along with the proportion of variance in the module data explained by the 
# summary profile:
properties[["cohort1"]][["1"]][["coherence"]]

# The same information in the test dataset:
properties[["cohort2"]][["1"]][["summary"]]
properties[["cohort2"]][["1"]][["coherence"]]

