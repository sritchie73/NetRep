## ---- echo=FALSE, cache=FALSE---------------------------------------------------------------------
options(width = 100)

## -------------------------------------------------------------------------------------------------
library("NetRep")
data("NetRep")
ls()

## ---- message=FALSE-------------------------------------------------------------------------------
# Convert the data to the 'bigMatrix' format:
discovery_data <- as.bigMatrix(discovery_data)
discovery_correlation <- as.bigMatrix(discovery_correlation)
discovery_network <- as.bigMatrix(discovery_network)
test_data <- as.bigMatrix(test_data)
test_correlation <- as.bigMatrix(test_correlation)
test_network <- as.bigMatrix(test_network)

## -------------------------------------------------------------------------------------------------
data_list <- list(discovery=discovery_data, test=test_data)
correlation_list <- list(discovery=discovery_correlation, test=test_correlation)
network_list <- list(discovery=discovery_network, test=test_network)

## -------------------------------------------------------------------------------------------------
# Each row corresponds to a module
preservation$observed
preservation$p.values

## -------------------------------------------------------------------------------------------------
max_pval <- apply(preservation$p.value, 1, max)
max_pval

## ---- dev="png", fig.height=7, fig.width=7, fig.align="center"------------------------------------
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,2,3,4),
  discovery="discovery", test="discovery"
)

## ---- dev="png", fig.height=7, fig.width=7, fig.align="center"------------------------------------
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,2,3,4),
  discovery="discovery", test="test"
)

## ---- dev="png", fig.height=7, fig.width=7, fig.align="center"------------------------------------
plotModule(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=c(1,2,3,4),
  discovery="discovery", test="test", orderNodesBy="test", orderModules=FALSE
)

## ---- dev="png", fig.height=7, fig.width=7, fig.align="center"------------------------------------
par(mar=c(6,6,6,4)+0.1)
plotCorrelation(
  data=data_list, correlation=correlation_list, network=network_list, 
  moduleAssignments=module_labels, modules=0:4, discovery="discovery",
  test="discovery", symmetric=TRUE, orderModules=FALSE
)

## -------------------------------------------------------------------------------------------------
# The summary profile of module 1 in the discovery dataset:
properties[["discovery"]][["1"]][["summary"]]
# Along with the proportion of variance in the module data explained by the 
# summary profile:
properties[["discovery"]][["1"]][["coherence"]]

# The same information in the test dataset:
properties[["test"]][["1"]][["summary"]]
properties[["test"]][["1"]][["coherence"]]

