## ------------------------------------------------------------------------
library("NetRep")
data("NetRep")
ls()

## ------------------------------------------------------------------------
# Convert the data to the 'bigMatrix' format:
discovery_data <- as.bigMatrix(discovery_data)
discovery_correlation <- as.bigMatrix(discovery_correlation)
discovery_network <- as.bigMatrix(discovery_network)
test_data <- as.bigMatrix(test_data)
test_correlation <- as.bigMatrix(test_correlation)
test_network <- as.bigMatrix(test_network)

## ------------------------------------------------------------------------
data_list <- list(discovery=discovery_data, test=test_data)
correlation_list <- list(discovery=discovery_correlation, test=test_correlation)
network_list <- list(discovery=discovery_network, test=test_network)

## ------------------------------------------------------------------------
# Each row corresponds to a module
preservation$observed
preservation$p.values

