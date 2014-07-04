
#' Are my networks preserved and reproducible?
#' 
#' @param datSets A list of \code{\link{big.matrix}} the underlying data the networks are 
#'        calculated from.
#' @param netSets A list of \code{\link{big.matrix}}
#' @param nodeLabelSets A list of vectors, one for each discovery data set, 
#'        assigning each node in the discovery network to a cluster
#' @param discoverySets A numeric vector
#' @param testSets A numeric vector
#' @param nPerm number of permutations to use.
#' @param verbose logical indicating whether to print detailed output about the
#'        function's execution
#' @param simplify logical. If \code{FALSE} the structure of the returned list
#'        will be regular: that is, a list of lists, where the top level 
#'        indicates a summary network
#"
networkReplication <- function(
    datSets, netSets, 
    discoverySets, testSets,
    nPerm = 10000, verbose=TRUE, simplify=TRUE
  ) {
  require(abind)
  require(foreach)
  require(itertools)
  
  # Sanity Check Input, and set up data for the five different cases we have to
  # deal with:
  #  1. Underlying Data and Networks provided for all.
  #  2. Only the networks are provided, with no underlying data.
  #  3. Only the discovery network provided, and underlying data 
  #     for all others.
  #  4. Only the discovery network provided, and underlying data provided 
  #     for all.
  #  5. Only the underlying data is provided.
  #
  # TODO: provide user-friendly error messages
  if (!missing(datSets) & !missing(netSets) & 
        all(!is.null(datSets)) & all(!is.null(netSets))) {
    case = 1 # case is used for simplifying output
    # Make sure our data is the right class
    stopifnot(all(sapply(datSets, class)) == "big.matrix")
    stopifnot(all(sapply(netSets, class)) == "big.matrix")
    # And that order matches between networks and underlying data
    stopifnot(length(datSets) == length(netSets))
    stopifnot(names(datSets) == names(netSets))
     
    nNets <- length(netSets)
  } else if (!missing(netSets) & missing(datSets)) {
    case = 2
    nNets <- length(netSets)
  }
  

  
  nCores <- getDoParWorkers()
  vCat(verbose, "Running on", nCores, "cores.")
  
  # Iterate pairwise over data-sets, comparing those marked "discovery"
  # with each marked as "replication".
  res <- foreach(disc=1:nNets) %:% foreach(test=1:nNets) %do% {
    if ((disc %in% discIndices) & (test %in% testIndices) & (test != disc)) {
      vCat(verbose, "Calculating preservation of networks from dataset", 
           discNames[disc], ", in dataset ", testNames[test])
      # Get a list of the discovery network nodes that are also in the test
      # network
      if (case == 1) {
        discNodes <- rownames(netSets[[test]])
      }
      
    }
  }
  if (simplify) {
    # TODO
    return(res)
  } else {
    return(res)
  }
  
}