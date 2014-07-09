


#' Are my networks preserved and reproducible?
#' 
#' @param datSets A list of \code{\link{big.matrix}} the underlying data the 
#'  networks are calculated from.
#' @param netSets A list of \code{\link{big.matrix}}
#' @param discoverySets A numeric vector
#' @param testSets A numeric vector
#' @param nodeLabelSets A list of vectors, one for each discovery data set, 
#'  assigning each node in the discovery dataset to a subnetwork.
#' @param ignoreSets An optional list of vectors, one for each discovery data 
#'  set, specifying which subnetworks to ignore. 
#' @param nPerm number of permutations to use.
#' @param verbose logical indicating whether to print detailed output about the
#'  function's execution.
#' @param simplify logical. If \code{FALSE} the structure of the returned list
#'  will be regular: that is, a list of lists, where the top level indicates a
#'  summary network.
#"
netRep <- function(
    datSets, netSets, discoverySets, testSets, nodeLabelSets, ignoreSets,
    
    nPerm = 10000, verbose=TRUE, simplify=TRUE
  ) {

  
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
    checkNetSets(netSets)
    checkDatSets(datSets)
    
    # Check concordance between networks and underlying data
    stopifnot(length(datSets) == length(netSets))
    stopifnot(sort(names(datSets)) == sort(names(netSets)))

    nNets <- length(netSets)
  } else if (!missing(netSets) & missing(datSets)) {
    case = 2
    nNets <- length(netSets)
  }
  

  if (simplify) {
    # TODO
    return(res)
  } else {
    return(res)
  }
  
}

#' Core Function for Assessing Network Replication
#' 
#' This function provides the main functionality, but requires well formed 
#' input. A wrapper function has been provided for usability purposes, see 
#' \code{\link{networkReplication}} for details.
#' 
#' @details
#'  Any function argument ending in "Sets" expects a list, where each element
#'  corresponds to a dataset. The order of the datasets should match across all
#'  arguments.  
#'  
#'  If \code{datSets} is \code{NULL}, then this function assumes a full list of
#'  \code{netSets} has been provided for each dataset of interest. 
#'  
#'  \code{ignoreSets} and \code{includeSets} should not overlap.
#' 
#' @param datSets \code{NULL}, or a list of \code{\link{big.matrix}}: one for 
#'   each dataset. 
#' @param netSets \code{NULL}, or a list of \code{\link{big.matrix}}: one for 
#'   each dataset. Alternatively, if the list elements are \code{NULL}, then the
#'   networks will be dynamically calculated using the provided 
#'   \code{buildNetFun} on the corresponding \code{datSet}.
#' @param nodeLabelSets a list, whose elements are \code{NULL} for each 
#'   \emph{test} dataset, and a vector for each \emph{discovery} dataset 
#'   assigning each node to a sub-network/cluster/module/component.
#' @param discovery a numeric vector indicating which elements of
#'   \code{datSets} and/or \code{netSets} are to be treated as the 
#'   \emph{discovery} datasets.
#' @param test a numeric vector indicating which elements of \code{datSets}
#'   and/or \code{netSets} are to be treated as the \emph{test} datasets.
#' @param buildNetFun A function for constructing a network from the 
#'   \code{datSets}. This is used only where a pre-constructed network is not 
#'   provided in \code{netSets}.
#' @param ignoreSets An optional list, where the elments for each 
#'   \emph{discovery} data set are vectors specifying which sub-networks to 
#'   skip.
#' @param includeSets An optional list, where the elments for each 
#'   \emph{discovery} data set are vectors specifying which sub-networks to 
#'   include. 
#'
#' @import foreach
#' @importFrom itertools isplitIndices
netRep.core <- function(
    datSets=NULL, netSets=NULL, nodeLabelSets, discovery, test, 
    buildNetFun, ignoreSets, includeSets) {

  
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
        discNodes <- rownames(netSets[[disc]]) %sub_in% rownames(netSets[[test]])
      }
      
      if (length(discNodes)) {
        warning("No nodes in dataset ", discNames[disc], 
                " are present in dataset ", testNames[test], ", skipping.")
        return(NULL)
      }
      
      # Get information about the sub-networks.
      sizes <- table(nodeLabelSets[[disc]][discNodes])
      overlap <- sizes/table(nodeLabelSets[[disc]])
      subnets <- unique(nodeLabelSets[[disc]][discNodes])
      
      # Create dictionary of sub-nets to -> nodes.
      moduleIndices <- sapply(subnets, function(name) {
        which(nodeLabelSets[[disc]][discNodes] == name)
      })
      
      # Calculate observed network statistics across datasets for each 
      # sub-network
      vCat(verbose, "Calculating observed statistics...")
      observed <- foreach(net=subnets, .combine=rbind) %do% {
        calculatePreservation(datSets[[disc]], netSets[[disc]],
                              datSets[[disc]], netSets[[test]],
                              moduleIndices[[module]])
      }
      rownames(observed) <- moduleNames
      
    }
  }
}

# Helper function for checking input sanity
checkNetSets <- function(netSets) {
  stopifnot(sapply(netSets, class) %in% c("big.matrix", "NULL"))
  stopifnot(sapply(netSets[!is.null(netSets)], function(x) !is.null(rownames(x))))
  stopifnot(sapply(netSets[!is.null(netSets)], function(x) !is.null(colnames(x))))
}

# Helper function for checking input sanity
checkDatSets <- function(datSets) {
  stopifnot(sapply(datSets, class) %in% c("big.matrix", "NULL"))
  stopifnot(sapply(datSets[!is.null(datSets)], function(x) !is.null(rownames(x))))
}