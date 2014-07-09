


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
#' @param nPerm number of permutations to use.
#' @param buildNetFun A function for constructing a network from the 
#'   \code{datSets}. This is used only where a pre-constructed network is not 
#'   provided in \code{netSets}.
#' @param ignoreSets An optional list, where the elments for each 
#'   \emph{discovery} data set are vectors specifying which sub-networks to 
#'   skip.
#' @param includeSets An optional list, where the elments for each 
#'   \emph{discovery} data set are vectors specifying which sub-networks to 
#'   include.
#' @param verbose logical; print output while the function is running 
#'   (\code{TRUE} by default).
#' @param ident numeric; a positive value indicating the indent level to start
#'   the output at. Defaults to 0. Each indent level adds two spaces to the 
#'   start of each line of output.
#'
#' @import foreach
#' @importFrom itertools isplitIndices
netRep.core <- function(
    datSets=NULL, netSets=NULL, nodeLabelSets, discovery, test, nPerm=10000,
    buildNetFun, ignoreSets=NULL, includeSets=NULL,
    verbose=TRUE, indent=0
  ) {
  
  nCores <- getDoParWorkers()
  vCat(verbose, "Running on", nCores, "cores.", indent)
  
  nNets <- length(nodeLabelSets)
  
  # Iterate pairwise over data-sets, comparing those marked "discovery"
  # with each marked as "replication".
  foreach(di=1:nNets) %:% foreach(ti=1:nNets) %do% {
    if ((di %in% discovery) & (ti %in% test) & (di != ti)) {
      # di: index of the discovery network
      # ti: index of the test network
      vCat(verbose, "Calculating preservation of network subsets from dataset", 
           discNames[di], ", in dataset ", testNames[ti], indent)
      
      # Get a vector of nodes which are present in both datasets. Depends on 
      # the combination of data input provided.
      if (is.null(netSets[[di]])) {
        if(is.null(netSets[[ti]])) {
          oNodes <- rownames(datSets[[di]]) %sub_in% rownames(datSets[[ti]])
          tnodes
        } else {
          oNodes <- rownames(datSets[[di]]) %sub_in% rownames(netSets[[ti]])
        }
      } else {
        if(is.null(netSets[[ti]])) {
          oNodes <- rownames(netSets[[di]]) %sub_in% rownames(datSets[[ti]])
        } else {
          oNodes <- rownames(netSets[[di]]) %sub_in% rownames(netSets[[ti]])
        }
      }
      if (length(oNodes) == 0) {
        warning("No nodes in dataset", discNames[di], 
                "are present in dataset ", testNames[ti], ", skipping.")
        return(NULL)
      }
      
      # TODO: Handle cases where we are building the networks on the fly
      if (is.null(netSets[[di]])) {
        stop("not implemented yet")
      }
      if (is.null(netSets[[ti]])) {
        stop("not implemented yet")
      }
      
      # Compute information about the network subsets, their size, and what 
      # proportion of each is overlapping the test network.
      dSizes <- table(nodeLabelSets[[di]])
      dSubsets <- unique(nodeLabelSets[[di]])
      oSubsets <- unique(nodeLabelSets[[di]][oNodes])
      # Only look at subsets of interest, if specified:
      if (!is.null(ignoreSets[[di]])) {
        dSubsets <- dSubsets %sub_nin% ignoreSets[[di]]
        oSubsets <- oSubsets %sub_nin% ignoreSets[[di]]
      }
      if (!is.null(includeSets[[di]])) {
        dSubsets <- dSubsets %sub_in% includeSets[[di]]
        oSubsets <- oSubsets %sub_in% includeSets[[di]]
      }
      # most of this code is to deal with subsets who have no nodes in the test
      # network.
      oSizes <- table(nodeLabelSets[[di]][oNodes])
      oSizes <- c(oSizes, rep(0, length(dSubsets %sub_nin% oSubsets)))
      names(oSizes) <- c(names(oSizes) %sub_nin% "", dSubsets %sub_nin% oSubsets)
      overlap <- oSizes[names(dSizes)]/dSizes
      
      # Create dictionary of network subsets labels to individual nodes. While
      # we only care about nodes in both networks, we need to know where those
      # nodes are located in both networks/datasets.
      subsetDict <- lapply(oSubsets, function(ss) {
        # get all overlapping nodes in this subset
        ons <- names(which(nodeLabelSets[[di]][oNodes] == ss))
        # Now get a mapping for these nodes to their respective indices in each
        # dataset/network
        datTypeDict <- list(
            discData = match(ons, rownames(datSets[[di]])),
            discNet = match(ons, rownames(netSets[[di]])),
            testData = match(ons, rownames(datSets[[ti]])),
            testNet = match(ons, rownames(netSets[[ti]]))
          )
        datTypeDict <- lapply(datTypeDict, `names<-`, ons)
      })
      names(subsetDict) <- oSubsets
      
      # Calculate observed network statistics across datasets for each 
      # sub-network
      vCat(verbose, indent+1, "Calculating observed test statistics...")
      
      observed <- foreach(ss=subsets, .combine=rbind) %do% {
        calcReplStats(datSets[[di]], netSets[[di]], datSets[[ti]], netSets[[ti]], 
                      subsetDict[[ss]])
      }
      rownames(observed) <- subsets
      vCat(verbose,  indent+1, "Done!")
      
      vCat(verbose,  indent+1, "Calculating null distributions with", 
           nPerm, "permutations...")
      # Calculate the null distribution for each of the statistics
      nulls <- foreach(chunk=isplitIndices(nPerm, chunks=nCores), .combine=abind3) %dopar% {
        foreach(i=chunk, .combine=abind3) %:% foreach(ss=subsets, .combine=rbind) %do% {
          # Randomise module assignments for each run.
          permuted <- sample(oNodes, size=oSizes[ss])
          subsetDict[[ss]][c("testData", "testNet")] <- list(
              match(permuted, rownames(datSets[[ti]])),
              match(permuted, rownames(netSets[[ti]]))
            )
          calcReplStats(datSets[[di]], netSets[[di]], datSets[[ti]], 
                        netSets[[ti]], subsetDict[[ss]])
        }
      }
      dimnames(nulls)[[1]] <- subsets
      dimnames(nulls)[[2]] <- colnames(observed)
      dimnames(nulls)[[3]] <- paste("permutation", seq_len(nPerm), sep=".")
      
      # Calculate the p-value for the observed statistic based on the null 
      # distribution
      vCat(verbose, indent+1, "Calculating P-values...")
      p.values <- foreach(j=seq_len(ncol(observed)), .combine=cbind) %:% 
        foreach(i=seq_len(nrow(observed)), .combine=c) %do% {
          pperm(permuted[i,j,], observed[i,j], lower.tail=FALSE)
      }
      dimnames(p.values) <- dimnames(observed)
      vCat(verbose, ident+1, "Done!")
      
      # Collate results
      return(list(observed=observed, null=null, p.value=p.values,
                  overlap=overlap, overlapSize=oSizes))
    } else {
      # We are not currently comparing these two datasets.
      return(NULL)
    }
  }
}
