# Are my networks preserved and reproducible?
# 
# @param datSets A list of \code{\link{big.matrix}} the underlying data the 
#  networks are calculated from.
# @param adjSets A list of \code{\link{big.matrix}}
# @param discoverySets A numeric vector
# @param testSets A numeric vector
# @param nodeLabelSets A list of vectors, one for each discovery data set, 
#  assigning each node in the discovery dataset to a subnetwork.
# @param ignoreSets An optional list of vectors, one for each discovery data 
#  set, specifying which subnetworks to ignore. 
# @param nPerm number of permutations to use.
# @param verbose logical indicating whether to print detailed output about the
#  function's execution.
# @param simplify logical. If \code{FALSE} the structure of the returned list
#  will be regular: that is, a list of lists, where the top level indicates a
#  summary network.

#' TODO
#' 
#' @export
netRep <- function() {
  NULL
}

# netRep <- function(
#     datSets, adjSets, discoverySets, testSets, nodeLabelSets, ignoreSets,
#     
#     nPerm = 10000, verbose=TRUE, simplify=TRUE
#   ) {
# 
#   
#   # Sanity Check Input, and set up data for the five different cases we have to
#   # deal with:
#   #  1. Underlying Data and Networks provided for all.
#   #  2. Only the networks are provided, with no underlying data.
#   #  3. Only the discovery network provided, and underlying data 
#   #     for all others.
#   #  4. Only the discovery network provided, and underlying data provided 
#   #     for all.
#   #  5. Only the underlying data is provided.
#   #
#   # TODO: provide user-friendly error messages
#   if (!missing(datSets) & !missing(adjSets) & 
#         all(!is.null(datSets)) & all(!is.null(adjSets))) {
#     case = 1 # case is used for simplifying output 
#     checkadjSets(adjSets)
#     checkDatSets(datSets)
#     
#     # Check concordance between networks and underlying data
#     stopifnot(length(datSets) == length(adjSets))
#     stopifnot(sort(names(datSets)) == sort(names(adjSets)))
# 
#     nNets <- length(adjSets)
#   } else if (!missing(adjSets) & missing(datSets)) {
#     case = 2
#     nNets <- length(adjSets)
#   }
#   
# 
#   if (simplify) {
#     # TODO
#     return(res)
#   } else {
#     return(res)
#   }
#   
# }

#' Assessing Network Replication
#' 
#' This function provides the main functionality, but requires well formed 
#' input. See details for instructions on usage. A wrapper function has been 
#' provided for usability purposes, see \code{\link{netRep}} for details.
#' 
#' @details
#'  Any function argument ending in \code{"Sets"} expects a list, where each 
#'  element corresponds to a dataset. \code{adjSets} should contain a list of 
#'  adjacency matrices, one for each dataset. \code{datSets} should contain a 
#'  list of \code{big.matrix} objects, which hold the underlying data the 
#'  corresponding lement of the \code{adjSets} list (for example, a matrix of 
#'  gene expression data). Here, we expect the rows of each element of 
#'  \code{datSets} to correspond to a node in the network.
#'
#'  Providing \code{datSets} is optional, if not provided, only statistics 
#'  computable on the adjacency matrices will be calculated. If the underlying
#'  data is provided for any dataset, the corresponding element \code{adjSet}
#'  can be left as \code{NULL}, and the adjacency matrix will be constructed 
#'  using the provided \code{buildNetFun} on the corresponding element of 
#'  \code{datSet}.
#'
#'  The elements of \code{nodeLabelSets} should be \code{NULL}, if its 
#'  corresponding dataset is not a \emph{"discovery"} dataset, or contain a 
#'  vector assigning each node in the network to a cluster/module/subnetwork 
#'  etc. \code{netRep} will test each of these network subsets in the 
#'  
#'  The order of the datasets should match across all \code{"Sets"}.
#'  
#'  If \code{datSets} is \code{NULL}, then this function assumes a full list of
#'  \code{adjSets} has been provided for each dataset of interest. 
#'  
#'  \code{ignoreSets} and \code{includeSets} should not overlap.
#' 
#' @references
#'   \enumerate{
#'     \item{
#'       Langfelder, P., Luo, R., Oldham, M. C. & Horvath, S. \emph{Is my 
#'       network module preserved and reproducible?} PLoS Comput. Biol. 
#'       \strong{7}, e1001057 (2011). 
#'     }
#'     \item{
#'       Knijnenburg, T. A., Wessels, L. F. A., Reinders, M. J. T. & Shmulevich, 
#'       I. \emph{Fewer permutations, more accurate P-values}. Bioinformatics 
#'       \strong{25}, i161–8 (2009). 
#'     }
#'   }
#' 
#' 
#' @param datSets \code{NULL}, or a list of \code{\link{big.matrix}}: one for 
#'   each dataset. 
#' @param adjSets \code{NULL}, or a list of \code{\link{big.matrix}}: one for 
#'   each dataset, corresponding to the adjacency matrix of edge weights between
#'   each pair of nodes in the network. Alternatively, if the list elements 
#'   are \code{NULL}, then the networks will be dynamically calculated using the 
#'   provided \code{buildNetFun} on the corresponding \code{datSet}.
#' @param nodeLabelSets a list, whose elements are \code{NULL} for each 
#'   \emph{test} dataset, and a vector for each \emph{discovery} dataset 
#'   assigning each node to a sub-network/cluster/module/component.
#' @param discovery a numeric vector indicating which elements of
#'   \code{datSets} and/or \code{adjSets} are to be treated as the 
#'   \emph{discovery} datasets.
#' @param test a numeric vector indicating which elements of \code{datSets}
#'   and/or \code{adjSets} are to be treated as the \emph{test} datasets.
#' @param nPerm number of permutations to use.
#' @param buildNetFun A function for constructing a network from the 
#'   \code{datSets}. This is used only where a pre-constructed network is not 
#'   provided in \code{adjSets}.
#' @param ignoreSets An optional list, where the elments for each 
#'   \emph{discovery} data set are vectors specifying which sub-networks to 
#'   skip.
#' @param includeSets An optional list, where the elments for each 
#'   \emph{discovery} data set are vectors specifying which sub-networks to 
#'   include.
#' @param ignoreDiag logical; ignore the diagonal of the adjacency matrix when
#'   assessing network subset preservation. Defaults to \code{TRUE}.
#' @param null the type of null model, either "overlap" or "all". If "overlap"
#'   (default) only the nodes present in both the discovery and test networks
#'   will be used to draw the null distribution for each statistic. If "all", 
#'   all nodes in the test network will be used to draw the null distribution.
#' @param tailApprox logical; if \code{TRUE}, the tail approximation method from
#'   (\emph{2}) is used to obtain an estimated p-value for network statistics
#'   which are more extreme than the null distribution obtained through the 
#'   permutation procedure.
#' @param verbose logical; print output while the function is running 
#'   (\code{TRUE} by default).
#' @param indent numeric; a positive value indicating the indent level to start
#'   the output at. Defaults to 0. Each indent level adds two spaces to the 
#'   start of each line of output.
#'
#' @import foreach
#' @importFrom itertools isplitIndices
netRep.core <- function(
    datSets=NULL, adjSets=NULL, nodeLabelSets, discovery, test, nPerm=10000,
    buildNetFun, ignoreSets=NULL, includeSets=NULL, ignoreDiag=TRUE, 
    null="overlap", tailApprox=TRUE, verbose=TRUE, indent=0
  ) {
  # The following declarations are for iterators declared inside each foreach 
  # loop. Declarations are required to satisfy NOTES generated by R CMD check, 
  # and also serve as useful documentation for the reader of the source code.
  di <- NULL # discovery dataset index in the "Sets" lists.
  ti <- NULL # test dataset index in the "Sets" lists.
  ss <- NULL # an individual network subset label or index.
  chunk <- NULL # chunk of permutations to run on a core.
  i <- NULL # current permutation number.
  stat <- NULL # column index for each statistic computed in the null 3D array.
  
  # Identify the null model to use 
  nullModels <- c("overlap", "all")
  if (is.na(pmatch(null, nullModels))) {
    stop("overlap must match one of ", paste(nullModels, collapse=" "))
  }
  model <- pmatch(null, nullModels)
  
  nCores <- getDoParWorkers()
  vCat(verbose, indent, "Running on", nCores, "cores.")
  
  nNets <- length(nodeLabelSets)
  
  # Set up names for more informative output, if the datasets and their 
  # corresponding networks have been given names.
  if (is.null(names(nodeLabelSets))) {
    setNames <- 1:nNets
  } else {
    setNames <- names(nodeLabelSets)
  }
  
  # Iterate pairwise over data-sets, comparing those marked "discovery"
  # with each marked as "replication".
  foreach(di=1:nNets) %:% foreach(ti=1:nNets) %do% {
    if ((di %in% discovery) & (ti %in% test) & (di != ti)) {
      vCat(verbose, indent, "Calculating preservation of network subsets from",
           "dataset", setNames[di], ", in dataset ", setNames[ti])
      
      # Get a vector of nodes which are present in both datasets. Depends on 
      # the combination of data input provided.
      if (is.null(adjSets[[di]])) {
        if(is.null(adjSets[[ti]])) {
          oNodes <- rownames(datSets[[di]]) %sub_in% rownames(datSets[[ti]])
        } else {
          oNodes <- rownames(datSets[[di]]) %sub_in% rownames(adjSets[[ti]])
        }
      } else {
        if(is.null(adjSets[[ti]])) {
          oNodes <- rownames(adjSets[[di]]) %sub_in% rownames(datSets[[ti]])
        } else {
          oNodes <- rownames(adjSets[[di]]) %sub_in% rownames(adjSets[[ti]])
        }
      }
      if (length(oNodes) == 0) {
        warning("No nodes in dataset", setNames[di],  "are present in dataset ", 
                setNames[ti], ", skipping.")
        return(NULL)
      }
      
      # TODO: 
      # Set up empty big.matices where we are building the networks on the fly
      if (is.null(adjSets[[di]])) {
        stop("not implemented yet")
      }
      if (is.null(adjSets[[ti]])) {
        stop("not implemented yet")
      }
      
      # Set the diagonals to NA if we're ignoring them in our calculations
      if (ignoreDiag) {
        oldDiags <- list(diag(adjSets[[di]]), diag(adjSets[[ti]]))
        diag(adjSets[[di]]) <- NA
        diag(adjSets[[ti]]) <- NA
        on.exit({
          diag(adjSets[[di]]) <- oldDiags[[1]]
          diag(adjSets[[ti]]) <- oldDiags[[2]]
        })
      }
      
      # Compute information about the network subsets, their size, and what 
      # proportion of each is overlapping the test network.
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
      # Get the size of each of the subsets of interest in the discovery dataset
      dSizes <- table(nodeLabelSets[[di]])
      dSizes <- dSizes[names(dSizes) %in% dSubsets] 
      # Get the size of the overlap of each of the subsets of interest between
      # the discovery and test datasets.
      oSizes <- table(nodeLabelSets[[di]][oNodes])
      oSizes <- oSizes[names(oSizes) %in% oSubsets]
      # Need to handle the case where a subset of interest has no overlap.
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
            discNet = match(ons, rownames(adjSets[[di]])),
            testData = match(ons, rownames(datSets[[ti]])),
            testNet = match(ons, rownames(adjSets[[ti]]))
          )
        datTypeDict <- lapply(datTypeDict, `names<-`, ons)
      })
      names(subsetDict) <- oSubsets
      
      # Calculate observed network statistics across datasets for each 
      # sub-network
      vCat(verbose, indent+1, "Calculating observed test statistics...")
      
      observed <- foreach(ss=oSubsets, .combine=rbind) %do% {
        calcReplStats(datSets[[di]], adjSets[[di]], datSets[[ti]], adjSets[[ti]], 
                      subsetDict[[ss]])
      }
      rownames(observed) <- oSubsets
      vCat(verbose, indent+1, "Done!")
      
      # Calculate the null distribution for each of the statistics.
      vCat(verbose, indent+1, "Calculating null distributions with", nPerm, 
           "permutations...")
      # To log progress, we will write our progress to a file for each chunk,
      # which can be monitored from the terminal, or another R session 
      dir.create("run-progress") 
      nulls <- foreach(chunk=isplitIndices(nPerm, chunks=nCores), .combine=abind3) %dopar% {
        pb <- setupParProgressLogs(chunk)
        on.exit(close(pb))
        foreach(i=chunk, .combine=abind3) %:% foreach(ss=oSubsets, .combine=rbind) %do% {
          on.exit(updateParProgress(pb, i))
          # Randomise module assignments for each run.
          if (model == "overlap") {
            permuted <- sample(oNodes, size=oSizes[ss])
          } else {
            permuted <- sample(tNodes, size=oSizes[ss])
          }
          
          subsetDict[[ss]][c("testData", "testNet")] <- list(
              match(permuted, rownames(datSets[[ti]])),
              match(permuted, rownames(adjSets[[ti]]))
            )
          calcReplStats(datSets[[di]], adjSets[[di]], datSets[[ti]], 
                        adjSets[[ti]], subsetDict[[ss]])
        }
      }
      dimnames(nulls)[[1]] <- oSubsets
      dimnames(nulls)[[2]] <- colnames(observed)
      dimnames(nulls)[[3]] <- paste("permutation", seq_len(nPerm), sep=".")
      unlink("run-progress", recursive=TRUE) # Remove progress logs
      
      # Calculate the p-value for the observed statistic based on the null 
      # distribution
      vCat(verbose, indent+1, "Calculating P-values...")
      p.values <- foreach(stat=seq_len(ncol(observed)), .combine=cbind) %:% 
        foreach(ss=seq_len(nrow(observed)), .combine=c) %do% {
          pperm(permuted[ss,stat,], observed[ss,stat], tailApprox,
                lower.tail=FALSE)
      }
      dimnames(p.values) <- dimnames(observed)
      vCat(verbose, indent+1, "Done!")
      
      # Collate results
      return(list(observed=observed, null=nulls, p.value=p.values,
                  overlap=overlap, overlapSize=oSizes))
    } else {
      # We are not currently comparing these two datasets.
      return(NULL)
    }
  }
}
