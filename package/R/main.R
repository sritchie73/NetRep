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
#'  If \code{verbose} is \code{TRUE} and a parallel backend has been registered,
#'  then one of the cores will be reserved for reporting the function's progress.
#'  You may register more "cores" than are physically have on the machine, since
#'  this master core has relatively low overhead. Note that if only 2 cores have
#'  been registered, this is the same as running sequentially, since
#'  
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
#'       \strong{25}, i161-8 (2009). 
#'     }
#'   }
#' 
#' 
#' @param datSets \code{NULL}, or a list of \code{\link{big.matrix}}: one for 
#'   each dataset. 
#' @param varNameSets \code{NULL}, or a list of character vectors giving the 
#'   \code{rownames} of \code{datSets}.
#' @param adjSets \code{NULL}, or a list of \code{\link{big.matrix}}: one for 
#'   each dataset, corresponding to the adjacency matrix of edge weights between
#'   each pair of nodes in the network. Alternatively, if the list elements 
#'   are \code{NULL}, then the networks will be dynamically calculated using the 
#'   provided \code{buildNetFun} on the corresponding \code{datSet}.
#' @param nodeNameSets \code{NULL}, or a list of character vectors giving the 
#'   \code{rownames} of \code{adjSets}.
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
#' @param simplify logical; if \code{TRUE}, simplify the structure of the output
#'  list if possible.
#'
#' @import foreach
#' @importFrom itertools isplitIndices
#' @export
netRepMain <- function(
  datSets=NULL, varNameSets=NULL, adjSets=NULL, nodeNameSets=NULL, 
  nodeLabelSets, discovery, test, nPerm=10000,
  buildNetFun, ignoreSets=NULL, includeSets=NULL, 
  null="overlap", tailApprox=FALSE, verbose=TRUE, indent=0, 
  simplify=TRUE
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
  
  # Determine set up of worker nodes (if applicable)
  nCores <- getDoParWorkers()
  if (verbose & nCores > 1) {
    vCat(verbose, indent, "Verbose output is TRUE; reserving one core as the", 
         "master, to report the progress of the worker cores who do the bulk",
         "of the calculation.")
    nWorkers <- nCores - 1
  } else {
    nWorkers <- nCores
  }
  vCat(verbose, indent, "Running with", nWorkers, "worker cores.")
  
  nNets <- length(nodeLabelSets)
  
  # Set up names for more informative output, if the datasets and their 
  # corresponding networks have been given names.
  if (is.null(names(nodeLabelSets))) {
    setNames <- 1:nNets
  } else {
    setNames <- names(nodeLabelSets)
  }
  
  # Set up temporary directory for working big.matrix objects
  scaledSets <- NULL
  if (!is.null(datSets)) {
    dir.create(".temp-objects", showWarnings=FALSE)
    on.exit({
      unlink(".temp-objects", recursive=TRUE)
    }, add=TRUE)
    
    scaledSets <- rep(list(NULL), length(datSets))
  }
  
  # Set up return list 
  res <- rep(list(NULL), nNets)
  res <- lapply(res, function(x) { 
    l <- rep(list(NULL), nNets); 
    names(l) <- names(nodeLabelSets) 
  })
  names(res) <- names(nodeLabelSets)
  
  # Iterate pairwise over datasets, comparing those marked "discovery"
  # with each marked as "test".
  for (di in seq_len(nNets)) {
    for (ti in seq_len(nNets)) {
      tryCatch({
        if ((di %in% discovery) & (ti %in% test) & (di != ti)) {
        # Set up return list
        res[[di]][[ti]] <- rep(list(NULL), 5)
        
        # Output messages
        vCat(verbose, indent, sep="", 
             "Calculating preservation of network subsets from dataset ",
             setNames[di], ", in dataset ", setNames[ti], ".")

        # Attach relevant matrices
        vCat(verbose, indent+1, "Attaching and checking matrices...")
        if (is.null(adjSets[[di]])) {
          stop("not implemented yet")
        } else {
          discAdj <- attach.big.matrix(adjSets[[di]])
          checkFinite(discAdj)
        }
        if (is.null(adjSets[[ti]])) {
          stop("not implemented yet")
        } else {
          testAdj <- attach.big.matrix(adjSets[[ti]])
          checkFinite(testAdj)
        }
        if (!is.null(datSets[[di]])) {
          discDat <- attach.big.matrix(datSets[[di]])
          checkFinite(discDat)
        } else {
          discDat <- NULL
        }
        if (!is.null(datSets[[ti]])) {
          testDat <- attach.big.matrix(datSets[[ti]])
          checkFinite(testDat)
        } else {
          testDat <- NULL
        }
        
        # Create scaled data 
        if (!is.null(discDat)) {
          vCat(verbose, indent+1, "Creating temporary scaled datasets...")
          if (is.null(scaledSets[[di]])) {
            descriptor <- paste0("scaled", di, ".desc")
            backing <- paste0("scaled", di, ".bin")
            scaledDisc <- scaleBigMatrix(
              discDat, backing, ".temp-objects", descriptor
            )
            scaledSets[[di]] <- file.path(".temp-objects", descriptor)
          } else {
            scaledDisc <- attach.big.matrix(scaledSets[[di]])
          }
        }
        if (!is.null(testDat)) {
          if (is.null(scaledSets[[ti]])) {
            descriptor <- paste0("scaled", ti, ".desc")
            backing <- paste0("scaled", ti, ".bin")
            scaledTest <- scaleBigMatrix(
              testDat, backing, ".temp-objects", descriptor
            )
            scaledSets[[ti]] <- file.path(".temp-objects", descriptor)
          } else {
            scaledTest <- attach.big.matrix(scaledSets[[ti]])
          }
        }
        
        # Get a vector of nodes which are present in both datasets. Depends on 
        # the combination of data input provided.
        vCat(verbose, indent+1, "Extracting information about node overlap...")
        if (is.null(testAdj)) {
          if(is.null(discAdj)) {
            oNodes <- varNameSets[[di]] %sub_in% varNameSets[[ti]]
          } else {
            oNodes <- nodeNameSets[[di]] %sub_in% varNameSets[[ti]]
          }
          tNodes <- varNameSets[[ti]]
        } else {
          if(is.null(discAdj)) {
            oNodes <- varNameSets[[di]] %sub_in% nodeNameSets[[ti]]
          } else {
            oNodes <- nodeNameSets[[di]] %sub_in% nodeNameSets[[ti]]
          }
          tNodes <- nodeNameSets[[ti]]
        }
        # Force the node indices to be strings, even if the provided identifiers
        # are integers.
        oNodes <- as.character(oNodes)
        tNodes <- as.character(tNodes)
        
        if (length(oNodes) == 0) {
          warning("No nodes in dataset ", setNames[di],  
                  " are present in dataset ", setNames[ti], ", skipping.")
          return(NULL)
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
        
        vCat(verbose, indent+1, "Calculating observed test statistics...")
        poke(discAdj, discDat, scaledDisc, testAdj, testDat, scaledTest)
        # Obtain the topological properties for each network subset in the
        # discovery dataset, we only want to calculate these once!
        discProps <- foreach(ss=oSubsets) %do% {
          subsetNodes <- names(which(nodeLabelSets[[di]][oNodes] == ss))
          # get the indices in the underlying data and adjacency matrices for 
          # the subset nodes. Sorted, because sequential memory access is faster.
          discDatInd <- match(subsetNodes, varNameSets[[di]])
          discAdjInd <- match(subsetNodes, nodeNameSets[[di]])
          subsetProps(discAdj, discAdjInd, discDat, scaledDisc, discDatInd)
        }
        names(discProps) <- oSubsets
        
        # Now calculate the observed value for each network statistic
        observed <- foreach(ss=oSubsets, .combine=rbind) %do% {
          subsetNodes <- names(which(nodeLabelSets[[di]][oNodes] == ss))
          testDatInd <- match(subsetNodes, varNameSets[[ti]])
          testAdjInd <- match(subsetNodes, nodeNameSets[[ti]])
          discAdjInd <- match(subsetNodes, nodeNameSets[[di]])
          testProps <- subsetProps(
            testAdj, testAdjInd, testDat, scaledTest, testDatInd
          )
          return(c(
            calcSplitTestStats(discProps[[as.character(ss)]], testProps),
            calcSharedTestStats(discAdj, discAdjInd, testAdj, testAdjInd)
          ))
        }
        rownames(observed) <- oSubsets
        
        # Calculate the null distribution for each of the statistics.
        vCat(verbose, indent+1, "Calculating null distributions with", nPerm, 
             "permutations...")
        if(verbose) {
          # To log progress, we will write our progress to a file for each chunk
          dir.create("run-progress", showWarnings=FALSE)
        }
        nulls <- foreach(
          chunk=ichunkTasks(verbose, nPerm, nCores),
          .combine=abind3
        ) %maybe_do_par% {
          if (verbose & length(chunk) == 1) {
            if (chunk == -1) {
              monitorProgress(nWorkers, indent+2)
              NULL
            }
          } else {
            poke(discAdj, discDat, scaledDisc, testAdj, testDat, scaledTest)
            if (verbose) {
              conns <- setupParProgressLogs(chunk, nWorkers, indent+2)
              progressBar <- conns[[1]]
              on.exit(lapply(conns, close))
            } 
            foreach(i=chunk, .combine=abind3) %do% { 
              # Update the progress at the end of the loop.
              if (verbose) {
                on.exit({
                  updateParProgress(progressBar, i)
                  if (nCores == 1) reportProgress(indent+2)
                })                
              }
              foreach(ss=oSubsets, .combine=rbind) %do% {
                # Select a random subset of nodes of the same size as the subset 
                # ss, depending on our null model.
                if (model == "overlap") {
                  permNames <- sample(oNodes, size=oSizes[ss])
                } else {
                  permNames <- sample(tNodes, size=oSizes[ss])
                }
                permDatInd <- match(permNames, varNameSets[[ti]])
                permAdjInd <- match(permNames, nodeNameSets[[ti]])
                discAdjInd <- match(subsetNodes, nodeNameSets[[di]])
                
                testProps <- subsetProps(
                  testAdj, permAdjInd, testDat, scaledTest, permDatInd
                )
                return(c(
                  calcSplitTestStats(discProps[[as.character(ss)]], testProps),
                  calcSharedTestStats(discAdj, discAdjInd, testAdj, permAdjInd)
                ))
              }
            }
          }
        }
        dimnames(nulls)[[1]] <- oSubsets
        dimnames(nulls)[[2]] <- colnames(observed)
        dimnames(nulls)[[3]] <- paste("permutation", seq_len(nPerm), sep=".")
        
        # Calculate the p-value for the observed statistic based on the null 
        # distribution
        vCat(verbose, indent+1, "Calculating P-values...")
        p.values <- foreach(stat=seq_len(ncol(observed)), .combine=cbind) %:% 
          foreach(ss=seq_len(nrow(observed)), .combine=c) %do% {
            pperm(nulls[ss,stat,], observed[ss,stat], tailApprox, 
                  lower.tail=FALSE)
          }
        dimnames(p.values) <- dimnames(observed)
        
        # Collate results
        # First order output nicely
        tryCatch({
          # For modules that are integer coded, make sure they're numerically
          # ordered, not alphabetically.
          arrOrder <- order(as.integer(rownames(null)))
          vOrder <- order(as.integer(names(overlap)))
        }, warning = function(w) {
          # If we can't cast to an integer, sort normally.
          arrOrder <- order(rownames(null))
          vOrder <- order(overlap)
        })
        res[[di]][[ti]][[1]] <- nulls[arrOrder,,]
        res[[di]][[ti]][[2]] <- observed[arrOrder,]
        res[[di]][[ti]][[3]] <- p.values[arrOrder,]
        res[[di]][[ti]][[4]] <- overlap[vOrder]
        res[[di]][[ti]][[5]] <- oSizes[vOrder]
        names(res[[di]][[ti]]) <- c(
          "null", "observed", "p.value", "overlapProp", "overlapSize"
        )
        
        vCat(verbose, indent+1, "Cleaning up temporary objects...")
        unlink("run-progress", recursive=TRUE)
        rm(discDat, scaledDisc, discAdj, testDat, scaledTest, testAdj)
        gc()
        vCat(verbose, indent, "Done!")
      }
      }, error = function(e) {
        warning(e$message)
        vCat(TRUE, indent+1, "Pair Failed with error: ", e$message, 
              "\nSkipping")
      })
    }
  }
  if (simplify) {
    # remove doubly-nested list structure where possible
    if (length(discovery) == 1 && length(test) == 1) {
      res <- res[[discovery]][[test]]
    } else if (length(discovery) == 1) {
      res <- res[[discovery]][test] 
    } else if (length(test) == 1) {
      res <- lapply(res[discovery], `[[`, test)
    }
  }
  res
}