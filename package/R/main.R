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
#'  gene expression data). Here, we expect the columns of each element of 
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
#' @param corSets a list of \code{\link{big.matrix}} objects: one for each 
#'   dataset, corresponding to the pairwise-gene correlation. 
#' @param adjSets a list of \code{\link{big.matrix}} objects: one for 
#'   each dataset, corresponding to the adjacency matrix of edge weights between
#'   each pair of nodes in the network.
#' @param nodeNameSets \code{NULL}, or a list of character vectors giving the 
#'   node names for the corresponding element of \code{datSets}, \code{corSets},
#'   and \code{adjSets}.
#' @param nodeLabelSets a list, whose elements are \code{NULL} for each 
#'   \emph{test} dataset, and a vector for each \emph{discovery} dataset 
#'   assigning each node to a sub-network/cluster/module/component.
#' @param discovery a numeric vector indicating which elements of
#'   \code{datSets} and/or \code{adjSets} are to be treated as the 
#'   \emph{discovery} datasets.
#' @param test a numeric vector indicating which elements of \code{datSets}
#'   and/or \code{adjSets} are to be treated as the \emph{test} datasets.
#' @param nPerm number of permutations to use.
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
#' @param verbose logical; print output while the function is running 
#'   (\code{TRUE} by default).
#' @param indent numeric; a positive value indicating the indent level to start
#'   the output at. Defaults to 0. Each indent level adds two spaces to the 
#'   start of each line of output.
#' @param simplify logical; if \code{TRUE}, simplify the structure of the output
#'  list if possible.
#'
#' @import foreach
#' @import RhpcBLASctl
#' @export
netRepMain <- function(
  datSets=NULL, corSets, adjSets, nodeNameSets=NULL, 
  nodeLabelSets, discovery, test, nPerm=10000, ignoreSets=NULL, 
  includeSets=NULL, null="overlap", verbose=TRUE, indent=0, simplify=TRUE
) {
  # The following declarations are for iterators declared inside each foreach 
  # loop. Declarations are required to satisfy NOTES generated by R CMD check, 
  # and also serve as useful documentation for the reader of the source code.
  di <- NULL # discovery dataset index in the "Sets" lists.
  ti <- NULL # test dataset index in the "Sets" lists.
  ii <- NULL # iterator over the subsets
  jj <- NULL # iterator over the statistics
  kk <- NULL # iterator over the permutations
  
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
  
  # Since we expect the user to explicitly handle the number of parallel threads,
  # we will disable the potential implicit parallelism on systems where R has
  # been compiled against a multithreaded BLAS, e.g. OpenBLAS. 
  omp_set_num_threads(1)
  blas_set_num_threads(1)
  
  # Note to self: if we do any linear algebra on the full networks in 
  # preparation at a later date, we'll first need to set the number of threads
  # to `nCores`.
  
  nNets <- length(nodeLabelSets)
  
  # Set up names for more informative output, if the datasets and their 
  # corresponding networks have been given names.
  if (is.null(names(nodeLabelSets))) {
    setNames <- 1:nNets
  } else {
    setNames <- names(nodeLabelSets)
  }
  
  # Set up temporary directory for working big.matrix objects, and permutation
  # results.
  dir.create(".temp-objects", showWarnings=FALSE)
  on.exit({
    unlink(".temp-objects", recursive=TRUE)
  }, add=TRUE)
  
  scaledSets <- NULL
  if (!is.null(datSets)) {
    scaledSets <- rep(list(NULL), length(datSets))
  }
  
  # Set up return list 
  res <- rep(list(NULL), nNets)
  res <- lapply(res, function(x) { 
    l <- rep(list(NULL), nNets)
    names(l) <- names(nodeLabelSets)
    l
  })
  names(res) <- names(nodeLabelSets)
  
  # Iterate pairwise over datasets, comparing those marked "discovery"
  # with each marked as "test".
  for (di in seq_len(nNets)) {
    for (ti in seq_len(nNets)) {
      if ((di %in% discovery) & (ti %in% test) & (di != ti)) {
        # Set up return list
        res[[di]][[ti]] <- rep(list(NULL), 5)
        
        # Output messages
        vCat(verbose, indent, sep="", 
             "Calculating preservation of network subsets from dataset ",
             setNames[di], ", in dataset ", setNames[ti], ".")

        # Attach relevant matrices
        vCat(verbose, indent+1, "Attaching matrices...")
        discCor <- attach.big.matrix(corSets[[di]])
        discAdj <- attach.big.matrix(adjSets[[di]])
        testCor <- attach.big.matrix(corSets[[ti]])
        testAdj <- attach.big.matrix(adjSets[[ti]])
        if (!is.null(datSets[[di]]) && !is.null(datSets[[ti]])) {
          discDat <- attach.big.matrix(datSets[[di]])
          testDat <- attach.big.matrix(datSets[[ti]])
        } else {
          discDat <- NULL
          testDat <- NULL
        }        
        
        vCat(verbose, indent+1, "Checking matrices...")
        checkFinite(discCor)
        checkFinite(discAdj)
        checkFinite(testCor)
        checkFinite(testAdj)
        if (!is.null(datSets[[di]]) && !is.null(datSets[[ti]])) {
          checkFinite(discDat)
          checkFinite(testDat)
        }
        
        # Create scaled data 
        if (!is.null(discDat) && !is.null(testDat)) {
          vCat(verbose, indent+1, "Creating temporary scaled datasets for kME",
               "calculations...")
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
        } else {
          scaledDisc <- NULL
          scaledTest <- NULL
        }
        
        # Get a vector of nodes which are present in both datasets. Depends on 
        # the combination of data input provided.
        vCat(verbose, indent+1, "Extracting information about node overlap...")
        oNodes <- nodeNameSets[[di]] %sub_in% nodeNameSets[[ti]]
        tNodes <- nodeNameSets[[ti]]
        
        # Force the node indices to be strings, even if the provided identifiers
        # are integers.
        oNodes <- as.character(oNodes)
        tNodes <- as.character(tNodes)
        
        if (length(oNodes) == 0) {
          warning("No nodes in dataset ", setNames[di],  
                  " are present in dataset ", setNames[ti], ", skipping.")
          next
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
        nodesPres <- table(nodeLabelSets[[di]][oNodes])
        nodesPres <- nodesPres[names(nodesPres) %in% oSubsets]
        # Need to handle the case where a subset of interest has no overlap.
        nodesPres <- c(nodesPres, rep(0, length(dSubsets %sub_nin% oSubsets)))
        names(nodesPres) <- c(names(nodesPres) %sub_nin% "", dSubsets %sub_nin% oSubsets)
        propNodesPres <- nodesPres[names(dSizes)]/dSizes
        
        # How many network properites and statistics will we return? 
        # Numbers required for data structure allocation
        nStats <- ifelse(is.null(discDat), 4, 7)
        nSubsets <- length(oSubsets)
        
        # Calculate some basic cross-tabulation statistics so we can assess
        # which modules in both datasets map to each other, if module detection
        # has also been performed for the test network
        if (!is.null(nodeLabelSets[[ti]])) {
          # Get total number of nodes from each discovery subset in each test subset 
          subsetOverlap <- table(
            nodeLabelSets[[di]][oNodes], 
            nodeLabelSets[[ti]][oNodes]
          )
          # filter on subsets the user cares about
          subsetOverlap <- subsetOverlap[dSubsets,]
          # order the tables
          tryCatch({
            # For modules that are integer coded, make sure they're numerically
            # ordered, not alphabetically.
            rOrder <- order(as.integer(rownames(subsetOverlap)))
            cOrder <- order(as.integer(colnames(subsetOverlap)))
          }, warning = function(w) {
            # If we can't cast to an integer, sort normally.
            rOrder <- order(rownames(subsetOverlap))
            cOrder <- order(colnames(subsetOverlap))
          })
          subsetOverlap <- subsetOverlap[rOrder, cOrder]
          # add information about sizes of both the discovery and test subsets
          dSubSizes <- table(nodeLabelSets[[di]][oNodes])[rOrder]
          tSubSizes <- table(nodeLabelSets[[ti]][oNodes])[cOrder]
          subsetOverlap <- cbind(dSubSizes, subsetOverlap)
          subsetOverlap <- rbind(c(NA, tSubSizes), subsetOverlap)
          rownames(subsetOverlap)[1] <- "size"
          colnames(subsetOverlap)[1] <- "size"
        } else {
          subsetOverlap <- NULL
        }
        
        # Obtain the topological properties for each network subset in the
        # discovery dataset, we only want to calculate these once!
        vCat(verbose, indent+1, "Calculating observed test statistics...")
        poke(discCor, discAdj, discDat, scaledDisc, testCor, testAdj, testDat, scaledTest)
        discProps <- rep(list(NULL), nSubsets)
        for (ii in seq_along(oSubsets)) {
          sNodes <- names(which(nodeLabelSets[[di]][oNodes] == oSubsets[ii]))
          # get the indices in the underlying data and adjacency matrices for 
          # the subset nodes.
          subsetInd <- match(sNodes, nodeNameSets[[di]])
          props <- subsetProps(
            discAdj, subsetInd, discDat, scaledDisc
          )
          discProps[[ii]] <- props
        }
        names(discProps) <- oSubsets
        
        # Now calculate the observed value for each network statistic
        observed <- matrix(NA, nrow=nSubsets, ncol=nStats)
        rownames(observed) <- oSubsets
        for (ii in seq_along(oSubsets)) {
          ss <- as.character(oSubsets[ii])
          sNodes <- names(which(nodeLabelSets[[di]][oNodes] == ss))
          testInd <- match(sNodes, nodeNameSets[[ti]])
          discInd <- match(sNodes, nodeNameSets[[di]])
          testProps <- subsetProps(
            testAdj, testInd, testDat, scaledTest
          )
          stats <- c(
            calcSplitTestStats(discProps[[ss]], testProps),
            calcSharedTestStats(discCor, discInd, testCor, testInd)
          )
          observed[ii,] <- stats
          colnames(observed) <- names(stats)
        }
        
        # Calculate the null distribution for each of the statistics.
        vCat(verbose, indent+1, "Calculating null distributions with", nPerm, 
             "permutations...")
        if(verbose) {
          # To log progress, we will write our progress to a file for each chunk
          dir.create("run-progress", showWarnings=FALSE)
        }
        foreach(chunk=ichunkTasks(verbose, nPerm, nCores)) %maybe_do_par% {
          if (verbose & length(chunk) == 1) {
            if (chunk == -1) {
              monitorProgress(nWorkers, indent+2)
              NULL
            }
          } else {
            poke(discCor, discAdj, discDat, scaledDisc, testCor, testAdj, testDat, scaledTest)
            if (verbose) {
              conns <- setupParProgressLogs(chunk, nWorkers, indent+2)
              progressBar <- conns[[1]]
              on.exit(lapply(conns, close))
            } 
            chunkStats <- array(NA, dim=c(nSubsets, nStats, length(chunk)))
            dimnames(chunkStats)[1:2] <- dimnames(observed)
            dimnames(chunkStats)[[3]] <- paste0("permutation.", chunk)
            for (kk in seq_along(chunk)) {
              gc() # avoid memory leaks by garbage collecting at each iteration.
              for (ii in seq_along(oSubsets)) {
                ss <- as.character(oSubsets[ii])
                # Select a random subset of nodes of the same size as the subset 
                # ss, depending on our null model.
                if (model == "overlap") {
                  permNames <- sample(oNodes, size=nodesPres[ss])
                } else {
                  permNames <- sample(tNodes, size=nodesPres[ss])
                }
                
                permInd <- match(permNames, nodeNameSets[[ti]])
                sNodes <- names(which(nodeLabelSets[[di]][oNodes] == ss))
                discInd <- match(sNodes, nodeNameSets[[di]])
                
                tryCatch({ 
                  testProps <- subsetProps(
                    testAdj, permInd, testDat, scaledTest
                  )
                  stats <- c(
                    calcSplitTestStats(discProps[[ss]], testProps),
                    calcSharedTestStats(discCor, discInd, testCor, permInd)
                  )
                  chunkStats[ii,,kk] <- stats
                }, error = function(e) {
                  warning(
                    "Calculation for subset ", oSubsets[ii], " failed on ",
                    "permutation ", chunk[kk], " with error message:\n",
                    e$message
                  )
                })
              }
              # Update the progress at the end of the loop.
              if (verbose) {
                updateParProgress(progressBar, chunk[kk])
                if (nCores == 1) {
                  reportProgress(indent+2)
                  if (chunk[kk] == nPerm) {
                    cat("\n")
                  }
                }
              }
            }
            chunkNum <- ceiling(chunk[1]/length(chunk))
            permFile <- paste0("chunk", chunkNum, "permutations.rds")
            saveRDS(chunkStats, file.path(".temp-objects", permFile))
          }
        }
        # Load in results
        nulls <- array(NA, dim=c(nSubsets, nStats, nPerm))
        dimnames(nulls)[[3]] <- rep("", dim(nulls)[3])
        chunkFiles <- list.files(".temp-objects", "chunk[0-9]*permutations.rds")
        offset <- 1
        for (cf in chunkFiles) {
          chunk <- readRDS(file.path(".temp-objects", cf))
          nCPerm <- dim(chunk)[3]
          nulls[,,offset:(offset+nCPerm-1)] <- chunk
          dimnames(nulls)[1:2] <- dimnames(chunk)[1:2]
          dimnames(nulls)[[3]][offset:(offset+nCPerm-1)] <- dimnames(chunk)[[3]]
          offset <- offset + nCPerm
        }        
        
        # Calculate the p-value for the observed statistic based on the null 
        # distribution
        vCat(verbose, indent+1, "Calculating P-values...")
        p.values <- matrix(NA, nrow=nSubsets, ncol=nStats)
        dimnames(p.values) <- dimnames(observed)
        for (ii in seq_along(oSubsets)) {
          for (jj in seq_len(nStats)) {
            if (colnames(observed)[jj] %in% c("mean.adj", "propVarExpl")) {
              alternative <- "greater"
              order <- FALSE
            } else {
              alternative <- "two.sided"
              order <- TRUE
            }
            p.values[ii, jj] <- perm.test(
              nulls[ii, jj, ], observed[ii, jj], 
              nodesPres[rownames(p.values)[ii]], length(oNodes),
              order=order, alternative=alternative
            )
          }
        }

        # Collate results
        # First order output nicely
        tryCatch({
          # For modules that are integer coded, make sure they're numerically
          # ordered, not alphabetically.
          arrOrder <- order(as.integer(rownames(nulls)))
          vOrder <- order(as.integer(names(propNodesPres)))
        }, warning = function(w) {
          # If we can't cast to an integer, sort normally.
          arrOrder <- order(rownames(nulls))
          vOrder <- order(propNodesPres)
        })
        
        # Order statistics: First density stats, then connectivity
        if (is.null(discDat)) {
          statOrder <- c("mean.adj", "cor.kIM", "cor.cor", "mean.cor")
        } else {
          statOrder <- c(
            "mean.adj", "propVarExpl", "cor.cor", "cor.kIM", "cor.kME",
            "mean.cor", "mean.kME"
          ) 
        }
        
        # Collate results
        res[[di]][[ti]][[1]] <- nulls[arrOrder, statOrder,]
        res[[di]][[ti]][[2]] <- observed[arrOrder, statOrder]
        res[[di]][[ti]][[3]] <- p.values[arrOrder, statOrder]
        res[[di]][[ti]][[4]] <- nodesPres[vOrder]
        res[[di]][[ti]][[5]] <- propNodesPres[vOrder]
        
        if(!is.null(subsetOverlap)) {
          res[[di]][[ti]][[6]] <- subsetOverlap
          names(res[[di]][[ti]]) <- c(
            "nulls", "observed", "p.values", "nodesPresent", "propNodesPresent",
            "contigency"
          )
        } else {
          names(res[[di]][[ti]]) <- c(
            "nulls", "observed", "p.values", "nodesPresent", "propNodesPresent"
          )
        }
        
        vCat(verbose, indent+1, "Cleaning up temporary objects...")
        unlink("run-progress", recursive=TRUE)
        rm(discDat, scaledDisc, discAdj, testDat, scaledTest, testAdj)
        gc()
        vCat(verbose, indent, "Done!")
      }
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

