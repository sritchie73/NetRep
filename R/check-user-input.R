#' Process, check, and format input
#' 
#' Checks user input for consistency, errors, and unifies the data structures.
#' 
#' @param discovery user input for the 'discovery' argument.
#' @param test user input for the 'test' argument.
#' @param network user input for the 'network' argument.
#' @param correlation user input for the 'correlation' argument.
#' @param data user input for the 'data' argument.
#' @param moduleAssignments user input for the 'moduleAssignments' argument.
#' @param modules user input for the 'modules' argument.
#' @param backgroundLabel user input for the 'backgroundLabel' argument.
#' @param verbose logical; should progress be reported? Default is \code{TRUE}.
#' @param plotFunction logical; are we checking for a plot function?
#' @param orderNodesBy user input for the 'orderNodesBy' argument in the 
#'  plotting functions.
#' @param orderSamplesBy user input for the 'orderSamplesBy' argument in the 
#'  plotting functions.
#' @param orderModules user input for the 'orderModules' argument in the 
#'  plotting functions.
#' 
#' @seealso
#' \code{\link{modulePreservation}}
#' \code{\link{plotModule}}, and
#' \code{\link{plotTopology}}
#' 
#' @return a list of containing the formatted user input
processInput <- function(
  discovery, test, network, correlation, data, moduleAssignments, modules, 
  backgroundLabel, verbose, plotFunction=FALSE, orderNodesBy=NA, 
  orderSamplesBy=NA, orderModules=NULL
) {
  # Where do we want to get:
  #   Each "argument" has a list of lists: at the top level, each element 
  #   corresponds to a "discovery" dataset, containing a list, where each element
  #   corresponds to a "test" dataset. 
  #
  # We want the user to be able to input data in whatever way makes the most 
  # sense for *their* analysis. So we need to process their input into the right
  # format.
  #  - Should not require the "list" wrapper if only element required (e.g.
  #    'moduleAssignments' if there is only 1 discovery dataset, or 'nCores'
  #    assuming we want that to be the same across arguments).
  #  - Datasets may be named, unamed, or a mixture of both. We want to infer the
  #    mapping between list elements across arguments where possible.
  #
  # Arguments to process for all user exposed functions:
  #  - 'data' 
  #  - 'correlation'
  #  - 'network'
  #  - 'moduleAssignments'
  #  - 'discovery'
  #  - 'test'
  #  - 'modules'
  #
  # Arguments we expect will usually only have 1 value, but may have more to 
  # be consistent with the above. These only occur in the 'modulePreservation'
  # function.
  #
  #  - 'nPerm'
  #

  # Because this is what I **actually** meant, obviously. Why does is.vector
  # return TRUE for lists???
  is.vector <- function(obj) {
    base::is.vector(obj) && !is.list(obj)
  }
  
  # ----------------------------------------------------------------------------
  # First, we need to know what the 'discovery' and 'test' datasets are.
  # ----------------------------------------------------------------------------
  
  # If both are not provided we assume the user is using one of the 
  # downstream analysis or plotting functions on a single dataset
  if (is.null(discovery) & is.null(test)) {
    discovery <- 1
    test <- 1
  }
  # If only one is provided, we assume the user meant to use one of the 
  # downstream analysis or plotting functions within the same dataset
  else if (is.null(discovery)) {
    discovery <- test
  } else if (is.null(test)) {
    test <- discovery
  }

  # Case 1: both discovery and test are vectors
  if (is.vector(discovery) && is.vector(test)) {
    # Are the datasets named?
    discNames <- NULL
    if (is.character(discovery))
      discNames <- discovery
    names(discovery) <- discNames
    
    if (is.numeric(discovery)) {
      tmp <- rep(list(NULL), max(discovery))
      tmp[discovery] <- rep(list(test), length(discovery))
      test <- tmp
    } else if (is.character(discovery)) {
      test <- rep(list(test), length(discovery))
      names(test) <- discNames
    } else {
      stop("'discovery' must be ", '"character" or "numeric"')
    }
  }
  # Case 2: discovery is a vector, test is a list (i.e. different test datasets
  # for each discovery dataset)
  else if (is.vector(discovery) && is.list(test)) {
    # Are the datasets named?
    discNames <- NULL
    if (is.character(discovery))
      discNames <- discovery
    names(discovery) <- discNames
    
    if (length(test) != length(discovery)) {
      stop("mismatch between 'discovery' and 'test' arguments")
    } else if (is.numeric(discovery)) {
      tmp <- rep(list(NULL), max(discovery))
      tmp[discovery] <- test
      test <- tmp
    }

    if (is.null(names(test))) {
      names(test) <- discNames
    }
      
    if (!is.null(names(test))) {
      if (names(test) %nin% discNames) {
        stop("mismatch between 'discovery' and 'test' arguments")
      }
    }
  } else {
    stop("incorrect data structures provided for 'discovery' or 'test' arguments")
  }
  
  # What are the dataset names, and/or what is the total number of datasets?
  dataNames <- NULL
  nDatasets <- 1
  matchByIndice <- FALSE # will we ever match by indice?
  if (is.character(discovery)) {
    dataNames <- unique(c(dataNames, discovery))
    nDatasets <- length(dataNames)
  } else if (is.numeric(discovery)) {
    nDatasets <- max(discovery)
    matchByIndice <- TRUE
  } else {
    stop("unable to match 'discovery' to provided datasets")
  }
  for (tv in test) {
    if (is.character(tv)) {
      dataNames <- unique(c(dataNames, tv))
      nDatasets <- length(dataNames)
    } else if (is.numeric(tv)) {
      nDatasets <- max(tv)
      matchByIndice <- TRUE
    } else if (!is.null(tv)) {
      stop("unable to match 'test' to provided datasets")
    }
  }
  
  # Make sure test and discovery are ordered the same way
  if (!is.null(names(discovery)))
    test <- test[names(discovery)]
  
  # Plots can only be generated within a single dataset at a time
  if (plotFunction) { 
    if ((!is.vector(discovery) || length(discovery) > 1) ||
        (!is.vector(test[[discovery]]) || length(test[[discovery]]) > 1)) {
      stop("only one 'discovery' and 'test' dataset can be specified when plotting")
    }
  }

  
  # ----------------------------------------------------------------------------
  # Next, process the 'correlation' and 'network' arguments
  # ----------------------------------------------------------------------------
  
  if (!is.list(correlation))
    correlation <- list(correlation)
  if (!is.list(network))
    network <- list(network)
  
  # Check if the data is in an appropriate format
  lapply(correlation, checkIsMatrix)
  lapply(network, checkIsMatrix)
  
  # Add any datasets names that are not in dataNames
  dataNames <- c(dataNames, names(correlation))
  dataNames <- c(dataNames, names(network))
  dataNames <- unique(dataNames)
  
  nDatasets <- max(c(nDatasets, length(dataNames), length(network)))
  
  # Check that we can match 'discovery' and 'test' to the provided matrices. 
  correlation <- verifyDatasetOrder(correlation, "correlation", dataNames, nDatasets)
  network <- verifyDatasetOrder(network, "network", dataNames, nDatasets)
  
  if (any(dataNames %nin% names(network)) || nDatasets != length(network)) {
    stop("mismatch between 'discovery', 'test', and the datasets provided")
  }
  
  # ----------------------------------------------------------------------------
  # Next, process the 'data' argument
  # ----------------------------------------------------------------------------
  
  # Handle special case where 'data' can be 'NULL'
  if (is.null(data)) {
    data <- rep(list(NULL), nDatasets)
    names(data) <- names(network)
  }

  # Otherwise check as per 'network' and 'correlation'
  if (!is.list(data))
    data <- list(data)

  # Check data is in appropriate format
  lapply(data, checkIsMatrix)

  # Check that we can match 'discovery' and 'test' to the provided matrices. 
  data <- verifyDatasetOrder(data, "data", dataNames, nDatasets)
  
  # Or if dataset names exist in one of the other input lists.
  if (!is.null(names(network))) {
    # If this fails, ignore for now and throw appropriate error when checking
    # data consistency.
    tryCatch({
      names(data) <- names(network)
    }, error = function(e) {}, warning=function(w) {})
  }
  
  # ----------------------------------------------------------------------------
  # Next, process the 'moduleAssignments' argument
  # ----------------------------------------------------------------------------
  
  # Handle cases where moduleAssignments is not provided.
  # Assume the user just wants to look at all nodes as a whole
  if (is.null(moduleAssignments)) {
    # Discovery datasets are named
    if (is.character(discovery)) {
      moduleAssignments <- lapply(discovery, function(di) {
        nodes <- colnames(network[[di]])
        structure(rep("1", length(nodes)), names=nodes)
      })
    }
    # Discovery datasets are referred to by index
    else if (is.numeric(discovery)) {
      moduleAssignments <- rep(list(NULL), max(discovery))
      for (di in discovery) {
        nodes <- colnames(network[[di]])
        moduleAssignments[[di]] <- structure(rep("1", length(nodes)), names=nodes)
      }
    }
    else {
      stop("unexpected error when automatically constructing",
           " 'moduleAssignments' object")
    }
  }
  
  # Handle cases where moduleAssignments assumed to be for discovery dataset
  if (!is.list(moduleAssignments)) 
    moduleAssignments <- list(moduleAssignments)
  
  if (length(moduleAssignments) < length(network)) {
    tmp <- rep(list(NULL), length(network))
    names(tmp) <- names(network)
    if (length(discovery) != length(moduleAssignments))
      stop("must have a 'moduleAssignments' vector for each 'discovery' dataset")
    tmp[discovery] <- moduleAssignments
    moduleAssignments <- tmp
  }
    
  # Make sure that its a list of vectors
  for (ii in seq_along(moduleAssignments)) {
    if (!is.null(moduleAssignments[[ii]]) && !is.vector(moduleAssignments[[ii]]))
      stop("expecting a list of vectors for 'moduleAssignments'")
  }
  
  # Check that we can match 'discovery' to the provided 'moduleAssignments'
  discNames <- names(discovery)
  if(!is.null(discNames) && is.null(names(moduleAssignments))) {
    stop("cannot match dataset names in 'discovery' to the provided ",
         "'moduleAssignments' list")
  }
  if (!is.null(discNames) && any(discNames %nin% names(moduleAssignments))) {
    stop("cannot match dataset names in 'discovery' to the provided ",
         "'moduleAssignments' list")
  }
  if (is.null(discNames) && (length(moduleAssignments) < length(unique(discovery)))) {
    stop("expecting ", length(discovery), " 'moduleAssignment' vectors, ", 
         nDatasets, " provided")
  }

  # Make sure numeric labels are transformed to character labels for matching
  # with 'modules'
  for (ii in seq_along(moduleAssignments)) {
    modVec <- moduleAssignments[[ii]]
    if (!is.null(modVec)) {
      moduleAssignments[[ii]] <- structure(as.character(modVec), names=names(modVec))
    }
  }
  
  # ----------------------------------------------------------------------------
  # Next, process the 'backgroundLabel' argument
  # ----------------------------------------------------------------------------
  
  # User had to explicitly turn off, they probably mean "don't ignore this 
  # module"
  if (is.null(backgroundLabel))
    backgroundLabel <- vector()

  # Make one for each discovery dataset if its a vector
  if (is.vector(backgroundLabel)) {
    tmp <- rep(list(NULL), length(moduleAssignments))
    names(tmp) <- names(moduleAssignments)
    for (ii in seq_along(moduleAssignments)) {
      if (!is.null(moduleAssignments[[ii]])) {
        tmp[[ii]] <- backgroundLabel
      }
    }
    backgroundLabel <- tmp
  }
  
  # if passed as a list, they probably mean one for each discovery dataset 
  if (length(backgroundLabel) < length(moduleAssignments)) {
    tmp <- rep(list(NULL), length(moduleAssignments))
    names(tmp) <- names(moduleAssignments)
    if (length(discovery) != length(backgroundLabel))
      stop("must have a 'backgroundLabel' vector for each 'discovery' dataset")
    tmp[discovery] <- backgroundLabel
    backgroundLabel <- tmp
  }
  
  # Make sure that its a list of vectors
  for (ii in seq_along(backgroundLabel)) {
    if (!is.null(backgroundLabel[[ii]]) && !is.vector(backgroundLabel[[ii]]))
      stop("expecting a list of vectors for 'backgroundLabel'")
  }
  
  # Check that we can match 'discovery' to the provided 'backgroundLabel'
  discNames <- names(discovery)
  if(!is.null(discNames) && is.null(names(backgroundLabel))) {
    stop("cannot match dataset names in 'discovery' to the provided ",
         "'backgroundLabel' list")
  }
  if (!is.null(discNames) && any(discNames %nin% names(backgroundLabel))) {
    stop("cannot match dataset names in 'discovery' to the provided ",
         "'backgroundLabel' list")
  }
  if (is.null(discNames) && (length(backgroundLabel) < length(unique(discovery)))) {
    stop("expecting ", length(discovery), " 'moduleAssignment' vectors, ", 
         nDatasets, " provided")
  }
  
  # If any nodes are missing module assignments, set them to the background 
  # label
  for (di in discovery) {
    unlabelled <- colnames(network[[di]]) %sub_nin% names(moduleAssignments[[di]])
    if (length(unlabelled) > 0) {
      if (length(backgroundLabel[[di]]) > 0) {
        bglabel <- backgroundLabel[[di]][1]
      } else {
        # We need a label that doesn't conflict with any existing module labels.
        # This is the simplest way to do it.
        bglabel <- as.integer(Sys.time())
      }
      
      bgnodes <- rep(bglabel, length(unlabelled))
      names(bgnodes) <- unlabelled
      moduleAssignments[[di]] <- c(moduleAssignments[[di]], bgnodes)
    }
  }
  
  # ----------------------------------------------------------------------------
  # Next, process the 'modules' argument
  # ----------------------------------------------------------------------------

  # If not specified, run for all modules except the network background
  if (is.null(modules)) {
    # Discovery datasets are named
    if (is.character(discovery)) {
      modules <- lapply(moduleAssignments, function(ma) {
        mods <- names(table(ma))
        mods[mods %nin% backgroundLabel]
      })
    }
    # Discovery datasets are referred to by index
    else if (is.numeric(discovery)) {
      modules <- rep(list(NULL), length(moduleAssignments))
      for (ii in seq_along(modules)) {
        if (!is.null(moduleAssignments[[ii]])) {
          mods <- names(table(moduleAssignments[[ii]]))
          mods <- mods[mods %nin% backgroundLabel]
          modules[[ii]] <- mods
        }
      }
    }
    else {
      stop("unexpected error when automatically constructing 'modules' object")
    }
    
    # If modules labels are numeric, sort them
    modules <- lapply(modules, function(modVec) {
      tryCatch({
        numVec <- as.numeric(modVec)
        return(modVec[order(numVec)])
      }, warning=function(w) {
        # Cant cast to numeric, return as is
        return(modVec)
      })
    })
    
    # Name if dataset names exist
    if (!is.null(names(moduleAssignments)))
      names(modules) <- names(moduleAssignments)
    
    # Or if dataset names exist in one of the other input lists.
    if (!is.null(names(network))) {
      # If this fails, ignore for now and throw appropriate error when checking
      # data consistency.
      tryCatch({
        names(modules) <- names(network)
      }, error = function(e) {}, warning=function(w) {})
    }
  }
  
  # Handle cases where modules assumed to be for discovery dataset
  if (!is.list(modules)) 
    modules <- list(modules)
  
  if (length(modules) < length(network)) {
    tmp <- rep(list(NULL), length(network))
    names(tmp) <- names(network)
    if (length(discovery) != length(modules))
      stop("must have a 'modules' vector for each 'discovery dataset")
    tmp[discovery] <- modules
    modules <- tmp
  }
  
  # Make sure that its a list of vectors
  for (ii in seq_along(modules)) {
    if (!is.null(modules[[ii]]) && !is.vector(modules[[ii]]))
      stop("expecting a list of vectors for 'modules'")
  }
  
  # Check that we can match 'discovery' to the provided 'moduleAssignments'
  discNames <- names(discovery)
  if(!is.null(discNames) && is.null(names(modules))) {
    stop("cannot match dataset names in 'discovery' to the provided ",
         "'modules' list")
  }
  if (!is.null(discNames) && any(discNames %nin% names(modules))) {
    stop("cannot match dataset names in 'discovery' to the provided ",
         "'modules' list")
  }
  if (is.null(discNames) && (length(modules) < nDatasets)) {
    stop("expecting ", nDatasets, " 'modules' vectors, ", 
         length(modules), " provided")
  }
  
  # ----------------------------------------------------------------------------
  # Next, process the plot function arguments
  # ----------------------------------------------------------------------------
  
  if (plotFunction) {
    if (!(
      is.null(orderNodesBy) ||
      is.vector(orderNodesBy) && is.numeric(orderNodesBy) ||
      is.vector(orderNodesBy) && is.character(orderNodesBy) ||
      is.vector(orderNodesBy) && length(orderNodesBy) == 1 && is.na(orderNodesBy)
    )) {
      stop("'orderNodesBy' must be a vector of dataset names or indices, 'NA' or",
           " 'NULL'")
    }
    
    if (!(
      is.null(orderSamplesBy) ||
      is.vector(orderSamplesBy) && length(orderSamplesBy) == 1 && 
      (is.numeric(orderSamplesBy) || is.character(orderSamplesBy) || is.na(orderSamplesBy))
    )) {
      stop("'orderSamplesBy' must be a vector containing a single dataset name, ",
           "or index, 'NA' or 'NULL'")
    }
    
    if (is.null(orderNodesBy))
      orderNodesBy <- discovery
    if (is.null(orderSamplesBy)) {
      if (is.null(data[[discovery]])) {
        orderSamplesBy <- NA
      } else {
        orderSamplesBy <- test[[discovery]]
      }
    }

    ti <- test[[discovery]]
    # Are the datasets specified in 'orderNodesBy' and 'orderSamplesBy' valid?
    if (length(orderNodesBy) > 1 || !is.na(orderNodesBy)) {
      if (is.character(orderNodesBy) && any(orderNodesBy %nin% names(network))) {
        stop("unable to match datasets in 'orderNodesBy' to provided datasets")
      } else if (is.numeric(orderNodesBy) && any(orderNodesBy > nDatasets) || 
                 any(orderNodesBy < 1)) {
        stop("unable to match datasets in 'orderNodesBy' to provided datasets")
      }
    }
    if (!is.na(orderSamplesBy)) {
      if (is.character(orderSamplesBy) && any(orderSamplesBy %nin% names(network))) {
        stop("unable to match datasets in 'orderSamplesBy' to provided datasets")
      } else if (is.numeric(orderSamplesBy) && (orderSamplesBy > nDatasets || orderSamplesBy < 1)) {
        stop("unable to match datasets in 'orderSamplesBy' to provided datasets")
      }
    }
    
    if (!is.na(orderSamplesBy) && is.null(data[[orderSamplesBy]])) {
      stop("'data' not provided for 'orderSamplesBy' dataset") 
    }
    
    # Check that data is provided for the 'orderNodesBy' dataset(s) if 
    # 'orderModules' is true.
    if ((orderModules && length(modules) > 1) && 
        (length(orderNodesBy) > 1 || !is.na(orderNodesBy)) && 
        any(sapply(data[orderNodesBy], is.null))) {
      stop("'data' not provided for 'orderNodesBy' dataset(s) and ",
           "'orderModules' = 'TRUE'") 
    }
  }
  
  # ----------------------------------------------------------------------------
  # Construct consistent names to use internally for each dataset
  # ----------------------------------------------------------------------------
  # Make sure the input data can be sensibly accessed using the provided 
  # 'discovery' and 'test' inputs.
  
  # If one, but not all, input lists are named, throw an error -- we don't know
  # what the datasets are!
  datNames <- names(data)
  corNames <- names(correlation)
  netNames <- names(network)
  modLabNames <- names(moduleAssignments)
  modNames <- names(modules)
  
  hasNames <- c(!is.null(datNames), !is.null(corNames), !is.null(netNames), 
                !is.null(modLabNames), !is.null(modNames))
  
  if (sum(hasNames) > 0 && !all(hasNames)) {
    stop("cannot match dataset names across all input arguments")
  }
  
  # Construct the dataset names to use as indices
  if (!is.null(names(network))) {
    datasetNames <- names(network)
  } else {
    datasetNames <- paste0("Dataset", seq_len(nDatasets))
    names(data) <- datasetNames
    names(modules) <- datasetNames
    names(correlation) <- datasetNames
    names(network) <- datasetNames
    names(moduleAssignments) <- datasetNames
  }
  names(datasetNames) <- datasetNames
  
  # Convert indices to dataset names for the plot functions
  if (plotFunction) {
    if (is.numeric(orderNodesBy)) {
      orderNodesBy <- datasetNames[orderNodesBy]
    }
    if (is.numeric(orderSamplesBy)) {
      orderSamplesBy <- datasetNames[orderSamplesBy]
    }
  }

  # ----------------------------------------------------------------------------
  # Check for data consistency
  # ----------------------------------------------------------------------------
  # Construct an iterator that includes only the datasets we're analysing:
  if (is.character(discovery)) {
    iterator <- match(discovery, datasetNames)
  } else {
    iterator <- discovery
  }
  for (tv in test) {
    if (is.character(tv)) {
      iterator <- c(iterator, match(tv, datasetNames))
    } else {
      iterator <- c(iterator, tv)
    }
  }
  if (plotFunction) {
    if (orderModules) {
      if (is.character(orderNodesBy)) {
        iterator <- c(iterator, match(orderNodesBy, datasetNames))
      } else if (is.numeric(orderNodesBy)) {
        iterator <- c(iterator, orderNodesBy)
      }
    }
    if (is.character(orderSamplesBy)) {
      iterator <- c(iterator, match(orderSamplesBy, datasetNames))
    } else if (is.numeric(orderSamplesBy)) {
      iterator <- c(iterator, orderSamplesBy)
    }
  }
  iterator <- datasetNames[unique(iterator)]
  
  # We want to iterate over the first discovery dataset last, so that we
  # can skip loading it in a second time in the calling function
  tokeep <- discovery[1]
  iterator <- c(iterator[-which(iterator == tokeep)], tokeep)
  
  # We need a list of nodes present in each dataset independent of having
  # the datasets loaded into RAM.
  nodelist <- rep(list(NULL), length(iterator))
  names(nodelist) <- datasetNames[iterator]
  
  # If plotting, we need to check that samples from the 'orderSamplesBy' 
  # dataset are present in the 'test' dataset to be drawn.
  if (plotFunction) {
    pIdx <- test[[discovery]]
    sIdx <- orderSamplesBy
    pSamples <- NULL
    sSamples <- NULL
  }
  
  anyDM <- do.call("any.disk.matrix", c(data[iterator], correlation[iterator], network[iterator]))
  vCat(verbose && !anyDM, 1, "Checking matrices for problems...")
  for (ii in iterator) {
    # First, we need to load in matrices into RAM if they are 'disk.matrix' 
    # objects
    anyDM <- any.disk.matrix(data[[ii]], correlation[[ii]], network[[ii]])
    vCat(verbose && anyDM, 1, 'Loading matrices of dataset "', 
         datasetNames[ii], '" into RAM...', sep="")
    dataLoaded <- loadIntoRAM(data[[ii]])
    correlationLoaded <- loadIntoRAM(correlation[[ii]])
    networkLoaded <- loadIntoRAM(network[[ii]])
    
    vCat(verbose && anyDM, 1, "Checking matrices for problems...")
    
    # If plotting, we need to check that samples from the 'orderSamplesBy' 
    # dataset are present in the 'test' dataset to be drawn.
    if (plotFunction) {
      if (ii == pIdx) pSamples <- rownames(dataLoaded)
      if (ii == sIdx) sSamples <- rownames(dataLoaded)
      if (!is.null(pSamples) && !is.null(sSamples)) {
        if(length(intersect(pSamples, sSamples)) == 0) {
          stop("no samples in the dataset specified by 'orderSamplesBy' ", 
               "are in the 'test' dataset to be drawn.")
        }
      }
    }
    
    # Make sure matrices are (a) actually matrices, and (b) contain numeric 
    # data
    if (class(networkLoaded) != "matrix" || 
        typeof(networkLoaded) %nin% c("double", "integer")) {
      stop("'network' for dataset ", '"', ii, '"', 
           " is not a numeric matrix")
    }
    if (class(correlationLoaded) != "matrix" || 
        typeof(correlationLoaded) %nin% c("double", "integer")) {
      stop("'correlation' for dataset ", '"', ii, '"', 
           " is not a numeric matrix")
    }
    
    if (!is.null(dataLoaded) && (class(dataLoaded) != "matrix" || 
                               typeof(dataLoaded) %nin% c("double", "integer"))) {
      stop("'data' for dataset ", '"', ii, '"', " is not a numeric matrix")
    }
    
    # Make sure the 'correlation' and 'network' matrices are square
    if (nrow(networkLoaded) != ncol(networkLoaded)) {
      stop("'network' for dataset ", '"', ii, '"', " is not square")
    }
    if (nrow(correlationLoaded) != ncol(correlationLoaded)) {
      stop("'correlation' for dataset ", '"', ii, '"', " is not square")
    }
    # And that they have the same dimensions
    if ((nrow(correlationLoaded) != nrow(networkLoaded)) ||
        (!is.null(dataLoaded) && (ncol(dataLoaded) != ncol(networkLoaded)))) {
      stop("'correlation', 'network', and 'data' have a different number of ",
           'nodes for dataset "', ii, '"')
    }
    
    # Make sure the matrices have dimension names
    if (is.null(rownames(networkLoaded)) || is.null(rownames(correlationLoaded)) ||
        (!is.null(dataLoaded) && is.null(rownames(dataLoaded)))) {
      stop("supplied matrices must have row and column names")      
    }
    
    # Make sure the 'correlation' and 'network' matrices are symmetric 
    if (any(rownames(networkLoaded) != colnames(networkLoaded))) {
      stop("mismatch between row and column names in 'network' for dataset ", 
           '"', ii, '"')
    }
    if (any(rownames(correlationLoaded) != colnames(correlationLoaded))) {
      stop("mismatch between row and column names in 'network' for dataset ",
           '"', ii, '"')
    }
    # Make sure the ordering of nodes is the same between 'correlation', 
    # 'network' and 'data'.
    if (any(colnames(networkLoaded) != colnames(correlationLoaded)) |
        (!is.null(dataLoaded) && any(colnames(networkLoaded) != colnames(dataLoaded)))) {
      stop("mismatch in node order between 'correlation' and 'network' for",
           'dataset "', ii, '"')
    }
    
    # Make sure the 'moduleAssignments' have the same nodes as the 'correlation'
    # etc.
    if (!is.null(moduleAssignments[[ii]])) {
      if (any(names(moduleAssignments[[ii]]) %nin% colnames(networkLoaded))) {
        stop("module assigments are present for nodes that are not in the",
             " 'network' inferred from dataset ", '"', ii, '"')
      }
    }
    
    # Make sure all module labels are in 'moduleAssignments'
    if (any(modules[[ii]] %nin% moduleAssignments[[ii]])) {
      stop("some 'modules' specified for dataset ", '"', ii, '"', 
           " are not in 'moduleAssignments' for dataset ", '"', ii, '"')
    }
    
    # Check matrices for non-finite values: these will cause the calculation
    # of network properties and module preservation statistics to hang. 
    if (!is.null(dataLoaded))
      CheckFinite(dataLoaded)
    CheckFinite(correlationLoaded)
    CheckFinite(networkLoaded)
    
    # Store the node names for later
    nodelist[[datasetNames[ii]]] <- colnames(networkLoaded)
    
    # Free up memory if any objects are big matrices, but return the last
    # dataset to pass to the calling function.
    if (ii != tokeep) {
      vCat(verbose && anyDM, 1, "Unloading dataset from RAM...")
      rm(dataLoaded, correlationLoaded, networkLoaded)
      gc()
    }
  }

  return(list(
    data=data, correlation=correlation, network=network, discovery=discovery,
    test=test, moduleAssignments=moduleAssignments, modules=modules,
    nDatasets=nDatasets, datasetNames=datasetNames,
    orderNodesBy=orderNodesBy, orderSamplesBy=orderSamplesBy,
    nodelist=nodelist, loadedIdx=tokeep, dataLoaded=dataLoaded, 
    correlationLoaded=correlationLoaded, networkLoaded=networkLoaded
  ))
}

#' Verify a 'list' input ordering
#' 
#' Check and order an input list: 
#' 
#' @param tocheck list of data to check.
#' @param errname name to print in error messages.
#' @param dataNames names of the datasets.
#' @param nDatasets number of datasets.
#' 
#' @return ordered 'tocheck' by dataset.
verifyDatasetOrder <- function(tocheck, errname, dataNames, nDatasets) {
  # Check that we can match 'discovery' and 'test' to the provided matrices
  if(!is.null(dataNames) && is.null(names(tocheck))) {
    stop("cannot match dataset names in 'discovery' and 'test' to the provided ",
         "'", errname, "' matrices")
  }
  if (!is.null(dataNames) && any(dataNames %nin% names(tocheck))) {
    stop("cannot match dataset names in 'discovery' and 'test' to the provided" ,
         "'", errname, "' matrices")
  }
  if (is.null(dataNames) && (length(tocheck) < nDatasets)) {
    stop("expecting ", nDatasets, "'", errname, "' matrices ", length(tocheck),
         " provided")
  }
  
  return(tocheck)
}

#' Check whether an object is a 'matrix' or a 'disk.matrix'
#' 
#' @param object object to check.
#' 
#' @return 
#'   throws an error or returns silently
checkIsMatrix <- function(object) {
  if (!is.null(object) && !is.matrix(object) && class(object) != "disk.matrix") { 
    stop('Input data must be a "matrix" or "disk.matrix"')
  }
}
  
#' Validate plot function arguments
#' 
#' Simple typechecking for the extensive plot arguments
#' 
#' @param orderModules user input for the corresponding argument in the plot functions.
#' @param plotNodeNames user input for the corresponding argument in the plot functions.
#' @param plotSampleNames user input for the corresponding argument in the plot functions.
#' @param plotModuleNames user input for the corresponding argument in the plot functions.
#' @param main user input for the corresponding argument in the plot functions.
#' @param drawBorders user input for the corresponding argument in the plot functions.
#' @param lwd user input for the corresponding argument in the plot functions.
#' @param naxt.line user input for the corresponding argument in the plot functions.
#' @param saxt.line user input for the corresponding argument in the plot functions.
#' @param maxt.line user input for the corresponding argument in the plot functions.
#' @param xaxt.line user input for the corresponding argument in the plot functions.
#' @param yaxt.line user input for the corresponding argument in the plot functions.
#' @param laxt.line user input for the corresponding argument in the plot functions.
#' @param xlab.line user input for the corresponding argument in the plot functions.
#' @param ylab.line user input for the corresponding argument in the plot functions.
#' @param main.line user input for the corresponding argument in the plot functions.
#' @param xaxt.tck user input for the corresponding argument in the plot functions.
#' @param yaxt.tck user input for the corresponding argument in the plot functions.
#' @param laxt.tck user input for the corresponding argument in the plot functions.
#' @param plotLegend user input for the corresponding argument in the plot functions.
#' @param legend.position user input for the corresponding argument in the plot functions.
#' @param legend.main user input for the corresponding argument in the plot functions.
#' @param legend.main.line input for the corresponding argument in the plot functions.
#' @param symmetric user input for the corresponding argument in the plot functions.
#' @param horizontal user input for the corresponding argument in the plot functions.
#' @param dataCols user input for the corresponding argument in the plot functions.
#' @param dataRange user input for the corresponding argument in the plot functions.
#' @param corCols user input for the corresponding argument in the plot functions.
#' @param corRange user input for the corresponding argument in the plot functions.
#' @param netCols user input for the corresponding argument in the plot functions.
#' @param netRange user input for the corresponding argument in the plot functions.
#' @param degreeCol user input for the corresponding argument in the plot functions.
#' @param contribCols user input for the corresponding argument in the plot functions.
#' @param summaryCols user input for the corresponding argument in the plot functions.
#' @param naCol user input for the corresponding argument in the plot functions.
#' @param dryRun user input for the corresponding argument in the plot functions.
#' 
checkPlotArgs <- function(
  orderModules, plotNodeNames, plotSampleNames, plotModuleNames, main,
  drawBorders, lwd, naxt.line, saxt.line, maxt.line, xaxt.line, 
  yaxt.line, laxt.line, xaxt.tck, yaxt.tck, laxt.tck, xlab.line, ylab.line,
  main.line, plotLegend, legend.position, legend.main, legend.main.line, 
  symmetric, horizontal, dataCols, dataRange, corCols, corRange, netCols, 
  netRange, degreeCol, contribCols, summaryCols, naCol, dryRun
) {
  # Return TRUE only if a an object is a vector, not a list.
  is.vector <- function(obj) {
    base::is.vector(obj) && !is.list(obj)
  }
  
  # Makes sure the check does not throw a warning if the vector has length > 1.
  is.na <- function(obj) {
    is.vector(obj) && length(obj) == 1 && base::is.na(obj) 
  }
  
  # Return TRUE if an argument is a numeric vector of length 1.
  is.snum <- function(obj) {
    is.vector(obj) && length(obj) == 1 && is.numeric(obj) && !is.na(obj)
  }
  
  # Return TRUE if an argument is a character vector of length 1.
  is.schar <- function(obj) {
    is.vector(obj) && length(obj) == 1 && is.character(obj) && !is.na(obj)
  }
  
  # Return TRUE if an argument is a logical vector of length 1 
  is.slog <- function(obj) {
    is.vector(obj) && length(obj) == 1 && is.logical(obj) && !is.na(obj)
  }
  
  if (!(missing(orderModules) || is.slog(orderModules)))
    stop("'orderModules' must be one of 'TRUE' or 'FALSE'")
  
  if (!(missing(plotNodeNames) || is.slog(plotNodeNames)))
    stop("'plotNodeNames' must be one of 'TRUE' or 'FALSE'")
  
  if (!(missing(plotSampleNames) || is.slog(plotSampleNames)))
    stop("'plotNodeNames' must be one of 'TRUE' or 'FALSE'")
  
  if (!(missing(plotModuleNames) || is.null(plotModuleNames) 
        || is.slog(plotModuleNames)))
    stop("'plotModuleNames' must be one of 'TRUE', 'FALSE', or 'NULL'")
  
  if (!(missing(main) || is.null(main) || is.schar(main) || is.na(main)))
    stop("'main' must be 'NULL' or a character vector of length 1")
  
  if (!(missing(drawBorders) || is.slog(drawBorders)))
    stop("'drawBorders' must be one of 'TRUE' or 'FALSE'")
  
  if (!missing(lwd)) {
    if (!is.snum(lwd)) {
      stop("'lwd' must be a numeric vector of length 1")
    }
    if (lwd < 0) {
      stop("'lwd' must be greater than 0")
    }
    if (is.infinite(lwd)) {
      stop("'lwd' must be finite")
    }
  }
  
  if (!missing(naxt.line)) {
    if (!(is.snum(naxt.line) || is.na(naxt.line) || is.null(naxt.line))) {
      stop("'naxt.line' must be a numeric vector of length 1, 'NA', or 'NULL'")
    }
    if (is.snum(naxt.line) && is.infinite(naxt.line)) {
      stop("'naxt.line' must be finite")
    }
  }
  
  if (!missing(saxt.line)) {
    if (!(is.snum(saxt.line) || is.na(saxt.line) || is.null(saxt.line))) {
      stop("'saxt.line' must be a numeric vector of length 1, 'NA', or 'NULL'")
    }
    if (is.snum(saxt.line) && is.infinite(saxt.line)) {
      stop("'saxt.line' must be finite")
    }
  }
  
  if (!missing(maxt.line)) {
    if (!(is.snum(maxt.line) || is.na(maxt.line) || is.null(maxt.line))) {
      stop("'maxt.line' must be a numeric vector of length 1, 'NA', or 'NULL'")
    }
    if (is.snum(maxt.line) && is.infinite(maxt.line)) {
      stop("'maxt.line' must be finite")
    }
  }
  
  if (!missing(xaxt.line)) {
    if (!(is.snum(xaxt.line) || is.na(xaxt.line) || is.null(xaxt.line))) {
      stop("'xaxt.line' must be a numeric vector of length 1, 'NA', or 'NULL'")
    }
    if (is.snum(xaxt.line) && is.infinite(xaxt.line)) {
      stop("'xaxt.line' must be finite")
    }
  }
  
  if (!missing(yaxt.line)) {
    if (!(is.snum(yaxt.line) || is.na(yaxt.line) || is.null(yaxt.line))) {
      stop("'yaxt.line' must be a numeric vector of length 1, 'NA', or 'NULL'")
    }
    if (is.snum(yaxt.line) && is.infinite(yaxt.line)) {
      stop("'yaxt.line' must be finite")
    }
  }
  
  if (!missing(laxt.line)) {
    if (!(is.snum(laxt.line) || is.na(laxt.line) || is.null(laxt.line))) {
      stop("'laxt.line' must be a numeric vector of length 1, 'NA', or 'NULL'")
    }
    if (is.snum(laxt.line) && is.infinite(laxt.line)) {
      stop("'laxt.line' must be finite")
    }
  }
  
  if (!missing(xlab.line)) {
    if (!(is.snum(xlab.line) || is.na(xlab.line) || is.null(xlab.line))) {
      stop("'xlab.line' must be a numeric vector of length 1, 'NA', or 'NULL'")
    }
    if (is.snum(xlab.line) && is.infinite(xlab.line)) {
      stop("'xlab.line' must be finite")
    }
  }
  
  if (!missing(ylab.line)) {
    if (!(is.snum(ylab.line) || is.na(ylab.line) || is.null(ylab.line))) {
      stop("'ylab.line' must be a numeric vector of length 1, 'NA', or 'NULL'")
    }
    if (is.snum(ylab.line) && is.infinite(ylab.line)) {
      stop("'ylab.line' must be finite")
    }
  }
  
  if (!missing(main.line)) {
    if (!(is.snum(main.line) || is.na(main.line) || is.null(main.line))) {
      stop("'main.line' must be a numeric vector of length 1, 'NA', or 'NULL'")
    }
    if (is.snum(main.line) && is.infinite(main.line)) {
      stop("'main.line' must be finite")
    }
  }
  
  if (!missing(legend.main.line)) {
    if (!(is.snum(legend.main.line) || is.na(legend.main.line) || 
          is.null(legend.main.line))) {
      stop("'legend.main.line' must be a numeric vector of length 1, 'NA', or 
           'NULL'")
    }
    if (is.snum(legend.main.line) && is.infinite(legend.main.line)) {
      stop("'legend.main.line' must be finite")
    }
  }
  
  if (!missing(xaxt.tck)) {
    if (!(is.snum(xaxt.tck) || is.na(xaxt.tck) || is.null(xaxt.tck))) {
      stop("'xaxt.tck' must be a numeric vector of length 1, or 'NA'")
    }
    if (is.snum(xaxt.tck) && is.infinite(xaxt.tck)) {
      stop("'xaxt.tck' must be finite")
    }
  }
  
  if (!missing(yaxt.tck)) {
    if (!(is.snum(yaxt.tck) || is.na(yaxt.tck) || is.null(yaxt.tck))) {
      stop("'yaxt.tck' must be a numeric vector of length 1, or 'NA'")
    }
    if (is.snum(yaxt.tck) && is.infinite(yaxt.tck)) {
      stop("'yaxt.tck' must be finite")
    }
  }
  
  if (!missing(laxt.tck)) {
    if (!(is.snum(laxt.tck) || is.na(laxt.tck) || is.null(laxt.tck))) {
      stop("'laxt.tck' must be a numeric vector of length 1, or 'NA'")
    }
    if (is.snum(laxt.tck) && is.infinite(laxt.tck)) {
      stop("'laxt.tck' must be finite")
    }
  }
  
  if (!(missing(plotLegend) || is.slog(plotLegend)))
    stop("'plotLegend' must be on of 'TRUE' or 'FALSE'")
  
  if (!missing(legend.position)) {
    if (!(is.snum(legend.position) || is.na(legend.position) || is.null(legend.position))) {
      stop("'legend.position' must be a numeric vector of length 1, 'NA', or 'NULL'")
    }
    if (is.snum(legend.position) && is.infinite(legend.position)) {
      stop("'legend.position' must be finite")
    }
  }
  
  if (!(missing(legend.main) || is.schar(legend.main) || is.na(legend.main) 
        || is.null(legend.main)))
    stop("'legend.main' must be a character vector of length 1")
  
  if (!(missing(symmetric) || is.slog(symmetric)))
    stop("'symmetric' must be one of 'TRUE' or 'FALSE'")
  
  if (!(missing(horizontal) || is.slog(horizontal)))
    stop("'horizontal' must be one of 'TRUE' or 'FALSE'")
  
  if (!missing(dataCols)) {
    if (!(is.na(dataCols) || is.null(dataCols) || is.character(dataCols))) {
      stop("'dataCols' must be a character vector")
    } else if (any(!areColors(dataCols))) {
      stop("invalid colors found in 'dataCols':",  
        paste(paste0('"', dataCols[!areColors(dataCols)], '"'), collapse=", "))
    }
  }
  
  if (!missing(corCols)) {
    if (!(is.na(corCols) || is.null(corCols) || is.character(corCols))) {
      stop("'corCols' must be a character vector")
    } else if (any(!areColors(corCols))) {
      stop("invalid colors found in 'corCols':",  
           paste(paste0('"', corCols[!areColors(corCols)], '"'), collapse=", "))
    }
  }
  
  if (!missing(netCols)) {
    if (!(is.na(netCols) || is.null(netCols) || is.character(netCols))) {
      stop("'netCols' must be a character vector")
    } else if (any(!areColors(netCols))) {
      stop("invalid colors found in 'netCols':",  
           paste(paste0('"', netCols[!areColors(netCols)], '"'), collapse=", "))
    }
  }
  
  if (!missing(degreeCol)) {
    if (!(is.na(degreeCol) || is.null(degreeCol) || is.schar(degreeCol))) {
      stop("'degreeCol' must be a character vector of length 1")
    } else if (!areColors(degreeCol)) {
      stop('invalid color, "', degreeCol, '" for', " 'degreeCol'", sep="")  
    }
  }
  
  if (!missing(contribCols)) {
    if (!(is.na(contribCols) || is.null(contribCols) || is.character(contribCols))) {
      stop("'contribCols' must be a character vector")
    } else if  (is.character(contribCols) && length(contribCols) %nin% 1:2) {
      stop("'contribCols' must be of length 1 or 2")
    } else if (any(!areColors(contribCols))) {
      stop("invalid colors found in 'contribCols':",  
           paste(paste0('"', contribCols[!areColors(contribCols)], '"'), collapse=", "))
    }
  }
  
  if (!missing(summaryCols)) {
    if (!(is.na(summaryCols) || is.null(summaryCols) || is.character(summaryCols))) {
      stop("'summaryCols' must be a character vector")
    } else if  (is.character(summaryCols) && length(summaryCols) %nin% 1:2) {
      stop("'summaryCols' must be of length 1 or 2")
    } else if (any(!areColors(summaryCols))) {
      stop("invalid colors found in 'summaryCols':",  
           paste(paste0('"', summaryCols[!areColors(summaryCols)], '"'), collapse=", "))
    }
  }

  if (!missing(naCol)) {
    if (!(is.na(naCol) || is.null(naCol) || is.schar(naCol))) {
      stop("'naCol' must be a character vector of length 1")
    } else if (!areColors(naCol)) {
      stop('invalid color, "', naCol, '" for', " 'naCol'", sep="")  
    }
  }
  
  if (!missing(dataRange)) {
    if (!(missing(dataRange) || is.na(dataRange) || is.null(dataRange) || 
          (is.numeric(dataRange) && length(dataRange) == 2))) {
      stop("'dataRange' must be a numeric vector of length 2")
    } else if (is.numeric(dataRange) && any(is.infinite(dataRange))) {
      stop("infinite values found in 'dataRange'")
    }
  }
  
  if (!missing(corRange)) {
    if (!(missing(corRange) || (is.numeric(corRange) && length(corRange) == 2))) {
      stop("'corRange' must be a numeric vector of length 2")
    } else if (is.numeric(corRange) && any(is.infinite(corRange))) {
      stop("infinite values found in 'corRange'")
    }
  }
  
  if (!missing(netRange)) {
    if (!(missing(netRange) || is.na(netRange) || is.null(netRange) || 
          (is.numeric(netRange) && length(netRange) == 2))) {
      stop("'netRange' must be a numeric vector of length 2 or 'NA'")
    } else if (is.numeric(netRange) && any(is.infinite(netRange))) {
      stop("infinite values found in 'netRange'")
    }
  }
  
  if (!(missing(dryRun) || is.slog(dryRun)))
    stop("'dryRun' must be one of 'TRUE' or 'FALSE'")
  
}