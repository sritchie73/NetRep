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
#' @param tempdir temporary directory to save new objects in.
#' @param plotFunction logical; are we checking for a plot function?
#' @param orderNodesBy use input for the 'orderNodesBy' argument in the 
#'  plotting functions.
#' @param orderSamplesBy use input for the 'orderSamplesBy' argument in the 
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
  backgroundLabel, verbose, tempdir, plotFunction=FALSE, orderNodesBy=NA, 
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
  
  # Make sure they're 'bigMatrix' objects
  correlation <- lapply(correlation, dynamicMatLoad, tempdir, verbose)
  network <- lapply(network, dynamicMatLoad, tempdir, verbose)
  
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

  # Make sure they're 'bigMatrix' objects
  data <- lapply(data, dynamicMatLoad, tempdir, verbose)

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
        bglabel <- getUUID() 
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
    if (is.null(orderSamplesBy))
      orderSamplesBy <- test[[discovery]]
  }

  # ----------------------------------------------------------------------------
  # Check for data consistency
  # ----------------------------------------------------------------------------
  # Now that we have made sure the input data can be sensibly accessed using the
  # provided 'discovery' and 'test' inputs, we need to make sure the data is 
  # consistent across the input lists.

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

  # Do the names need to match?
  if (matchByIndice) {
    iterator <- seq_along(network)
    if (all(hasNames) && (any(datNames != corNames) || any(corNames != netNames) || 
        any(modNames != netNames) || any(modLabNames != modNames))) {
      stop("mismatch in dataset names across input arguments when matching ",
           "'discovery' or 'test' by index")
    }
  } else {
    iterator <- names(network)
  }
  
  for (ii in iterator) {
    # Make sure the 'correlation' and 'network' matrices are square
    if (nrow(network[[ii]]) != ncol(network[[ii]])) {
      stop("'network' for dataset ", '"', ii, '"', " is not square")
    }
    if (nrow(correlation[[ii]]) != ncol(correlation[[ii]])) {
      stop("'correlation' for dataset ", '"', ii, '"', " is not square")
    }
    # And that they have the same dimensions
    if ((nrow(correlation[[ii]]) != nrow(network[[ii]])) ||
        (!is.null(data[[ii]]) && (ncol(data[[ii]]) != ncol(network[[ii]])))) {
      stop("'correlation', 'network', and 'data' have a different number of ",
           'nodes for dataset "', ii, '"')
    }
    
    # Make sure the 'correlation' and 'network' matrices are symmetric 
    if (any(rownames(network[[ii]]) != colnames(network[ii]))) {
      stop("mismatch between row and column names in 'network' for dataset ", 
           '"', ii, '"')
    }
    if (any(rownames(correlation[[ii]]) != colnames(correlation[ii]))) {
      stop("mismatch between row and column names in 'network' for dataset ",
           '"', ii, '"')
    }
    # Make sure the ordering of nodes is the same between 'correlation', 
    # 'network' and 'data'.
    if (any(colnames(network[[ii]]) != colnames(correlation[ii])) |
        (!is.null(data[[ii]]) && any(colnames(network[[ii]]) != colnames(data[[ii]])))) {
      stop("mismatch in node order between 'correlation' and 'network' for",
           'dataset "', ii, '"')
    }
    
    # Make sure the 'moduleAssignments' have the same nodes as the 'correlation'
    # etc.
    if (!is.null(moduleAssignments[[ii]])) {
      if (any(names(moduleAssignments[[ii]]) %nin% colnames(network[[ii]]))) {
        stop("module assigments are present for nodes that are not in the",
             " 'network' inferred from dataset ", '"', ii, '"')
      }
    }
    
    # Make sure all module labels are in 'moduleAssignments'
    if (any(modules[[ii]] %nin% moduleAssignments[[ii]])) {
      stop("some 'modules' specified for dataset ", '"', ii, '"', 
          " are not in 'moduleAssignments' for dataset ", '"', ii, '"')
    }
  }
  
  if (plotFunction) {
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
    
    # Check that samples from the 'orderSamplesBy' dataset are present in the 
    # 'test' dataset to be drawn and that data for that dataset is present
    if (!is.na(orderSamplesBy) && is.null(data[[orderSamplesBy]])) {
      stop("'data' not provided for 'orderSamplesBy' dataset") 
    }
    
    if (!is.na(orderSamplesBy) && 
        sum(rownames(data[[orderSamplesBy]]) %in% rownames(data[[ti]])) == 0) {
      stop("no samples in the dataset specified by 'orderSamplesBy' are in the", 
           " 'test' dataset to be drawn.")
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
  
  # Sanity check input data for values that will cause the calculation of 
  # network properties and statistics to hang. This can take a while, so 
  # we only want to check datasets that we're analysing (especially for plotting)
  vCat(verbose, 1, "Checking matrices for non-finite values...")
  
  # Construct an iterator that includes only the datasets we're analysing:
  if (is.character(discovery)) {
    iterator <- match(discovery, names(network))
  } else {
    iterator <- discovery
  }
  for (tv in test) {
    if (is.character(tv)) {
      iterator <- c(iterator, match(tv, names(network)))
    } else {
      iterator <- c(iterator, tv)
    }
  }
  if (plotFunction) {
    if (orderModules) {
      if (is.character(orderNodesBy)) {
        iterator <- c(iterator, match(orderNodesBy, names(network)))
      } else if (is.numeric(orderNodesBy)) {
        iterator <- c(iterator, orderNodesBy)
      }
    }
    if (is.character(orderSamplesBy)) {
      iterator <- c(iterator, match(orderSamplesBy, names(network)))
    } else if (is.numeric(orderSamplesBy)) {
      iterator <- c(iterator, orderSamplesBy)
    }
  }
  iterator <- unique(iterator)

  # Now check finite
  for (ii in iterator) {
    if (!is.null(data[[ii]])) 
      checkFinite(data[[ii]])
    checkFinite(correlation[[ii]])
    checkFinite(network[[ii]])
  }
  
  # Temporarily create scaled data set for the calculation of the
  # summary expression profile. Again, to save space and time, only do this for
  # datasets we're analysing
  scaledData <- rep(list(NULL), length(data))
  names(scaledData) <- names(data)
  for (ii in iterator) {
    if (!is.null(data[[ii]])) {
      scaledData[[ii]] <- scaleBigMatrix(data[[ii]], tempdir)
    }
  }
  
  if (!is.null(names(network))) {
    datasetNames <- names(network)
  } else {
    datasetNames <- paste0("Dataset", seq_len(nDatasets))
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

  return(list(
    data=data, correlation=correlation, network=network, discovery=discovery,
    test=test, moduleAssignments=moduleAssignments, modules=modules,
    nDatasets=nDatasets, datasetNames=datasetNames, scaledData=scaledData,
    orderNodesBy=orderNodesBy, orderSamplesBy=orderSamplesBy
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

#' Dynamically detect and load a bigMatrix object depending on input type
#' 
#' @param object user input object
#' @param tempdir temporary directory to save objects in
#' @param verbose logical; should progress be reported? 
#' @param ... additional arguments to pass to read.bigMatrix
#'   
#' @return
#'  A 'bigMatrix' object or error.
dynamicMatLoad <- function(object, tempdir, verbose, ...) {
  basename <- paste0("tmp", getUUID())
  if (is.null(object)) {
    return(NULL)
  } else if (is.bigMatrix(object)) {
    return(object)
  } else if (class(object) == "character") {
    if (!file.exists(object))
      stop("file", object, "does not exist")
    
    # Is this file a big.matrix descriptor?
    if (readLines(object, 1) == "new(\"big.matrix.descriptor\"") {
      return(load.bigMatrix(object))
    } else {
      vCat(
        verbose, 1,
        "Creating new 'bigMatrix' in a temporary directory for file ", object,
        ". This could take a while."
      )
      backingfile <- file.path(tempdir, basename)
      return(read.bigMatrix(file=object, backingfile=backingfile, ...))
    }
    
  } else if (class(object) == "matrix") {
    vCat(
      verbose, 1,
      "Matrix encountered. Creating new 'bigMatrix' in a temporary directory.",
      " This could take a while."
    )
    backingfile <- file.path(tempdir, basename)
    return(as.bigMatrix(object, backingfile=backingfile, ...))
  } 
  stop("unable to load object of type ", class(object), " as a bigMatrix!")
}

#' Validate plot function arguments
#' 
#' Simple typechecking for the extensive plot arguments
#' 
#' @param orderModules user input for the corresponding argument in the plot functions.
#' @param plotNodesNames user input for the corresponding argument in the plot functions.
#' @param plotSampleNames user input for the corresponding argument in the plot functions.
#' @param plotModuleNames user input for the corresponding argument in the plot functions.
#' @param main user input for the corresponding argument in the plot functions.
#' @param drawBorders user input for the corresponding argument in the plot functions.
#' @param border.width user input for the corresponding argument in the plot functions.
#' @param gaxt.line user input for the corresponding argument in the plot functions.
#' @param saxt.line user input for the corresponding argument in the plot functions.
#' @param maxt.line user input for the corresponding argument in the plot functions.
#' @param legend.tick.size user input for the corresponding argument in the plot functions.
#' @param laxt.line user input for the corresponding argument in the plot functions.
#' 
checkPlotArgs <- function(
  orderModules, plotNodeNames, plotSampleNames, plotModuleNames, main, 
  drawBorders, border.width, gaxt.line, saxt.line, maxt.line, legend.tick.size, 
  laxt.line
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
  
  if (!is.slog(orderModules))
    stop("'orderModules' must be one of 'TRUE' or 'FALSE'")
  
  if (!is.slog(plotNodeNames)) 
    stop("'plotNodeNames' must be one of 'TRUE' or 'FALSE'")
  
  if (!is.slog(plotSampleNames)) 
    stop("'plotNodeNames' must be one of 'TRUE' or 'FALSE'")
  
  if (!(is.null(plotModuleNames) || is.slog(plotModuleNames)))
    stop("'plotModuleNames' must be one of 'TRUE', 'FALSE', or 'NULL'")
  
  if (!(is.null(main) || is.schar(main) || is.na(main)))
    stop("'main' must be 'NULL' or a character vector of length 1")
  
  if (!(is.slog(drawBorders)))
    stop("'drawBorders' must be one of 'TRUE' or 'FALSE'")
  
  if (!is.snum(border.width))
    stop("'border.width' must be a numeric vector of length 1")
  if (border.width < 0)
    stop("'border.width' must be greater than 0")
  
  if (!(is.snum(gaxt.line) || is.na(gaxt.line)))
    stop("'gaxt.line' must be a numeric vector of length 1 or 'NA'")
  
  if (!(is.snum(saxt.line) || is.na(saxt.line)))
    stop("'saxt.line' must be a numeric vector of length 1 or 'NA'")
  
  if (!(is.null(maxt.line) || is.snum(maxt.line) || is.na(maxt.line)))
    stop("'maxt.line' must be a numeric vector of length 1, 'NA', or 'NULL'")
  
  if (!is.snum(laxt.line) || is.na(laxt.line))
    stop("'laxt.line' must be a numeric vector of length 1 or 'NA'")
  
  if (!(is.snum(legend.tick.size) || is.na(legend.tick.size)))
    stop("'legend.tick.size' must be a numeric vector of length 1 or 'NA'")
  
  
}