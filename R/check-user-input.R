#' Check input data for consistency
#' 
#' Make sure all objects are the right class, have the same dimensions, and
#' are ordered the same way across datasets.
#' 
#' @inheritParams common_params
#' @param discovery name or index denoting which dataset the module of interest
#'   was discovered in.
#' @param test name or index denoting which dataset to apply the function to.
#'  
checkSets <- function(
  data, correlation, network, moduleAssignments, discovery, test
) {
  if (class(discovery) %nin% c("character", "numeric", "integer"))
    stop("'discovery' must be a vector of dataset names or indices")
  if (class(test) %nin% c("character", "numeric", "integer"))
    stop("'test' must be a vector of dataset names or indices")
  
  classes <- c(
    class(correlation),
    class(network),
    class(moduleAssignments)
  )
  if (!is.null(data)) {
    classes <- c(class(data), classes)
  }
  
  # Make sure the objects are lists, lists of bigMatrix objects, and have the
  # same length.
  if (!all(classes == "list")) {
    stop("Expecting lists of 'bigMatrix' objects.")
  }
  classes2 <- c(
    sapply(correlation, class),
    sapply(network, class)
  )
  if (!is.null(data)) {
    geClasses <- sapply(data, class)
    geClasses <- geClasses[geClasses != "NULL"]
    classes2 <- c(geClasses, classes2)      
  }
  if(!all(classes2 == "bigMatrix")) {
    stop(
      "expecting 'bigMatrix' objects, or a list of 'bigMatrix' objects for ",
      "the 'data', correlation', and 'network' arguments."
    )
  }
  nDatasets <- c(length(correlation), length(network))
  if (!is.null(data)) {
    nDatasets <- c(length(data), nDatasets)
  }
  if (length(unique(nDatasets)) > 1) {
    stop(
      "different number of datasets encountered across the",
      " 'data', correlation', and 'network' arguments."
    )
  }
  
  # If we're subsetting by dataset index, make sure the ordering of datasets
  # matches.
  if (class(discovery) != "character" || class(test) != "character") {
    if (
      any(names(correlation) != names(network)) &&
      any(names(correlation) != names(moduleAssignments)) &&
      (!is.null(data) && any(names(data) != names(correlation)))
    ) {
      stop(
        "'discovery' or 'test' provided as a vector of indices, but ",
        "ordering of dataset is not the same between the ", 
        ifelse(is.null(data), "", "'data'"),
        " 'correlation', 'network', and 'moduleAssignments' lists."
      )
    }
  }
  
  if (!is.null(names(correlation))) {
    datasets <- names(correlation)
  } else {
    datasets <- seq_along(correlation)
  }
  
  # Check dimensions of data is ok
  for (dd in datasets) {
    if (nrow(correlation[[dd]]) != ncol(correlation[[dd]]))
      stop("'correlation' matrix for dataset ", dd, " is not square.")
    if (nrow(network[[dd]]) != ncol(network[[dd]]))
      stop("'network' matrix for dataset ", dd, " is not square.")
    
    ncols <- c(ncol(correlation[[dd]]), ncol(network[[dd]]))
    if (!is.null(data)) {
      ncols <- c(ncol(data[[dd]]), ncols)
    }
    if (!is.null(moduleAssignments[[dd]])) {
      ncols <- c(ncols, length(moduleAssignments[[dd]]))
    }

    if (length(unique(ncols)) > 1) {
      stop(
        "Input matrices for dataset ", dd, " do not contain the same number of",
        " variables (columns)."
      )
    }
    
    # Check variable name order is the same within and across input matrices
    if (any(rownames(correlation[[dd]]) != colnames(correlation[[dd]]))) {
      stop(
        "'correlation' matrix for dataset ", dd, " does not have consistent", 
        "row and column name ordering."
      )
    }
    if (any(rownames(network[[dd]]) != colnames(network[[dd]]))) {
      stop(
        "'network' matrix for dataset ", dd, " does not have consistent", 
        "row and column name ordering."
      )
    }
    if (any(rownames(network[[dd]]) != rownames(correlation[[dd]]))) {
      stop(
        "'network' and 'correlation' matrices for dataset ", dd, " do not have",
        " consistent ordering of variables (column names)."
      )
    }
    if (!is.null(data) & !is.null(data[[dd]])) {
      if (any(colnames(data[[dd]]) != colnames(correlation[dd]))) {
        stop(
          "'data' and 'correlation' matrices for dataset ", dd, " do not have",
          " consistent ordering of variables (column names)."
        )
      }
    }
    if (!is.null(moduleAssignments[[dd]])) {
      if (any(rownames(network[[dd]]) != names(moduleAssignments[[dd]]))) {
        stop(
          "'moduleAssignments' has differing variable ordering to the ",
          "'correlation' and 'network' matrices for dataset ", dd, "."
        )
      }
    }
    for (di in discovery) {
      if (is.null(moduleAssignments[[di]])) {
        stop(
          "No 'moduleAssignments' present for discovery dataset ", di, "." 
        )
      }
    }
  }
} 

#' Format moduleAssignments list
#' 
#' If module discovery has not been performed for all datasets, it may be
#' easier for the user to provide a simplified list structure. This function
#' unifies it to the same structure as the other datasets.
#' 
#' @param moduleAssignments input provided by the user
#' @param discovery vector of discovery networks provided by the user
#' @param nDatasets number of datasets
#' @param datasetNames names of the datasets
#' @param nDiscVars number of variables in the discovery dataset. 
#'  Only required if \code{moduleAssignments} is missing.
#' @param discVarNames varaibles in the discovery dataset. Only 
#'  required if \code{moduleAssignments} is missing.
#' 
#' @return
#'  List structure for moduleAssignments for internal calculations.
#'
formatModuleAssignments <- function(
  moduleAssignments, discovery, nDatasets, datasetNames, 
  nDiscVars, discVarNames
) {
  if (class(discovery) %nin% c("character", "numeric", "integer"))
    stop("'discovery' must be a vector of dataset names or indices")
  
  # Handle case where moduleAssignments is not provided
  if (missing(moduleAssignments)) {
    moduleAssignments <- rep(list(NULL), nDatasets)
    names(moduleAssignments) <- datasetNames
    moduleAssignments[[discovery]] <- rep("1", nDiscVars)
    names(moduleAssignments[[discovery]]) <- discVarNames
  }
  
  # If we're subsetting by dataset index, make sure the ordering of datasets
  # matches.
  if (
    !is.null(moduleAssignments) &&
    (length(discovery) > nDatasets || 
    (class(discovery) == "character" && any(discovery %nin% datasetNames)) ||
    (length(discovery) > 1 && (
      !is.list(moduleAssignments) ||
      length(moduleAssignments) < length(discovery)
    )))
  ) {
    stop(
      "mismatch between number of 'discovery' specified and ",
      "'moduleAssignments' provided"
    )
  }
  
  if (!is.list(moduleAssignments))
    moduleAssignments <- list(moduleAssignments)
  
  # Now structure the moduleAssignments list sensibly.
  if (length(moduleAssignments) < nDatasets) {
    tmp <- rep(list(NULL), nDatasets)
    names(tmp) <- datasetNames
    tmp[discovery] <- moduleAssignments
    moduleAssignments <- tmp
  }
  moduleAssignments
}

#' Format the list of modules to include.
#' 
#' @param include input provided by the user
#' @param discovery vector of discovery networks provided by the user
#' @param nDatasets number of datasets
#' @param datasetNames names of the datasets
#' 
#' @return
#'  List structure to match the rest of the input data
#'
formatInclude <- function(
  include, discovery, nDatasets, datasetNames
) {
  if (missing(include)) {
    include <- rep(list(NULL), length(discovery))
    names(include) <- discovery
  }
  
  if (class(discovery) %nin% c("character", "numeric", "integer"))
    stop("'discovery' must be a vector of dataset names or indices")
  
  msg <- paste0(
    "mismatch between number of 'discovery' datasets specified and ",
    "number of 'include' vectors provided"
  )
  if (!is.list(include)) {
    # If not a list and there's only one discovery dataset, we can wrap in a 
    # list with the appropriate name
    if (length(discovery) > 1) {
      stop(msg)
    } else {
      tmp <- rep(list(NULL), nDatasets)
      names(tmp) <- datasetNames
      tmp[[discovery]] <- include
      include <- tmp
    }
  } else {
    if (!is.null(names(include)) && !all(names(include) %in% datasetNames)) {
      stop(
        "Unable to match list element names of 'include' to dataset names")
    }
    if (length(include) < nDatasets) {
      tmp <- rep(list(NULL), nDatasets)
      names(tmp) <- datasetNames
      if (length(include) == length(discovery)) {
        if (is.null(names(include))) {
          tmp[discovery] <- include
        } else {
          tmp[names(include)] <- include
        }
      } else {
        if (is.null(names(include))) {
          stop(msg)
        } else {
          tmp[names(include)] <- include
        }
      }
      include <- tmp
    }
  }
  include
}

#' Format exclude list
#' 
#' @param exclude input provided by the user
#' @param discovery vector of discovery networks provided by the user
#' @param nDatasets number of datasets
#' @param datasetNames names of the datasets
#' 
#' @return
#'  List structure to match the rest of the input data
#'
formatExclude <- function(
  exclude, discovery, nDatasets, datasetNames
) {
  if (missing(exclude)) {
    exclude <- rep(list(NULL), length(discovery))
    names(exclude) <- discovery
  }
  
  if (class(discovery) %nin% c("character", "numeric", "integer"))
    stop("'discovery' must be a vector of dataset names or indices")
  
  msg <- paste0(
    "mismatch between number of 'discovery' specified and ",
    "'exclude' provided"
  )
  if (!is.list(exclude)) {
    if (length(discovery) > 1) {
      stop(msg)
    } else {
      tmp <- rep(list(NULL), nDatasets)
      names(tmp) <- datasetNames
      tmp[[discovery]] <- exclude
      exclude <- tmp
    }
  } else {
    if (
      !is.null(names(exclude)) && 
      !all(names(exclude) %in% datasetNames)
    ) {
      stop("unable to match all names of 'exclude' to dataset names")
    }
    if (length(exclude) < nDatasets) {
      tmp <- rep(list(NULL), nDatasets)
      names(tmp) <- datasetNames
      if (length(exclude) == length(discovery)) {
        if (is.null(names(exclude))) {
          tmp[discovery] <- exclude
        } else {
          tmp[names(exclude)] <- exclude
        }
      } else {
        if (is.null(names(exclude))) {
          stop(msg)
        } else {
          tmp[names(exclude)] <- exclude
        }
      }
      exclude <- tmp
    }
  }
  exclude
}

#' Format the optional 'data' argument list
#' 
#' Allows the user to set \code{data} to \code{NULL} without the data
#' causing cascading input-check errors.
#' 
#' @param data user provided input to 'data' argument of package functions.
#' @param nDatasets total number of datasets
#' @param datasetNames names of the datasets
#' 
#' @return
#'  A formatted list of gene expression data
formatDataList <- function(data, nDatasets, datasetNames) {
  if (is.null(data) || identical(data, list(NULL))) {
    res <- rep(list(NULL), nDatasets)
    names(res) <- datasetNames
    res
  } else {
    data
  }
}

#' Dynamically detect and load a bigMatrix object depending on input type
#' 
#' @param object user input object
#' @param ... additional arguments to pass to read.bigMatrix
#'   
#' @return
#'  A 'bigMatrix' object or error.
dynamicMatLoad <- function(object, ...) {
  basename <- paste0("tmp", getUUID())
  if (is.null(object)) {
    return(NULL)
  } else if (is.list(object)) {
    ret <- lapply(object, dynamicMatLoad, ...)
    names(ret) <- names(object)
    return(ret)
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
        TRUE, 0,
        "Creating new 'bigMatrix' in a temporary directory for file ", object,
        ". This could take a while."
      )
      backingfile <- file.path(tempdir(), basename)
      return(read.bigMatrix(file=object, backingfile=backingfile, ...))
    }
    
  } else if (class(object) == "matrix") {
    vCat(
      TRUE, 0,
      "Matrix encountered. Creating new 'bigMatrix' in a temporary directory.",
      " This could take a while."
    )
    backingfile <- file.path(tempdir(), basename)
    return(as.bigMatrix(object, backingfile=backingfile, ...))
  } 
  stop("unable to load object of type ", class(object), " as a bigMatrix!")
}
