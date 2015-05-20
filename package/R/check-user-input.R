#' Check input data for consistency
#' 
#' Make sure all objects are the right class, have the same dimensions, and
#' are ordered the same way across datasets.
#' 
#' @param geneExpression a list of \code{bigMatrix} objects
#' @param coexpression a list of \code{bigMatrix} objects
#' @param adjacency a list of \code{bigMatrix} objects
#' @param moduleAssignments a list of either NULL, or vectors assigning genes to
#'  modules in the \code{discovery} datasets.
#' @param discovery name or index denoting which dataset the module of
#'  interest was discovered in. See details.
#' @param test name or index denoting which dataset to apply the function to.
#'  See details.
#'  
checkSets <- function(
  geneExpression, coexpression, adjacency, moduleAssignments, discovery, test
) {
  if (class(discovery) %nin% c("character", "numeric", "integer"))
    stop("'discovery' must be a vector of dataset names or indices")
  if (class(test) %nin% c("character", "numeric", "integer"))
    stop("'test' must be a vector of dataset names or indices")
  
  classes <- c(
    class(coexpression),
    class(adjacency),
    class(moduleAssignments)
  )
  if (!is.null(geneExpression)) {
    classes <- c(class(geneExpression), classes)
  }
  
  # Make sure the objects are lists, lists of bigMatrix objects, and have the
  # same length.
  if (!all(classes == "list")) {
    stop("Expecting lists of 'bigMatrix' objects")
  }
  classes2 <- c(
    sapply(coexpression, class),
    sapply(adjacency, class)
  )
  if (!is.null(geneExpression)) {
    classes2 <- c(sapply(geneExpression, class), classes2)      
  }
  if(!all(classes2 == "bigMatrix")) {
    stop(
      "expecting 'bigMatrix' objects, or a list of 'bigMatrix' objects for ",
      "the 'geneExpression', coexpression', and 'adjacency' arguments"
    )
  }
  nDatasets <- c(length(coexpression), length(adjacency))
  if (!is.null(geneExpression)) {
    nDatasets <- c(length(geneExpression), nDatasets)
  }
  if (length(unique(nDatasets)) > 1) {
    stop(
      "different number of datasets encountered across the",
      " 'geneExpression', coexpression', and 'adjacency' arguments"
    )
  }
  
  # If we're subsetting by dataset index, make sure the ordering of datasets
  # matches.
  if (class(discovery) != "character" || class(test) != "character") {
    if (
      any(names(coexpression) != names(adjacency)) &&
      any(names(coexpression) != names(moduleAssignments)) &&
      (!is.null(geneExpression) && any(names(geneExpression) != names(coexpression)))
    ) {
      stop(
        "'discovery' or 'test' provided as a vector of indices, but ",
        "ordering of dataset is not the same between the ", 
        ifelse(is.null(geneExpression), "", "'geneExpression'"),
        " 'coexpression', 'adjacency', and 'moduleAssignments'"
      )
    }
  }
  
  if (!is.null(names(coexpression))) {
    datasets <- names(coexpression)
  } else {
    datasets <- seq_along(coexpression)
  }
  for (dd in datasets) {
    if (nrow(coexpression[[dd]]) != ncol(coexpression[[dd]]))
      stop("'coexpression' of dataset ", dd, " is not square")
    if (nrow(adjacency[[dd]]) != ncol(adjacency[[dd]]))
      stop("'adjacency' of dataset ", dd, " is not square")
    
    ncols <- c(ncol(coexpression[[dd]]), ncol(adjacency[[dd]]))
    if (!is.null(geneExpression)) {
      ncols <- c(ncol(geneExpression[[dd]]), ncols)
    }
    if (!is.null(moduleAssignments[[dd]])) {
      ncols <- c(ncols, length(moduleAssignments[[dd]]))
    }

    if (length(unique(ncols)) > 1) {
      stop(
        "Different number of genes encountered across arguments for dataset ",
        dd
      )
    }
    
    if (any(rownames(coexpression[[dd]]) != colnames(coexpression[[dd]]))) {
      stop(
        "'coexpression' of dataset ", dd, " has differing row and column",
        " ordering"
      )
    }
    if (any(rownames(adjacency[[dd]]) != colnames(adjacency[[dd]]))) {
      stop(
        "'adjacency' of dataset ", dd, " has differing row and column",
        " ordering"
      )
    }
    if (any(rownames(adjacency[[dd]]) != rownames(coexpression[[dd]]))) {
      stop(
        "'adjacency' and 'coexpression' of dataset ", dd, " have differing",
        " gene ordering"
      )
    }
    if (!is.null(moduleAssignments[[dd]])) {
      if (any(rownames(adjacency[[dd]]) != names(moduleAssignments[[dd]]))) {
        stop(
          "'moduleAssignments' has differing gene ordering to the ",
          "'coexpression' and 'adjacency' for dataset ", dd
        )
      }
    }
    for (di in discovery) {
      if (is.null(moduleAssignments[[di]])) {
        stop(
          "No 'moduleAssignments' present for discovery dataset ", di  
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
#' @param nDiscGenes number of genes in the discovery dataset. Only required if
#'  moduleAssignments is missing.
#' @param discGeneNames names of the genes in the discovery dataset. Only 
#'  required if moduleAssignments is missing.
#' 
#' @return
#'  List structure for moduleAssignments for internal calculations.
#'
formatModuleAssignments <- function(
  moduleAssignments, discovery, nDatasets, datasetNames, 
  nDiscGenes, discGeneNames
) {
  if (class(discovery) %nin% c("character", "numeric", "integer"))
    stop("'discovery' must be a vector of dataset names or indices")
  
  # Handle case where moduleAssignments is not provided
  if (missing(moduleAssignments)) {
    moduleAssignments <- rep(list(NULL), nDatasets)
    names(moduleAssignments) <- datasetNames
    moduleAssignments[[discovery]] <- rep("1", nDiscGenes)
    names(moduleAssignments[[discovery]]) <- discGeneNames
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

#' Format includeModules list
#' 
#' @param includeModules input provided by the user
#' @param discovery vector of discovery networks provided by the user
#' @param nDatasets number of datasets
#' @param datasetNames names of the datasets
#' 
#' @return
#'  List structure to match the rest of the input data
#'
formatInclude <- function(
  includeModules, discovery, nDatasets, datasetNames
) {
  if (missing(includeModules)) {
    includeModules <- rep(list(NULL), length(discovery))
    names(includeModules) <- discovery
  }
  
  if (class(discovery) %nin% c("character", "numeric", "integer"))
    stop("'discovery' must be a vector of dataset names or indices")
  
  msg <- paste0(
    "mismatch between number of 'discovery' specified and ",
    "'includeModules' provided"
  )
  if (!is.list(includeModules)) {
    if (length(discovery) > 1) {
      stop(msg)
    } else {
      tmp <- rep(list(NULL), nDatasets)
      names(tmp) <- datasetNames
      tmp[[discovery]] <- includeModules
      includeModules <- tmp
    }
  } else {
    if (
      !is.null(names(includeModules)) && 
        !all(names(includeModules) %in% datasetNames)
    ) {
      stop("unable to match all names of 'includeModules' to dataset names")
    }
    if (length(includeModules) < nDatasets) {
      tmp <- rep(list(NULL), nDatasets)
      names(tmp) <- datasetNames
      if (length(includeModules) == length(discovery)) {
        if (is.null(names(includeModules))) {
          tmp[discovery] <- includeModules
        } else {
          tmp[names(includeModules)] <- includeModules
        }
      } else {
        if (is.null(names(includeModules))) {
          stop(msg)
        } else {
          tmp[names(includeModules)] <- includeModules
        }
      }
      includeModules <- tmp
    }
  }
  includeModules
}

#' Format excludeModules list
#' 
#' @param excludeModules input provided by the user
#' @param discovery vector of discovery networks provided by the user
#' @param nDatasets number of datasets
#' @param datasetNames names of the datasets
#' 
#' @return
#'  List structure to match the rest of the input data
#'
formatExclude <- function(
  excludeModules, discovery, nDatasets, datasetNames
) {
  if (missing(excludeModules)) {
    excludeModules <- rep(list(NULL), length(discovery))
    names(excludeModules) <- discovery
  }
  
  if (class(discovery) %nin% c("character", "numeric", "integer"))
    stop("'discovery' must be a vector of dataset names or indices")
  
  msg <- paste0(
    "mismatch between number of 'discovery' specified and ",
    "'excludeModules' provided"
  )
  if (!is.list(excludeModules)) {
    if (length(discovery) > 1) {
      stop(msg)
    } else {
      tmp <- rep(list(NULL), nDatasets)
      names(tmp) <- datasetNames
      tmp[[discovery]] <- excludeModules
      excludeModules <- tmp
    }
  } else {
    if (
      !is.null(names(excludeModules)) && 
      !all(names(excludeModules) %in% datasetNames)
    ) {
      stop("unable to match all names of 'excludeModules' to dataset names")
    }
    if (length(excludeModules) < nDatasets) {
      tmp <- rep(list(NULL), nDatasets)
      names(tmp) <- datasetNames
      if (length(excludeModules) == length(discovery)) {
        if (is.null(names(excludeModules))) {
          tmp[discovery] <- excludeModules
        } else {
          tmp[names(excludeModules)] <- excludeModules
        }
      } else {
        if (is.null(names(excludeModules))) {
          stop(msg)
        } else {
          tmp[names(excludeModules)] <- excludeModules
        }
      }
      excludeModules <- tmp
    }
  }
  excludeModules
}


#' Dynamically detect and load a bigMatrix object depending on input type
#' 
#' @param object user input object
#' @param ... additional arguments to pass to read.bigMatrix or as.bigMatrix,
#'   usually a temporary directory for the backingpath
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
      backingname <- basename(object)
      backingpath <- gsub(backingname, "", object)
      backingname <- gsub(".desc", "", backingpath)
      return(load.bigMatrix(backingname, backingpath))
    } else {
      vCat(
        TRUE, 0,
        "Creating new 'bigMatrix' in a temporary directory for file ", object,
        ". This could take a while."
      )
      return(read.bigMatrix(file=object, backingname=basename, ...))
    }
    
  } else if (class(object) == "matrix") {
    vCat(
      TRUE, 0,
      "Matrix encountered. Creating new 'bigMatrix' in a temporary directory.",
      " This could take a while."
    )
    return(as.bigMatrix(object, backingname=basename, ...))
  } 
  stop("unable to load object of type ", class(object), " as a bigMatrix!")
}
