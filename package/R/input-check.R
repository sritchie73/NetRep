#' Check input data for consistency
#' 
#' Make sure all objects are the right class, have the same dimensions, and
#' are ordered the same way across datasets.
#' 
#' #' @param geneExpression a \code{\link{bigMatrix}} object, or a list of 
#'  \code{bigMatrix} objects, see details.
#' @param coexpression a \code{\link{bigMatrix}} object, or a list of 
#'  \code{bigMatrix} objects, see details.
#' @param adjacency a \code{\link{bigMatrix}} object, or a list of 
#'  \code{bigMatrix} objects, see details.
#' @param moduleAssignments a vector assigning genes to modules, or a list of 
#'  such vectors. See details.
#'  
checkSets <- function(
  geneExpression, coexpression, adjacency, moduleAssignments
) {
  if (!is.list(moduleAssignments) && is.list(adjacency))
    moduleAssignments <- list(moduleAssignments)
  
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
  if (any(classes == "list") && !all(classes == "list")) {
    stop("Input error: some arguments contain multiple datasets, but not all")
  }
  if (all(classes == "list")) {
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
  } else if (!all(head(classes, -1) == "bigMatrix")) {
    stop(
      "expecting 'bigMatrix' objects, or a list of 'bigMatrix' objects for ",
      "the 'geneExpression', coexpression', and 'adjacency' arguments"
    )
  }
  
  # Unify data structures for the rest of the function
  if (!is.null(geneExpression) & !is.list(geneExpression))
    geneExpression <- list(geneExpression)
  if (!is.list(coexpression))
    coexpression <- list(coexpression)
  if (!is.list(adjacency))
    adjacency <- list(adjacency)
  if (!is.list(moduleAssignments))
    moduleAssignments <- list(moduleAssignments)
  
  if (length(moduleAssignments) < length(coexpression)) {
    tmp <- rep(list(NULL), length(coexpression))
    tmp[1:length(moduleAssignments)] <- moduleAssignments
    moduleAssignments <- tmp
  } else if (length(moduleAssignments) < length(coexpression)) {
    stop("Too many datasets in the 'moduleAssignments' argument")
  }
  
  for (dd in seq_along(coexpression)) {
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
  }
} 