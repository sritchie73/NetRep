#' Fast wrapper functions for Rcpp functions
#' 
#' @template wrapper_desc
#' @template lowmem_inputs
#' @template no_sanity
#' @inheritParams coexp_params
#' @inheritParams adj_param
#' @inheritParams ge_param
#' @inheritParams ind_param
#' 
#' @name wrappers
NULL

#' @rdname wrappers
#' @return 
#'  \code{coexpStats:} a list containing the \emph{cor.coexp} and 
#'  \emph{mean.coexp} statistics for the specified module across datasets.
coexpStats <- function(
  discCoexp, discIndices, testCoexp, testIndices
) {
  # Attach the big.matrix objects if not attached yet
  disc.attached <- discCoexp@attached
  test.attached <- testCoexp@attached
  if (!disc.attached)
    discCoexp <- attach.bigMatrix(discCoexp)
  if (!test.attached)
    testCoexp <- attach.bigMatrix(testCoexp)
  
  res <- CoexpStats(
    discCoexp@matrix@address, discIndices, 
    testCoexp@matrix@address, testIndices
  )
  
  # detach big.matrix objects if they were detached to begin with
  if (!disc.attached)
    discCoexp <- detach.bigMatrix(discCoexp)
  if (!test.attached)
    testCoexp <- detach.bigMatrix(testCoexp)

  lapply(res, as.vector)
}

#' @rdname wrappers
#' @return 
#'  \code{netProps:} a list of network properties calculated from an adjacency
#'  matrix. These properties can either be scalers (summarising the whole
#'  network subset), or vectors (characterising some property for each node in
#'  the network subset).
adjProps <- function(adj, moduleIndices) {
  # Attach the big.matrix object if not attached yet
  is.attached <- adj@attached
  if (!is.attached)
    adj <- attach.bigMatrix(adj)
  
  res <- AdjProps(adj@matrix@address, moduleIndices)
  
  # detach big.matrix objects if they were detached to begin with
  if (!is.attached)
    adj <- detach.bigMatrix(adj)
  
  lapply(res, as.vector)
}

#' @rdname wrappers
#' @return 
#'   \code{dataProps:} a list of properties quantifying the relationship between
#'   a network subset and the underlying data the adjacency matrix was 
#'   calculated from. These properties can either be scalers (summarising the 
#'   whole network subset), or vectors (characterising some property for each 
#'   node in the network subset).
dataProps <- function(sge, moduleIndices) {
  # Attach the big.matrix object if not attached yet
  is.attached <- sge@attached
  if (!is.attached)
    sge <- attach.bigMatrix(sge)
  
  res <- DataProps(sge@matrix@address, moduleIndices)
  
  # detach big.matrix objects if they were detached to begin with
  if (!is.attached)
    sge <- detach.bigMatrix(sge)
  
  lapply(res, as.vector)
}

#' Get the properties of a network module
#' 
#' @template api_inputs
#' 
#' @export
networkProperties <- function(
  geneExpression=NULL, coexpression, adjacency, moduleAssignments, modules,
  discovery=1, test=1
) {
  #checkSets(geneExpression, coexpression, adjacency, moduleAssignments)

  # Unify data structures for the rest of the function
  if (!is.null(geneExpression) & !is.list(geneExpression))
    geneExpression <- list(geneExpression)
  if (!is.list(coexpression))
    coexpression <- list(coexpression)
  if (!is.list(adjacency))
    adjacency <- list(adjacency)
  if (!is.list(moduleAssignments))
    moduleAssignments <- list(moduleAssignments)
  
  if (is.null(moduleAssignments[[discovery]]))
    stop("no module assignments in the discovery dataset")
  
  # Temporarily create scaled gene expression set for the calculation of the
  # summary expression profile
  sge <- NULL
  if (!is.null(geneExpression)) {
    tryCatch({
      checkFinite(geneExpression[[test]]) 
    }, error = function(e) {
      stop("Non-finite values encountered for the test dataset gene expression")
    })
    obj.dir <- paste0(".temp-objects", getUUID())
    sge <- scaleBigMatrix(geneExpression[[test]], obj.dir)
    on.exit({
      unlink(obj.dir, recursive=TRUE)
    }, add=TRUE)
  }
  
  # Get the properties for each module of interest 
  res <- lapply(modules, function(mod) {
    # Get the row/column indices of the module in the dataset of interest 
    sub <- moduleAssignments[[discovery]][moduleAssignments[[discovery]] == mod]
    modInds <- match(names(sub), rownames(coexpression[[test]]))
    na.inds <- which(is.na(modInds))
    modInds <- na.omit(modInds)
    
    if (length(modInds) == 0) {
      stop(
        "none of the genes for module ", mod, 
        " are present in the test dataset"
      )
    }
    
    tryCatch({
      checkFinite(adjacency[[test]]) 
    }, error = function(e) {
      stop("Non-finite values encountered for the test dataset adjacency")
    })
    
    # Get the properties calculated from the gene expression
    geProps <- NULL
    if (!is.null(sge)) {
      geProps <- dataProps(sge, modInds)
      # rename for clarity
      names(geProps) <- c("summaryExpression", "moduleMembership", "propVarExpl")
      geProps[[2]] <- insert.nas(geProps[[2]], na.inds)
      names(geProps[[1]]) <- rownames(sge)
      names(geProps[[2]]) <- names(sub)
    }
    
    # Get the properties calculated from the adjacency.
    adjProps <- adjProps(adjacency[[test]], modInds)
    names(adjProps) <- c("connectivity", "density")
    adjProps[[1]] <- insert.nas(adjProps[[1]], na.inds)
    names(adjProps[[1]]) <- names(sub)
    
    c(geProps, adjProps)
  })
  if (length(res) == 1) {
    res <- res[[1]]
  } else {
    names(res) <- modules
  } 
  res
}
