#' Order nodes and modules within a network.
#' 
#' Order nodes in descending order of \emph{weighted degree} and order 
#' modules by the similarity of their summary vectors.
#' 
#' @inheritParams common_params
#' @inheritParams simplify_param
#' @inheritParams orderModules_param
#' 
#' @param na.rm logical; If \code{TRUE}, nodes and modules present in the 
#'   \code{discovery} dataset but missing from the test dataset are excluded. If
#'   \code{FALSE}, missing nodes and modules are put last in the ordering.
#' @param mean logical; if \code{TRUE}, node order will be calculated for each
#'   \code{discovery} dataset by averaging the weighted degree and pooling 
#'   \emph{module summary} vectors across the specified \code{test} datasets.
#'   If \code{FALSE}, the node order is calculated separately in each test 
#'   dataset.
#'      
#' @details
#'  \subsection{Input data structure:}{
#'   The \link[=modulePreservation]{preservation of network modules} in a second
#'   dataset is quantified by measuring the preservation of topological
#'   properties between the \emph{discovery} and \emph{test} datasets. These 
#'   properties are calculated not only from the interaction networks inferred
#'   in each dataset, but also from the data used to infer those networks (e.g.
#'   gene expression data) as well as the correlation structure between 
#'   variables/nodes. Thus, all functions in the \code{NetRep} package have the 
#'   following arguments: 
#'   \itemize{
#'     \item{\code{network}:}{
#'       a list of interaction networks, one for each dataset.
#'     }
#'     \item{\code{data}:}{
#'       a list of data matrices used to infer those networks, one for each 
#'       dataset.
#'     }
#'     \item{\code{correlation}:}{
#'      a list of matrices containing the pairwise correlation coefficients 
#'      between variables/nodes in each dataset.
#'     } 
#'     \item{\code{moduleAssignments}:}{
#'      a list of vectors, one for each \emph{discovery} dataset, containing 
#'      the module assignments for each node in that dataset.
#'     }
#'     \item{\code{modules}:}{
#'      a list of vectors, one for each \emph{discovery} dataset, containing
#'      the names of the modules from that dataset to analyse.  
#'     }
#'     \item{\code{discovery}:}{
#'       a vector indicating the names or indices of the previous arguments' 
#'       lists to use as the \emph{discovery} dataset(s) for the analyses.
#'     }
#'     \item{\code{test}:}{
#'       a list of vectors, one vector for each \emph{discovery} dataset, 
#'       containing the names or indices of the \code{network}, \code{data}, and 
#'       \code{correlation} argument lists to use as the \emph{test} dataset(s) 
#'       for the analysis of each \emph{discovery} dataset.
#'     }
#'   }
#'   
#'   The formatting of these arguments is not strict: each function will attempt
#'   to make sense of the user input. For example, if there is only one 
#'   \code{discovery} dataset, then input to the \code{moduleAssigments} and 
#'   \code{test} arguments may be vectors, rather than lists. If the 
#'   \code{nodeOrder} are being calculate within the \emph{discovery} or
#'   \emph{test} datasets, then the \code{discovery} and \code{test} arguments do
#'   not need to be specified, and the input matrices for the \code{network},
#'   \code{data}, and \code{correlation} arguments do not need to be wrapped in
#'   a list.
#' }
#' 
#' @references
#' \enumerate{
#'    \item{
#'      Langfelder, P., Mischel, P. S. & Horvath, S. \emph{When is hub gene 
#'      selection better than standard meta-analysis?} PLoS One \strong{8}, 
#'      e61505 (2013).
#'    }
#' }
#' 
#' @return
#'  A nested list structure. At the top level, the list has one element per 
#'  \code{'discovery'} dataset. Each of these elements is a list that has one
#'  element per \code{'test'} dataset analysed for that \code{'discovery'} 
#'  dataset. Each of these elements is a list that has one element per 
#'  \code{'modules'} specified, containing a vector of node names for the
#'  requested module. When \code{simplify = TRUE} then the simplest possible 
#'  structure will be returned. E.g. if the node ordering are requested for 
#'  module(s) in only one dataset, then a single vector of node labels will
#'  be returned. 
#'  
#' When \code{simplify = FALSE} then a nested list of datasets will always be 
#' returned, i.e. each element at the top level and second level correspond to 
#' a dataset, and each element at the third level will correspond to modules 
#' discovered in the dataset specified at the top level if module labels are 
#' provided in the corresponding \code{moduleAssignments} list element. E.g. 
#' \code{results[["Dataset1"]][["Dataset2"]][["module1"]]} will contain the 
#' order of nodes calculated in "Dataset2", where "module1" was indentified in
#' "Dataset1". Modules and datasets for which calculation of the node order 
#' have not been requested will contain \code{NULL}.
#'  
#' @seealso \code{\link{networkProperties}}
#'  
#' @examples
#' \dontrun{
#' # load in example data, correlation, and network matrices for a discovery
#' # and test dataset:
#' data("NetRep")
#' 
#' # Set up input lists for each input matrix type across datasets. The list
#' # elements can have any names, so long as they are consistent between the
#' # inputs.
#' data_list <- list(discovery=discovery_data, test=test_data)
#' correlation_list <- list(discovery=discovery_correlation, test=test_correlation)
#' network_list <- list(discovery=discovery_network, test=test_network)
#' labels_list <- list(discovery=module_labels)
#' 
#' # Sort nodes by module similarity and node degree
#' nodes <- nodeOrder(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list
#' )
#' }
#' 
#' @name nodeOrder
#' @export
nodeOrder <- function(
  data=NULL, correlation, network, moduleAssignments=NULL, modules=NULL, 
  backgroundLabel="0", discovery=NULL, test=NULL, na.rm=FALSE, 
  orderModules=TRUE, mean=FALSE, simplify=TRUE, verbose=TRUE
) {
  #-----------------------------------------------------------------------------
  # Input processing and sanity checking
  #-----------------------------------------------------------------------------
  vCat(verbose, 0, "Validating user input...")
  
  if (!is.logical(na.rm) || is.na(na.rm) || length(na.rm) > 1) {
    stop("'na.rm' must be either 'TRUE' or 'FALSE'")
  }
  
  if (!is.logical(orderModules) || is.na(orderModules) || length(orderModules) > 1) {
    stop("'orderModules' must be either 'TRUE' or 'FALSE'")
  }
  
  # We can save time and space by not checking, scaling, or calculating 
  # properties from the 'data' if we're not ordering by module. 
  if (!orderModules)
    data <- NULL
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose)
  
  # If 'orderModules' is TRUE we need to make sure that there is data in each 
  # test dataset if there is more than one 'module'
  with(finput, {
    if (orderModules) {
      for (di in discovery) {
        for (ti in test[[di]]) {
          if (is.null(data[[ti]]) && length(modules[[di]]) > 1) {
            stop("'data' must be provided for all 'test' datasets ",
                 "if 'orderModules' is TRUE")
          }   
        }
      }
    }
  })
  
  vCat(verbose, 0, "User input ok!")
  
  # Similarly when 'orderModules = TRUE' but there is only one module (per test)
  # Then we can save time by not calculating the data-based network properties.
  # We can't detect this without sanity checking the user input first though, so
  # unnecessary scaled datasets may have been created
  with(finput, {
    for (di in discovery) {
      if (length(modules[[di]]) == 1) {
        data[di] <- list(NULL)
      }
    }
  })
  
  # Calculate the network properties
  props <- with(finput, { 
    netPropsInternal(
      data, correlation, network, moduleAssignments, 
      modules, discovery, test, nDatasets, datasetNames, verbose
    ) 
  })
  
  res <- nodeOrderInternal(
    props, orderModules, simplify, verbose, na.rm, mean
  )
  
  # Simplify the output data structure where possible
  if (simplify) {
    res <- simplifyList(res, depth=3)
  }
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
  res
  return(res)
}

#' Internal use function for calculating node order
#' 
#' Used by plotting functions: assumes user input has been sanitised already
#' 
#' @param props network properties to calculate order from.
#' @param orderModules logical; are we ordering the modules by similarity?
#' @param simplify logical; are we simplifying the output?
#' @param verbose logical; is verbose printing turned on?
#' @param na.rm logical; remove missing nodes?
#' @param mean logical; should the node order be calculated by averaging across
#'  test datasets?
#' 
#' @return list structure of ordered nodes.
nodeOrderInternal <- function(
  props, orderModules, simplify, verbose, na.rm, mean
) {
  vCat(verbose, 0, "Ordering nodes...")
  
  # Average weighted degree and pool summary profiles across datasets if 'mean'
  # is 'TRUE'
  if (mean) {
    avgProps <- foreach(ii = seq_along(props)) %do% {
      # Skip datasets if they have not had modules discovered in them 
      if (is.null(props[[ii]][[1]])) {
        return(list(avgProps=NULL))
      }
      warned <- FALSE
      # For each module in the discovery dataset, loop through the test datasets
      # to combine the properties.
      modProps <- foreach(mi = names(props[[ii]][[1]])) %do% {
        # If the order of this module has not been requested, return NULL
        notRequested <- foreach(jj = seq_along(props[[ii]]), .combine=c) %do% { 
          is.null(props[[ii]][[jj]][[mi]]) 
        }
        if (all(notRequested)) {
          return(NULL)
        }
        
        if (sum(!notRequested) == 1 & !warned) {
          warning("'mean' is 'TRUE' where only one 'test' dataset specified ", 
                  "for discovery dataset ", '"', names(props)[ii], '"')
          warned <- TRUE # suppress printing out many warnings across modules
        }
        
        # Get the node names
        nodeNames <- names(props[[ii]][!notRequested][[1]][[mi]][["degree"]])
        
        # Otherwise calculate the mean weighted degree and concatenate the 
        # summary profiles
        degreeMat <- foreach(jj = seq_along(props[[ii]]), .combine = rbind) %do% {
          degree <- props[[ii]][[jj]][[mi]][["degree"]]
        }
        # Gracefully handles case where the module is only requested across 1 
        # dataset.
        if(is.null(dim(degreeMat))) {
          degreeMat <- matrix(degreeMat, nrow=1)
        }
        avgDegree <- colMeans(degreeMat, na.rm=TRUE)
        names(avgDegree) <- nodeNames
        
        pooledSummary <- foreach(jj = seq_along(props[[ii]]), .combine=c) %do% {
          props[[ii]][[jj]][[mi]][["summary"]]
        }
        return(list(degree = avgDegree, summary = pooledSummary))
      }
      names(modProps) <- names(props[[ii]][[1]])
      return(list(avgProps=modProps))
    }
    names(avgProps) <- names(props)
    props <- avgProps
  }
  
  for (ii in seq_along(props)) {
    for (jj in seq_along(props[[ii]])) {
      hasProps <- which(!sapply(props[[ii]][[jj]], is.null))
      if (length(hasProps) == 0) {
        next
      }
      # First get module order
      if (length(hasProps) > 1 && orderModules) {
        modProps <- props[[ii]][[jj]][hasProps]
        
        # First get the summary profiles for each module
        summaries <- matrix(0, ncol=length(hasProps), 
                            nrow=length(modProps[[1]][["summary"]]))
        colnames(summaries) <- names(modProps)
        rownames(summaries) <- names(modProps[[1]][["summary"]])
        for (mi in seq_along(modProps)) {
          summaries[,mi] <- modProps[[mi]][["summary"]]
        }
        
        # Identify modules where no nodes are present in the test dataset
        na.inds <- which(apply(summaries, 2, function(s) all(is.na(s))))
        na.mods <- colnames(summaries)[na.inds]
        if (length(na.mods) > 0) {
          summaries <- summaries[,-na.mods]
        }
        
        # Cluster modules that do have summary profiles
        clusteredMods <- hclust(as.dist(1-cor(summaries)))
        hasProps <- colnames(summaries)[clusteredMods$order]
        
        # Add the 'NA' modules
        hasProps  <- c(hasProps , na.mods)
      }
      
      # order nodes within each module
      modProps <- props[[ii]][[jj]][hasProps]
      for (kk in rev(seq_along(modProps))) {
        nodeDegree <- modProps[[kk]][["degree"]]
        sortedNodes <- sort(nodeDegree, decreasing=TRUE, na.last=TRUE)
        # Remove missing nodes and modules
        if (na.rm) {
          sortedNodes <- na.omit(sortedNodes)
          if (length(sortedNodes) == 0) {
            sortedNodes <- NULL
          }
        }
        modProps[[kk]] <- names(sortedNodes)
      }
      if (simplify) {
        modProps <- unlist(modProps)
        names(modProps) <- NULL
      }
      props[[ii]][[jj]] <- modProps
    }
  }
  # remove the extra list level in 'avgProps': this was added so that the
  # above code would work regardless of the 'mean' argument.
  if (mean) { 
    for (ii in rev(seq_along(props))) {
      if (!is.null(props[[ii]][["avgProps"]])) {
        props[[ii]] <- props[[ii]][["avgProps"]]
      } else {
        props[ii] <- list(NULL)
      }
    }
  }
  return(props)
}
