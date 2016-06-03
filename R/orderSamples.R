#' Order samples within a network.
#' 
#' Get the order of samples within a module based on the module summary vector.
#' 
#' @inheritParams common_params
#' @inheritParams simplify_param
#' 
#' @param na.rm logical; If \code{TRUE} variables present in the 
#'   \code{discovery} dataset but missing from the \code{test} dataset are 
#'   excluded. If \code{FALSE} missing variables are put last in the ordering.
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
#'   \code{sampleOrder} are being calculate within the \emph{discovery} or
#'   \emph{test} datasets, then the \code{discovery} and \code{test} arguments do
#'   not need to be specified, and the input matrices for the \code{network},
#'   \code{data}, and \code{correlation} arguments do not need to be wrapped in
#'   a list.
#' }
#' 
#' @return
#'  A nested list structure. At the top level, the list has one element per 
#'  \code{'discovery'} dataset. Each of these elements is a list that has one
#'  element per \code{'test'} dataset analysed for that \code{'discovery'} 
#'  dataset. Each of these elements is a list that has one element per 
#'  \code{'modules'} specified, containing a vector of node names for the
#'  requested module. When \code{simplify = TRUE} then the simplest possible 
#'  structure will be returned. E.g. if the sample ordering are requested for 
#'  in only one dataset, then a single vector of node labels will be returned. 
#'  
#' When \code{simplify = FALSE} then a nested list of datasets will always be 
#' returned, i.e. each element at the top level and second level correspond to 
#' a dataset, and each element at the third level will correspond to modules 
#' discovered in the dataset specified at the top level if module labels are 
#' provided in the corresponding \code{moduleAssignments} list element. E.g. 
#' \code{results[["Dataset1"]][["Dataset2"]][["module1"]]} will contain the 
#' order of samples calculated in "Dataset2", where "module1" was indentified
#' in "Dataset1". Modules and datasets for which calculation of the sample
#' order have not been requested will contain \code{NULL}.
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
#' network_list <- list(discovery=discovery_network, test=test_network)
#' data_list <- list(discovery=discovery_data, test=test_data)
#' correlation_list <- list(discovery=discovery_correlation, test=test_correlation)
#' labels_list <- list(discovery=module_labels)
#' 
#' # Sort nodes within module 1 in descending order by module summary
#' samples <- sampleOrder(
#'   network=network_list, data=data_list, correlation=correlation_list,
#'   moduleAssignments=labels_list, modules="1" 
#' )
#' }
#' 
#' @seealso \code{\link{networkProperties}}
#' 
#' @name sampleOrder
#' @export
sampleOrder <- function(
  network, data, correlation, moduleAssignments=NULL, modules=NULL, 
  backgroundLabel="0", discovery=NULL, test=NULL, na.rm=FALSE, 
  simplify=TRUE, verbose=TRUE
) {
  #----------------------------------------------------------------------------
  # Input processing and sanity checking
  #----------------------------------------------------------------------------
  vCat(verbose, 0, "Validating user input...")
  
  if (!is.logical(na.rm) || is.na(na.rm) || length(na.rm) > 1) {
    stop("'na.rm' must be either 'TRUE' or 'FALSE'")
  }
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose)
  
  discovery <- finput$discovery
  test <- finput$test
  data <- finput$data
  
  # We need to make sure that there is data in each 'test' dataset.
  for (di in discovery) {
    for (ti in test[[di]]) {
      if (is.null(data[[ti]])) {
        stop("'data' must be provided for all 'test' datasets")
      }   
    }
  }
  
  vCat(verbose, 0, "User input ok!")
  
  # Calculate the network properties.
  props <- with(finput, {
    netPropsInternal(network, data, correlation, moduleAssignments, 
                     modules, discovery, test, nDatasets, datasetNames, verbose)
  })

  res <- sampleOrderInternal(props, verbose, na.rm)

  # Simplify the output data structure where possible
  if (simplify) {
    res <- simplifyList(res, depth=3)
  }
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
  res
  return(res)
}

#' Internal use function for calculating sample order
#' 
#' Used by plotting functions: assumes user input has been sanitised already
#' 
#' @param props network properties to calculate order from.
#' @param verbose logical; is verbose printing turned on?
#' @param na.rm logical; remove missing nodes?
#' 
#' @return list structure of ordered nodes.
sampleOrderInternal <- function(props, verbose, na.rm) {
  vCat(verbose, 0, "Ordering samples...")
  for (ii in seq_along(props)) {
    for (jj in seq_along(props[[ii]])) {
      for (kk in seq_along(props[[ii]][[jj]])) {
        if (!is.null(props[[ii]][[jj]][[kk]])) {
          summary <- props[[ii]][[jj]][[kk]][["summary"]]
          if (na.rm) {
            summary <- na.omit(summary)
          }
          summaryOrder <- order(summary, decreasing=TRUE, na.last=TRUE)
          
          if(is.null(names(summary))) {
            props[[ii]][[jj]][[kk]] <- summaryOrder
          } else {
            props[[ii]][[jj]][[kk]] <- names(summary)[summaryOrder]
          }
        }
      }
    }
  }
  return(props)
}
