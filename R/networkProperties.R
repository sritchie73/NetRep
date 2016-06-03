#' Calculate the topological properties for a network module
#' 
#' Calculates the network properties used to assess module preservation for one
#' or more modules in a user specified dataset.
#' 
#' @inheritParams common_params
#' @inheritParams simplify_param
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
#'   \code{networkProperties} are being calculate within the \emph{discovery} or
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
#'  \code{'modules'} specified. Each of these is a list containing the following  
#'  objects:
#'  \itemize{
#'    \item{\code{'degree'}:}{
#'      The weighted within-module degree: the sum of edge weights for each 
#'      node in the module.
#'    }
#'    \item{\code{'avgWeight'}:}{
#'      The average edge weight within the module.
#'    }
#'  }
#'  If the \code{'data'} used to infer the \code{'test'} network is provided  
#'  then the following are also returned:
#'  \itemize{
#'    \item{\code{'summary'}:}{
#'      A vector summarising the module across each sample. This is calculated 
#'      as the first eigenvector of the module from a principal component 
#'      analysis.
#'    }
#'    \item{\code{'contribution'}:}{
#'      The \emph{node contribution}: the similarity between each node and the
#'      \emph{module summary profile} (\code{'summary'}).
#'    }
#'    \item{\code{'coherence'}:}{
#'      The proportion of module variance explained by the \code{'summary'}
#'      vector.
#'    }
#'  }
#'  When \code{simplify = TRUE} then the simplest possible structure will be 
#'  returned. E.g. if the network properties are requested for only one module 
#'  in only one dataset, then the returned list will have only the above elements. 
#'  
#'  When \code{simplify = FALSE} then a nested list of datasets will always be 
#'  returned, i.e. each element at the top level and second level correspond to
#'  a dataset, and each element at the third level will correspond to modules 
#'  discovered in the dataset specified at the top level if module labels are 
#'  provided in the corresponding \code{moduleAssignments} list element. E.g. 
#'  \code{results[["Dataset1"]][["Dataset2"]][["module1"]]} will contain the 
#'  properties of "module1" as calculated in "Dataset2", where "module1" was 
#'  indentified in "Dataset1". Modules and datasets for which calculation of 
#'  the network properties have not been requested will contain \code{NULL}.
#'  
#' @examples
#' \dontrun{
#' # load in example data, correlation, and network matrices for a discovery and test dataset:
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
#' # Calculate the topological properties of all network modules in the discovery dataset
#' props <- networkProperties(
#'   network=network_list, data=data_list, correlation=correlation_list, 
#'   moduleAssignments=labels_list
#' )
#'   
#' # Calculate the topological properties in the test dataset for the same modules
#' test_props <- networkProperties(
#'   network=network_list, data=data_list, correlation=correlation_list,
#'   moduleAssignments=labels_list, discovery="discovery", test="test"
#' )
#' }
#' 
#' @seealso \link[=nodeOrder]{Getting nodes ordered by degree.}, and
#'   \link[=sampleOrder]{Ordering samples by module summary}
#' 
#' @rdname networkProperties
#' @export
networkProperties <- function(
  network, data, correlation, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, simplify=TRUE, 
  verbose=TRUE
) {
  #----------------------------------------------------------------------------
  # Input processing and sanity checking
  #----------------------------------------------------------------------------
  vCat(verbose, 0, "Validating user input...")
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose)
  
  vCat(verbose, 0, "User input ok!")
  
  # Calculate the network properties
  res <- with(finput, {
    netPropsInternal(network, data, correlation, moduleAssignments, 
                     modules, discovery, test, nDatasets, datasetNames, verbose)
  })
  
  # Simplify the output data structure where possible
  if (simplify) {
    res <- simplifyList(res, depth=3)
  }
  on.exit({vCat(verbose, 0, "Done!")}, add=TRUE)
  res
}

#' Internal function for calculating 'networkProperties'
#' 
#' This function is used by several package functions. It assumes that all user
#' input has been processed already by \code{\link{processInput}}. This allows
#' for function-specific checking (i.e. failing early where the \code{'data'} is
#' required), while avoiding duplication of time-intensive checks 
#' (e.g. \code{\link{CheckFinite}}).
#' 
#' @param network \code{'network'} after processing by \code{'processInput'}.
#' @param data \code{'data'} after processing by \code{'processInput'}.
#' @param correlation \code{'correlation'} after processing by
#'   \code{'processInput'}.
#' @param moduleAssignments \code{'moduleAssignments'} after processing by
#'   \code{'processInput'}.
#' @param modules \code{'modules'} after processing by \code{'processInput'}.
#' @param discovery \code{'discovery'} after processing by
#'   \code{'processInput'}.
#' @param test \code{'test'} after processing by \code{'processInput'}.
#' @param nDatasets a vector containing the total number of input datasets, 
#'  returned by \code{'processInput'}.
#' @param datasetNames a vector of dataset names returned by
#'   \code{'processInput'}.
#' @param verbose logical; should progress be reported? Default is \code{TRUE}.
#'   
netPropsInternal <- function(
  network, data, correlation, moduleAssignments, modules, discovery, test,
  nDatasets, datasetNames, verbose
) {
  # The following declarations are for iterators declared inside each foreach 
  # loop. Declarations are required to satisfy NOTES generated by R CMD check, 
  # and also serve as useful documentation for the reader of the source code.
  di <- NULL # discovery dataset 
  ti <- NULL # test dataset 
  mi <- NULL # iterator over the modules
  
  # Set up results list
  res <- foreach(di = seq_len(nDatasets)) %do% {
    res2 <- foreach(ti = seq_len(nDatasets)) %do% {
      if (!is.null(moduleAssignments[[di]])) {
        allMods <- sortModuleNames(unique(moduleAssignments[[di]]))
        res3 <- foreach(mi = seq_along(allMods)) %do% {} 
        names(res3) <- allMods
        return(res3)
      }
    }
    names(res2) <- datasetNames
    return(res2)
  } 
  names(res) <- datasetNames
  
  props <- foreach(di = discovery) %do% {
    foreach(ti = test[[di]]) %do% {
      
      # Load matrices into RAM if they are 'big.matrix' objects.
      #   Calls to the garbage collector are necessary so we don't have 
      #   a copy of the data in shared memory as well.
      anyBM <- any.big.matrix(data[[ti]], network[[ti]])
      if (anyBM) {
        vCat(verbose, 0, "Loading 'big.matrix' data for dataset", 
             datasetNames[ti], "into RAM...")
      }
      data_mat <- loadIntoRAM(data[[ti]])
      gc()
      network_mat <- loadIntoRAM(network[[ti]])
      gc()
      
      vCat(verbose, 0, "Calculating network properties in dataset", 
           datasetNames[ti], "...")
      
      if (is.null(data[[ti]])) {
        NetworkPropertiesNoData(
          network_mat, moduleAssignments[[di]], modules[[di]]
        )
      } else {
        NetworkProperties(
          data_mat, network_mat, moduleAssignments[[di]], modules[[di]]
        )
      }
    }
  }
  
  # We populate the results list separately since they cannot be assigned 
  # directly in a parallel loop.
  for (ii in seq_along(props)) {
    for (jj in seq_along(props[[ii]])) {
      for (kk in seq_along(props[[ii]][[jj]])) {
        di <- discovery[ii]
        ti <- test[[di]][jj]
        mi <- as.character(modules[[di]][kk])
        res[[di]][[ti]][[mi]] <- props[[ii]][[jj]][[kk]]
      }
    }
  }
  
  return(res)
}
