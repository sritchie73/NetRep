#' Calculate the topological properties for a network module
#' 
#' Calculates the network properties used to assess module preservation for one
#' or more modules in a user specified dataset.
#' 
#' @inheritParams common_params
#' @inheritParams simplify_param
#' @inheritParams par_param
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
#' \subsection{'bigMatrix' input data:}{
#'   Although the \code{data}, \code{correlation}, and \code{network} arguments 
#'   expect data to be provided in the \code{\link{bigMatrix}}, they can be 
#'   provided as regular \code{\link[base]{matrix}} objects, in which case they 
#'   will be temporarily converted to \code{bigMatrix} objects. \strong{This is 
#'   not recommended}, as each matrix will be copied into memory of each 
#'   parallel R session, resulting in much higher memory usage (see section on
#'   Memory usage ), and increased computation time for the conversion and
#'   copying processes. It is therefore strongly recommended that the user save
#'   their data separately as \code{\link{bigMatrix}} objects prior to running
#'   the permutation procedure or using any other package function. This is also
#'   useful for other analyses, as \code{bigMatrix} objects can be
#'   instantaneously loaded into any future R session.
#'   
#'   Alternatively, the \code{data}, \code{correlation}, and \code{network} 
#'   arguments will also accept file paths to tabular data or \code{bigMatrix} 
#'   backingfiles.
#' }
#' \subsection{Memory usage:}{
#'   Provided there are no additional objects in the R session the
#'   permutation procedure will use only the memory required to store each 
#'   matrix once, along with an additional 200 MB per core used by each vanilla
#'   R session.
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
#' # Convert them to the 'bigMatrix' format:
#' discovery_data <- as.bigMatrix(discovery_data)
#' discovery_correlation <- as.bigMatrix(discovery_correlation)
#' discovery_network <- as.bigMatrix(discovery_network)
#' test_data <- as.bigMatrix(test_data)
#' test_correlation <- as.bigMatrix(test_correlation)
#' test_network <- as.bigMatrix(test_network)
#' 
#' # Set up input lists for each input matrix type across datasets:
#' data_list <- list(discovery=discovery_data, test=test_data)
#' correlation_list <- list(discovery=discovery_correlation, test=test_correlation)
#' network_list <- list(discovery=discovery_network, test=test_network)
#' labels_list <- list(discovery=module_labels)
#' 
#' # Calculate the topological properties of all network modules in the discovery dataset
#' props <- networkProperties(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list
#' )
#'   
#' # Calculate the topological properties in the test dataset for the same modules
#' test_props <- networkProperties(
#'   data=data_list, correlation=correlation_list, network=network_list, 
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
  data=NULL, correlation, network, moduleAssignments=NULL, modules=NULL,
  backgroundLabel="0", discovery=NULL, test=NULL, nCores=NULL, simplify=TRUE, 
  verbose=TRUE
) {
  #-----------------------------------------------------------------------------
  # Input processing and sanity checking
  #-----------------------------------------------------------------------------
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=FALSE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef, par$oldOMPThreads, par$oldBLASThreads)
  }, add=TRUE)

  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir)
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
  vCat(verbose, 0, "User input ok!")
  
  # Calculate the network properties
  res <- with(finput, {
    netPropsInternal(data, correlation, network, moduleAssignments, 
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
#' (e.g. \code{\link{checkFinite}}).
#' 
#' @param data \code{'data'} after processing by \code{'processInput'}.
#' @param correlation \code{'correlation'} after processing by
#'   \code{'processInput'}.
#' @param network \code{'network'} after processing by \code{'processInput'}.
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
  data, correlation, network, moduleAssignments, modules, discovery, test,
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
  
  vCat(verbose, 0, 'Calculating network properties...\n')
  props <- foreach(di = discovery) %do% {
    foreach(ti = test[[di]]) %do% {
      if (is.null(data[[ti]])) {
        NetworkPropertiesNoData(network[[ti]][,], moduleAssignments[[di]], 
                                modules[[di]])
      } else {
        NetworkProperties(data[[ti]][,], network[[ti]][,], moduleAssignments[[di]], 
                          modules[[di]])
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


#' Order nodes and modules within a network.
#' 
#' Order nodes in descending order of \emph{weighted degree} and order 
#' modules by the similarity of their summary vectors.
#' 
#' @inheritParams common_params
#' @inheritParams simplify_param
#' @inheritParams par_param
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
#' \subsection{'bigMatrix' input data:}{
#'   Although the \code{data}, \code{correlation}, and \code{network} arguments 
#'   expect data to be provided in the \code{\link{bigMatrix}}, they can be 
#'   provided as regular \code{\link[base]{matrix}} objects, in which case they 
#'   will be temporarily converted to \code{bigMatrix} objects. \strong{This is 
#'   not recommended}, as each matrix will be copied into memory of each 
#'   parallel R session, resulting in much higher memory usage (see section on
#'   Memory usage ), and increased computation time for the conversion and
#'   copying processes. It is therefore strongly recommended that the user save
#'   their data separately as \code{\link{bigMatrix}} objects prior to running
#'   the permutation procedure or using any other package function. This is also
#'   useful for other analyses, as \code{bigMatrix} objects can be
#'   instantaneously loaded into any future R session.
#'   
#'   Alternatively, the \code{data}, \code{correlation}, and \code{network} 
#'   arguments will also accept file paths to tabular data or \code{bigMatrix} 
#'   backingfiles.
#' }
#' \subsection{Memory usage:}{
#'   Provided there are no additional objects in the R session the
#'   permutation procedure will use only the memory required to store each 
#'   matrix once, along with an additional 200 MB per core used by each vanilla
#'   R session.
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
#' # Convert them to the 'bigMatrix' format:
#' discovery_data <- as.bigMatrix(discovery_data)
#' discovery_correlation <- as.bigMatrix(discovery_correlation)
#' discovery_network <- as.bigMatrix(discovery_network)
#' test_data <- as.bigMatrix(test_data)
#' test_correlation <- as.bigMatrix(test_correlation)
#' test_network <- as.bigMatrix(test_network)
#' 
#' # Set up input lists for each input matrix type across datasets:
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
  backgroundLabel="0", discovery=NULL, test=NULL, nCores=NULL, na.rm=FALSE, 
  orderModules=TRUE, mean=FALSE, simplify=TRUE, verbose=TRUE
) {
  #-----------------------------------------------------------------------------
  # Input processing and sanity checking
  #-----------------------------------------------------------------------------
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  vCat(verbose, 0, "Validating user input...")
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=FALSE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef, par$oldOMPThreads, par$oldBLASThreads)
  }, add=TRUE)
  
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
                         verbose, tmp.dir)
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
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

#' Order samples within a network.
#' 
#' Get the order of samples within a module based on the module summary vector.
#' 
#' @inheritParams common_params
#' @inheritParams simplify_param
#' @inheritParams par_param
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
#' \subsection{'bigMatrix' input data:}{
#'   Although the \code{data}, \code{correlation}, and \code{network} arguments 
#'   expect data to be provided in the \code{\link{bigMatrix}}, they can be 
#'   provided as regular \code{\link[base]{matrix}} objects, in which case they 
#'   will be temporarily converted to \code{bigMatrix} objects. \strong{This is 
#'   not recommended}, as each matrix will be copied into memory of each 
#'   parallel R session, resulting in much higher memory usage (see section on
#'   Memory usage ), and increased computation time for the conversion and
#'   copying processes. It is therefore strongly recommended that the user save
#'   their data separately as \code{\link{bigMatrix}} objects prior to running
#'   the permutation procedure or using any other package function. This is also
#'   useful for other analyses, as \code{bigMatrix} objects can be
#'   instantaneously loaded into any future R session.
#'   
#'   Alternatively, the \code{data}, \code{correlation}, and \code{network} 
#'   arguments will also accept file paths to tabular data or \code{bigMatrix} 
#'   backingfiles.
#' }
#' \subsection{Memory usage:}{
#'   Provided there are no additional objects in the R session the
#'   permutation procedure will use only the memory required to store each 
#'   matrix once, along with an additional 200 MB per core used by each vanilla
#'   R session.
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
#' # Convert them to the 'bigMatrix' format:
#' discovery_data <- as.bigMatrix(discovery_data)
#' discovery_correlation <- as.bigMatrix(discovery_correlation)
#' discovery_network <- as.bigMatrix(discovery_network)
#' test_data <- as.bigMatrix(test_data)
#' test_correlation <- as.bigMatrix(test_correlation)
#' test_network <- as.bigMatrix(test_network)
#' 
#' # Set up input lists for each input matrix type across datasets:
#' data_list <- list(discovery=discovery_data, test=test_data)
#' correlation_list <- list(discovery=discovery_correlation, test=test_correlation)
#' network_list <- list(discovery=discovery_network, test=test_network)
#' labels_list <- list(discovery=module_labels)
#' 
#' # Sort nodes within module 1 in descending order by module summary
#' samples <- sampleOrder(
#'   data=data_list, correlation=correlation_list, network=network_list, 
#'   moduleAssignments=labels_list, modules="1" 
#' )
#' }
#' 
#' @seealso \code{\link{networkProperties}}
#' 
#' @name sampleOrder
#' @export
sampleOrder <- function(
  data=NULL, correlation, network, moduleAssignments=NULL, modules=NULL, 
  backgroundLabel="0", discovery=NULL, test=NULL, nCores=NULL, na.rm=FALSE, 
  simplify=TRUE, verbose=TRUE
) {
  #-----------------------------------------------------------------------------
  # Input processing and sanity checking
  #-----------------------------------------------------------------------------
  tmp.dir <- file.path(tempdir(), paste0(".NetRep", getUUID()))
  dir.create(tmp.dir, showWarnings=FALSE)
  
  # Register parallel backend. 
  par <- setupParallel(nCores, verbose, reporterCore=FALSE)
  nCores <- par$nCores
  on.exit({
    cleanupCluster(par$cluster, par$predef, par$oldOMPThreads, par$oldBLASThreads)
  }, add=TRUE)
  
  vCat(verbose, 0, "Validating user input...")
  
  if (!is.logical(na.rm) || is.na(na.rm) || length(na.rm) > 1) {
    stop("'na.rm' must be either 'TRUE' or 'FALSE'")
  }
  
  # Now try to make sense of the rest of the input
  finput <- processInput(discovery, test, network, correlation, data, 
                         moduleAssignments, modules, backgroundLabel,
                         verbose, tmp.dir)
  on.exit({
    vCat(verbose, 0, "Cleaning up temporary objects...")
    unlink(tmp.dir, recursive = TRUE)
  }, add = TRUE)
  
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
    netPropsInternal(data, correlation, network, moduleAssignments, 
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

#' Filter structure returned by 'netPropsInternal' to specified test datasets
#' 
#' Used in the plot functions to filter the structure to the 
#' \code{'orderNodesBy'} and \code{'orderSamplesBy'} \emph{test} datasets.
#' 
#' @param props nested list returned by \code{\link{netPropsInternal}}.
#' @param test a vector of datasets to filter to.
#' @param discovery a vector containing a single discovery dataset specified in
#'   the parent plot function. 
#' @param modules a vector of modules from the specified \code{discovery} 
#'   dataset to filter to. If \code{NULL}, then all modules are kept.
#'   
#' @return 
#'  a nested list structure identical to that returned by 
#'  \code{'netPropsInternal'}, but where only the entries for test datasets
#'  specified by the \code{test} argument and modules specified by the 
#'  \code{modules} argument contain non-\code{NULL} entries.
filterInternalProps <- function(props, test, discovery, modules=NULL) {
  ii <- NULL; jj <- NULL; mi <- NULL # suppress CRAN note
  fProps <- foreach(ii = seq_along(props)) %do% {
    fProps2 <- foreach(jj = seq_along(props[[ii]])) %do% {
      if (names(props)[ii] == discovery) {
        fProps3 <- foreach(mi = seq_along(props[[ii]][[jj]])) %do% {
          if (names(props[[ii]])[jj] %in% test && 
              (is.null(modules) || names(props[[ii]][[jj]])[mi] %in% modules)) {
            return(props[[ii]][[jj]][[mi]])
          }
        }
        names(fProps3) <- names(props[[ii]][[jj]])
        return(fProps3)
      }
    }
    names(fProps2) <- names(props[[ii]])
    return(fProps2)
  }
  names(fProps) <- names(props)
  return(fProps)
}

#' Get the network properties and order for a plot
#'
#' @param data data returned by \code{'processInput'}.
#' @param correlation list returned by \code{'processInput'}.
#' @param network list returned by \code{'processInput'}.
#' @param moduleAssignments list returned by \code{'processInput'}.
#' @param modules vector of modules to show on the plot.
#' @param di name of the discovery dataset.
#' @param ti name of the test dataset.
#' @param orderNodesBy vector returned by \code{'processInput'}.
#' @param orderSamplesBy vector returned by \code{'processInput'}.
#' @param orderModules vector returned by \code{'checkPlotArgs'}.
#' @param datasetNames vector returned by \code{'processInput'}.
#' @param nDatasets vector returned by \code{'processInput'}.
#' @param dryRun logical; are we just doing a dry run of the plot?
#' @param verbose logical; turn on verbose printing.
#'
plotProps <- function(
  data, correlation, network, moduleAssignments, modules, di,
  ti, orderNodesBy, orderSamplesBy, orderModules, datasetNames, nDatasets, 
  dryRun, verbose
) {
  mods <- modules[[di]]
  mi <- NULL # suppresses CRAN note
  
  if (dryRun) {
    # If doing a dry run just get the nodes and samples that will be shown on
    # the plot in any order.
    moduleOrder <- mods
    nodeOrder <- unlist(sapply(mods, function(mi) {
      names(moduleAssignments[[di]][moduleAssignments[[di]] == mi])
    }))
    if (identical(orderNodesBy, ti)) {
      nodeOrder <- intersect(nodeOrder, colnames(network[[ti]]))
    }
    
    if (is.null(orderSamplesBy)) {
      sampleOrder <- NULL
    } else if (!is.na(orderSamplesBy)) {
        sampleOrder <- rownames(data[[orderSamplesBy]])
    } else {
      sampleOrder <- rownames(data[[ti]])
    }
    testProps <- NULL
  } else {
    # Scenarios:
    # - No ordering of nodes + samples. We only need to calculate the network 
    #   properties for the 'test' dataset.
    # - Ordering of nodes only. We need to calculate the network properties in
    #   all datasets specified in 'orderNodesBy' (may be one or more) and in the
    #   'test' dataset (may or may not be specified in 'orderNodesBy').
    # - Ordering of samples only. We need to calculate the network properties in
    #   the 'orderSamplesBy' dataset, and in the 'test' dataset (which may or 
    #   may not be the same as 'orderSamplesBy').
    # - Ordering of both. We need to calculate the network properties in the
    #   'orderSamplesBy', 'orderNodesBy', and 'test' datasets.
    # this vector contains all datasets required for plotting
    
    plotDatasets <- list(unique(na.omit(c(ti, orderSamplesBy, orderNodesBy))))
    names(plotDatasets) <- datasetNames[di]
    
    # Calculate the network properties for all datasets required
    props <- netPropsInternal(
      data, correlation, network, moduleAssignments, modules, di,
      plotDatasets, nDatasets, datasetNames, FALSE
    )
    
    # Order nodes based on degree
    if (length(orderNodesBy) > 1 || !is.na(orderNodesBy)) {
      if (length(orderNodesBy) > 1) {
        mean <- TRUE
      } else {
        mean <- FALSE
      }
      
      # nodeOrderInternal will average acros all test datasets, so we need to 
      # filter just to those specified in 'orderNodesBy' while preserving the
      # structure of 'props'
      orderProps <- filterInternalProps(props, orderNodesBy, di)
      nodeOrder <- nodeOrderInternal(
        orderProps, orderModules, simplify=FALSE, verbose, na.rm=FALSE, mean
      )
      nodeOrder <- simplifyList(nodeOrder, depth=3)
      
      # The module order will be the names of the simplified list iff there are
      # multiple modules to render
      if (!is.list(nodeOrder)) {
        moduleOrder <- mods
        if (is.numeric(moduleOrder))
          moduleOrder <- as.character(moduleOrder)
      } else {
        moduleOrder <- names(nodeOrder)
      }
      
      # Now flatten the node order list
      nodeOrder <- unlist(nodeOrder)
    } else {
      hasProps <- !sapply(props[[di]][[ti]], is.null) 
      moduleOrder <- names(props[[di]][[ti]])[hasProps]
      nodeOrder <- foreach(mi = moduleOrder, .combine=c) %do% {
        names(props[[di]][[ti]][[mi]]$degree)
      }
    }
    
    if (is.null(orderSamplesBy)) {
      sampleOrder <- NULL
    } else if (!is.na(orderSamplesBy)) {
      orderProps <- filterInternalProps(props, orderSamplesBy, di, moduleOrder[1])
      sampleOrder <- sampleOrderInternal(orderProps, verbose, na.rm=FALSE)
      sampleOrder <- simplifyList(sampleOrder, depth=3)
    } else {
      sampleOrder <- rownames(data[[ti]])
    }
    
    # Just keep the properties we need for plotting
    testProps <- simplifyList(props[[di]][[ti]], depth=1)
    if (length(moduleOrder) == 1) {
      testProps <- list(testProps)
      names(testProps) <- moduleOrder
    }
  }
  
  #-----------------------------------------------------------------------------
  # Identify nodes and samples from the 'discovery' dataset not present in the 
  # 'test' dataset.
  #-----------------------------------------------------------------------------
  
  na.pos.x <- which(nodeOrder %nin% colnames(network[[ti]]))
  if (length(na.pos.x) > 0) {
    presentNodes <- nodeOrder[-na.pos.x]
  } else {
    presentNodes <- nodeOrder
  }
  
  if (is.null(sampleOrder)) {
    na.pos.y <- NULL
    presentSamples <- NULL
  } else if (!is.numeric(sampleOrder)) {
    na.pos.y <- which(sampleOrder %nin% rownames(data[[ti]]))
    if (length(na.pos.y) > 0) {
      presentSamples <- sampleOrder[-na.pos.y]
    } else {
      presentSamples <- sampleOrder
    }
  } else {
    na.pos.y <- vector()
    presentSamples <- sampleOrder
  }
  

  # Returned processed results
  return(list(
    testProps=testProps, nodeOrder=nodeOrder, moduleOrder=moduleOrder,
    sampleOrder=sampleOrder, na.pos.x=na.pos.x, na.pos.y=na.pos.y,
    presentNodes=presentNodes, presentSamples=presentSamples
  ))
}


