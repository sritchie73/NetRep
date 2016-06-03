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
#' @param network list returned by \code{'processInput'}.
#' @param data data returned by \code{'processInput'}.
#' @param correlation list returned by \code{'processInput'}.
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
  network, data, correlation, moduleAssignments, modules, di,
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
      network, data, correlation, moduleAssignments, modules, di,
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
