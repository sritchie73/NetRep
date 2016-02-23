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
#'  For example, \code{results[[1]][[2]][["blue"]][["degree"]]} is a vector
#'  containing the \emph{weighted node degree} for the "blue" module from the
#'  dataset 1, as calculated in dataset 2. module preservation p-values when
#'  assessing the preservation of modules from dataset 1 in dataset 2. If
#'  \code{simplify = TRUE} then the list structure will be simplified where
#'  possible.
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
    cleanupCluster(par$cluster, par$predef)
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
    netPropsInternal(scaledData, correlation, network, moduleAssignments, 
                     modules, discovery, test, datasetNames, verbose)
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
#' (e.g. \code{\link{checkFinite}}) and data duplication (e.g. through
#' \code{\link{scaleBigMatrix}}).
#' 
#' @param scaledData scaled data generated by \code{'processInput'}.
#' @param correlation \code{'correlation'} after processing by
#'   \code{'processInput'}.
#' @param network \code{'network'} after processing by \code{'processInput'}.
#' @param moduleAssignments \code{'moduleAssignments'} after processing by
#'   \code{'processInput'}.
#' @param modules \code{'modules'} after processing by \code{'processInput'}.
#' @param discovery \code{'discovery'} after processing by
#'   \code{'processInput'}.
#' @param test \code{'test'} after processing by \code{'processInput'}.
#' @param datasetNames a vector of dataset names returned by
#'   \code{'processInput'}.
#' @param verbose logical; should progress be reported? Default is \code{TRUE}.
#'   
netPropsInternal <- function(
  scaledData, correlation, network, moduleAssignments, modules, discovery, test,
  datasetNames, verbose
) {
  #-----------------------------------------------------------------------------
  # Set up variables for running in parallel
  #-----------------------------------------------------------------------------
  # The following declarations are for iterators declared inside each foreach 
  # loop. Declarations are required to satisfy NOTES generated by R CMD check, 
  # and also serve as useful documentation for the reader of the source code.
  di <- NULL # discovery dataset 
  ti <- NULL # test dataset 
  mi <- NULL # iterator over the modules
  
  vCat(verbose, 0, 'Calculating properties for:\n')
  res <- foreach(di = discovery) %:% 
           foreach(ti = test[[di]]) %:% 
             foreach(mi = modules[[di]]) %dopar% {
    vCat(
      verbose, 1, sep="", 'Module "', mi, '" from dataset "', di, 
      '" in dataset "', ti, '"\n'
    )
    
    # Get the row/column indices of the module in the dataset of interest 
    sub <- moduleAssignments[[di]][moduleAssignments[[di]] == mi]
    modInds <- match(names(sub), rownames(correlation[[ti]]))
    na.inds <- which(is.na(modInds))
    modInds <- na.omit(modInds)
    
    # Get the properties calculated from the underlying data used to infer the
    # network
    datProps <- NULL
    if (!is.null(scaledData[[ti]])) {
      # We need to handle the case where no module variables are present in the
      # test dataset differently.
      if (length(modInds) > 0) {
        datProps <- dataProps(scaledData[[ti]], modInds)
        names(datProps) <- c("summary", "contribution", "coherence")
        datProps[[2]] <- insert.nas(datProps[[2]], na.inds)
      } else {
        dataProps <- list(
          summary=rep(NA, nrow(scaledData[[ti]])),
          contribution=rep(NA, length(sub)),
          coherence=NA
        )
      }
      names(datProps[["summary"]]) <- rownames(scaledData[[ti]])
      names(datProps[["contribution"]]) <- names(sub)
    } else {
      dataProps <- NULL
    }
    
    if (!is.null(network)) {
      # Get the properties calculated from the network.
      if (length(modInds) > 0) {
        netProps <- netProps(network[[ti]], modInds)
        names(netProps) <- c("degree", "avgWeight")
        netProps[[1]] <- insert.nas(netProps[[1]], na.inds)
      } else {
        netProps <- list(
          degree=rep(NA, length(sub)),
          avgWeight=NA
        )
      }
      names(netProps[["degree"]]) <- names(sub)
    } else {
      netProps <- NULL
    }
    return(c(datProps, netProps))
  }
  # Now we need to name the output 
  names(res) <- datasetNames[discovery]
  for (di in discovery) {
    names(res[[di]]) <- datasetNames[test[[di]]]
    for (ti in test[[di]]) {
      names(res[[di]][[ti]]) <- modules[[di]]
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
#' @param na.rm logical; If \code{TRUE}, nodes and moduels present in the 
#'   \code{discovery} dataset but missing from the test dataset are excluded. If
#'   \code{FALSE}, missing nodes and modules are put last in the ordering.
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
#'  requested module. If \code{simplify = TRUE}, then there will be a single
#'  vector of node names for each \code{'test'} dataset.
#'
#' @seealso \code{\link{networkProperties}}
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
  orderModules=TRUE, simplify=TRUE, verbose=TRUE
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
    cleanupCluster(par$cluster, par$predef)
  }, add=TRUE)
  
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
          if (is.null(scaledData[[ti]]) && length(modules[[di]]) > 1) {
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
        scaledData[di] <- list(NULL)
      }
    }
  })
  
  # Calculate the network properties
  props <- with(finput, { 
    netPropsInternal(
      scaledData, correlation, network, moduleAssignments, 
      modules, discovery, test, datasetNames, verbose
    ) 
  })
  
  res <- nodeOrderInternal(
    props, orderModules, simplify, verbose, na.rm
  )

  # Simplify the output data structure where possible
  if (simplify) {
    res <- simplifyList(res, depth=2)
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
#' 
#' @return list structure of ordered nodes.
nodeOrderInternal <- function(props, orderModules, simplify, verbose, na.rm) {
  vCat(verbose, 0, "Ordering nodes...")
  res <- lapply(props, function(discProps) { # For each discovery dataset 
    r1 <- lapply(discProps, function(testProps) {  # For each test dataset
      # First get module order
      if (length(testProps) > 1 && orderModules) {
        # Order modules by similarity
        
        # First get the summary profiles for each module
        summaries <- matrix(0, ncol=length(testProps), 
                            nrow=length(testProps[[1]][["summary"]]))
        colnames(summaries) <- names(testProps)
        rownames(summaries) <- names(testProps[[1]][["summary"]])
        for (mi in seq_along(testProps)) {
          summaries[,mi] <- testProps[[mi]][["summary"]]
        }
        
        # Identify modules where no nodes are present in the test dataset
        na.inds <- which(apply(summaries, 2, function(s) all(is.na(s))))
        na.mods <- colnames(summaries)[na.inds]
        if (length(na.mods) > 0) {
          summaries <- summaries[,-na.mods]
        }
        
        # Cluster modules that do have summary profiles
        clusteredMods <- hclust(as.dist(1-cor(summaries)))
        moduleOrder <- colnames(summaries)[clusteredMods$order]
        
        # Add the 'NA' modules
        moduleOrder <- c(moduleOrder, na.mods)
      } else {
        # Otherwise order as is
        moduleOrder <- names(testProps)
      }
      
      # next order nodes within each module
      nodeOrder <- foreach(mi = moduleOrder) %do% {
        nodeDegree <- testProps[[mi]][["degree"]]
        sortedNodes <- sort(nodeDegree, decreasing=TRUE, na.last=TRUE)
        sortedNodes
      }
      names(nodeOrder) <- moduleOrder
      
      # Remove missing nodes and modules 
      if (na.rm) {
        for (ii in rev(seq_along(nodeOrder))) {
          nodeOrder[[ii]] <- na.omit(nodeOrder[[ii]])
          if (length(nodeOrder[[ii]]) == 0) {
            nodeOrder[[ii]] <- NULL
          }
        }
      }
      
      # Now just get the names
      for (ii in seq_along(nodeOrder)) {
        nodeOrder[[ii]] <- names(nodeOrder[[ii]])
      }
      
      if (simplify) {
        nodeOrder <- unlist(nodeOrder)
        names(nodeOrder) <- NULL
      }
      return(nodeOrder)
    })
  })
  return(res)
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
#'  \code{'modules'} specified, containing a vector of sample names or indices 
#'  for the requested module.'
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
    cleanupCluster(par$cluster, par$predef)
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
  
  # Calculate the network properties. We don't actually need the network-based
  # properties though, so we can speed things up by ignoring them
  props <- with(finput, {
    netPropsInternal(scaledData, correlation, NULL, moduleAssignments, 
                     modules, discovery, test, datasetNames, verbose)
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
  res <- lapply(props, function(discProps) {
    lapply(discProps, function(testProps) {
      lapply(testProps, function(modProps) {
        summary <- modProps[["summary"]]
        if (na.rm) {
          summary <- na.omit(summary)
        }
        summaryOrder <- order(summary, decreasing=TRUE, na.last=TRUE)
        
        if(is.null(names(summary))) {
          return(summaryOrder)
        } else {
          return(names(summary)[summaryOrder])
        }
      })
    })
  })
  return(res)
}
