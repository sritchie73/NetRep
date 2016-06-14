#' Combine results of multiple permutation procedures
#' 
#' This function takes the output from multiple runs of 
#' \code{\link{modulePreservation}}, combines their results, and returns a new
#' set of permutation test P-values. This is useful for parallelising 
#' calculations across multiple machines.
#' 
#' @param pres1,pres2 lists returned by \code{\link{modulePreservation}}.
#' 
#' @details
#'  The calls to 'modulePreservation' must have been identical for both input
#'  lists, with the exception of the number of threads used and the number of
#'  permutations calculated.
#' 
#' @return
#'  A nested list containing the same elements as 
#'  \code{\link{modulePreservation}}.
#' 
#' @examples 
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
#' pres1 <- modulePreservation(
#'  network=network_list, data=data_list, correlation=correlation_list,
#'  moduleAssignments=labels_list, nPerm=5000, discovery="discovery", 
#'  test="test", nThreads=2
#' )
#' 
#' pres2 <- modulePreservation(
#'  network=network_list, data=data_list, correlation=correlation_list, 
#'  moduleAssignments=labels_list, nPerm=5000, discovery="discovery", 
#'  test="test", nThreads=2
#' )
#' 
#' combined <- combineAnalyses(pres1, pres2)
#' 
#' @importFrom abind abind
#' @export
combineAnalyses <- function(pres1, pres2) {
  msg <- "module preservation analyses differ between 'pres1' and 'pres2'"
  
  if ("observed" %in% names(pres1)) {
    pres1 <- combineAnalysesInternal(pres1, pres2)
  } else {
    for (ii in seq_along(pres1)) {
      if ("observed" %in% names(pres1[[ii]])) {
        if (!is.null(pres1[[ii]]) && !is.null(pres2[[ii]])) {
          pres1[[ii]] <- combineAnalysesInternal(pres1[[ii]], pres2[[ii]])
        } else if (xor(is.null(pres1[[ii]]), is.null(pres2[[ii]]))) {
          stop(msg)
        }
        
      } else {
        for (jj in seq_along(pres1[[ii]])) {
          if (!is.null(pres1[[ii]][[jj]]) && !is.null(pres2[[ii]][[jj]])) {
            pres1[[ii]][[jj]] <- combineAnalysesInternal(pres1[[ii]][[jj]], pres2[[ii]][[jj]])
          } else if (xor(is.null(pres1[[ii]][[jj]]), is.null(pres2[[ii]][[jj]]))) {
            stop(msg)
          }
        }
      }
    }
  }
  
  return(pres1)
}

#' Combine results of multiple permutation procedures
#' 
#' This function takes the module preservation analyses for two dataset 
#' comparisons and combines them.
#' 
#' @param pres1 input from 'combineAnalyses'
#' @param pres2 input from 'combineAnalyses'
#' 
#' @return
#'  A combined module preservation analysis.
#' 
combineAnalysesInternal <- function(pres1, pres2) {
  msg <- "module preservation analyses differ between 'pres1' and 'pres2'"
  badInput <- FALSE
  tryCatch({
    if (nrow(pres1$observed) != nrow(pres2$observed) ||
        ncol(pres1$observed) != ncol(pres2$observed) ||
        pres1$alternative != pres2$alternative ||
        pres1$totalSize != pres2$totalSize ||
        pres1$propVarsPresent != pres2$propVarsPresent ||
        pres1$nVarsPresent != pres1$nVarsPresent ||
        (!is.null(pres1$contingency) && is.null(pres1$contingency)) ||
        (is.null(pres1$contingency) && !is.null(pres1$contingency)) ||
        (!is.null(pres1$contingency) && !is.null(pres1$contingency) && 
         !all.equal(pres1$contingency, pres2$contingency)) &&
        !all.equal(pres1$observed, pres2$observed)) {
      badInput <- TRUE
    }
  }, error=function(e) {
    stop("'pres1' and 'pres2' do not appear to be output from 'modulePreservation'")
  })
  if (badInput) {
    stop("module preservation analysis performed in 'pres1' and 'pres2'",
         " are not comparable")
  }
  
  res <- pres1
  res$nulls <- abind::abind(pres1$nulls, pres2$nulls, along=3)
  res$p.values <- permutationTest(res$nulls, res$observed, res$nVarsPresent,
                                  res$totalSize, res$alternative)
  return(res)
}