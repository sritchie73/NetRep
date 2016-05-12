#' Multi machine analysis
#' 
#' This function takes the output from multiple runs of 
#' \code{\link{modulePreservation}}, combines their results, and returns a new
#' set of permutation test P-values. This is useful for parallelising 
#' calculations across multiple machines.
#' 
#' @param pres1,pres2 lists returned by \code{\link{modulePreservation}}.
#' 
#' @details
#'  The input arguments must be for a single dataset comparison, i.e. the list
#'  must contain the elements \code{nulls} for the null distributions. If 
#'  multiple \code{discovery} or \code{test} datasets have been specified, you
#'  will need to extract the relevant list elements 
#' 
#' @return
#'  A nested list containing the same elements as 
#'  \code{\link{modulePreservation}}.
#' 
#' @examples 
#' \dontrun{
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
#' # Assess module preservation on two different machines to spread computational
#' # burden.
#' pres1 <- modulePreservation(
#'  data=data_list, correlation=correlation_list, network=network_list,
#'  moduleAssignments=labels_list, nPerm=500, discovery="discovery", 
#'  test="test"
#' )
#' 
#' pres2 <- modulePreservation(
#'  data=data_list, correlation=correlation_list, network=network_list,
#'  moduleAssignments=labels_list, nPerm=500, discovery="discovery", 
#'  test="test"
#' )
#' 
#' combined <- combineAnalyses(pres1, pres2)
#' 
#' }
#' 
#' @importFrom abind abind
#' @export
combineAnalyses <- function(pres1, pres2) {
  # Validate use input
  msg <- "module preservation analyses differ between 'pres1' and 'pres2'"
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
      stop("module preservation analysis performed in 'pres1' and 'pres2'",
           " are not comparable")
    }
  }, error=function(e) {
    stop("'pres1' and 'pres2' do not appear to be output from 'modulePreservation'")
  })
  
  res <- pres1
  res$nulls <- abind::abind(pres1$nulls, pres2$nulls, along=3)
  res$p.values <- permutationTest(res$nulls, res$observed, res$nVarsPresent,
                                  res$totalSize, res$alternative)
  return(res)
}