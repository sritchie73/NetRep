#' Permutation test P-values for module preservation statistics
#' 
#' Evaluates the statistical significance of each module preservation test 
#' statistic for one or more modules.
#'  
#' @details 
#'  Calculates exact p-values for permutation tests when permutations are 
#'  randomly drawn with replacement using the \code{\link[statmod]{permp}} 
#'  function in the \code{\link{statmod}} package. 
#'  
#'  This function may be useful for re-calculating permutation test P-values,
#'  for example when there are missing values due to sparse data. In this case
#'  the user may decide that these missing values should be assigned 0 so that
#'  P-values aren't signficant purely due to many incalculable statistics leading
#'  to low power.
#'  
#' @references 
#'   \enumerate{
#'     \item{
#'       Phipson, B. & Smyth, G. K. \emph{Permutation P-values should never be 
#'       zero: calculating exact P-values when permutations are randomly drawn.}
#'       Stat. Appl. Genet. Mol. Biol. \strong{9}, Article39 (2010). 
#'     }
#'   }
#'  
#' @param nulls a 3-dimension matrix where the columns correspond to module
#'   preservation statistics, rows correspond to modules, and the third 
#'   dimension to null distribution observations drawn from the permutation 
#'   procedure in \code{\link{modulePreservation}}.
#' @param observed a matrix of observed values for each module preservation
#'  statistc (columns) for each module (rows) returned from 
#'  \code{\link{modulePreservation}}.
#' @param nVarsPresent a vector containing the number of variables/nodes in each
#'  module that was present in the \emph{test} dataset. Returned as a list 
#'  element of the same name by \code{\link{modulePreservation}}.
#' @param totalSize the size of the test network used to perform the test. 
#'  Returned as a list element of the same name by 
#'  \code{\link{modulePreservation}}.
#' @param alternative a character string specifying the alternative hypothesis, 
#'  must be one of "greater" (default), "less", or "two.sided". 
#'  You can specify just the initial letter.
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
#' # Note that we recommend running at least 10,000 permutations to make sure 
#' # that the null distributions are representative.
#'
#' preservation <- modulePreservation(
#'  network=network_list, data=data_list, correlation=correlation_list, 
#'  moduleAssignments=labels_list, nPerm=nPerm, discovery="discovery", 
#'  test="test"
#' )
#' 
#' # Re-calculate the permutation test P-values
#' p.values <- permutationTest(
#'   preservation$nulls, preservation$observed, preservation$nVarsPresent,
#'   preservation$totalSize, preservation$alternative
#' )
#' 
#' @aliases permutation permuted
#' @name permutationTest
#' @importFrom statmod permp
#' @keywords internal
#' @export
permutationTest <- function(
  nulls, observed, nVarsPresent, totalSize, alternative="greater"
) {
  # Validate user input
  validAlts <- c("two.sided", "less", "greater")
  altMatch <- pmatch(alternative, validAlts)
  if (is.na(altMatch))
    stop("Alternative must be one of ", validAlts)
  
  if (!is.numeric(totalSize) || length(totalSize) > 1 || totalSize < 1)
    stop("'totalSize' must be a single number > 0")
  
  statNames <- c("avg.weight", "coherence", "cor.cor", "cor.degree", 
                 "cor.contrib", "avg.cor", "avg.contrib")
  
  if (!is.matrix(observed) || ncol(observed) %nin% c(4,7) || 
      any(colnames(observed) %nin% statNames) || !is.numeric(observed)) {
    stop("expecting 'observed' to be a numeric matrix output by the ", 
         "'modulePreservation' function")
  }
  if (!is.array(nulls) || ncol(nulls) %nin% c(4,7) || length(dim(nulls)) != 3 ||
      any(colnames(nulls) %nin% statNames) || !is.numeric(nulls)) {
    stop("expecting 'nulls' to be a numeric matrix output by the ", 
         "'modulePreservation' function")
  }
  if (any(colnames(nulls) != colnames(observed)) || 
      any(rownames(nulls) != rownames(observed)) ||
      any(colnames(observed) != colnames(nulls)) || 
      any(rownames(observed) != rownames(nulls))) {
    stop("mismatch in dimension names between 'nulls' and 'observed'")
  }
  
  if (!is.numeric(nVarsPresent) || length(nVarsPresent) != nrow(nulls) ||
      any(rownames(nulls) != names(nVarsPresent)) || 
      any(names(nVarsPresent) != rownames(nulls))) {
    stop("expecting 'nVarsPresent' to be a numeric vector output by the ",
         "'modulePreservation' function")
  }
  
  # Calculate module preservation statistic p-values
  p.values <- matrix(NA, nrow(nulls), ncol(nulls), dimnames=dimnames(observed))
  for (mi in seq_len(nrow(p.values))) {
    for (si in seq_len(ncol(p.values))) {
      # If the observed value is missing, leave the p-value missing.
      if (is.na(observed[mi, si])) {
        next
      }
      
      # Node order in the sampling does not matter when calculating the average
      # edge weight or module coherence
      if (colnames(observed)[si] %in% c("avg.weight", "coherence")) {
        order <- FALSE
      } else {
        order <- TRUE
      }
      
      # This in turn affects the total number of possible permutations
      if (order) {
        total.nperm = prod(totalSize:(totalSize - nVarsPresent[mi] + 1))
      } else {
        total.nperm = choose(totalSize, nVarsPresent[mi])
      }
      
      # Calculate necessary components to perform any of the alternative tests
      permuted <- sort(nulls[mi,si,])
      nPerm <- length(permuted)
      less.extreme <- length(permuted[permuted <= observed[mi, si]])
      more.extreme <- length(permuted[permuted >= observed[mi, si]])
      lower.pval <- permp(less.extreme, nPerm, total.nperm=total.nperm)
      upper.pval <- permp(more.extreme, nPerm, total.nperm=total.nperm)
      
      if (altMatch == 1L) {
        p.values[mi, si] <- min(lower.pval, upper.pval)*2
      } else if (altMatch == 2L) {
        p.values[mi, si] <- lower.pval
      } else if (altMatch == 3L) {
        p.values[mi, si] <- upper.pval
      } 
    }
  }
  # Check for missing values that aren't due to a module not being present
  missingMods <- apply(observed, 1, function(x) all(is.na(x)))
  if (any(is.na(observed[!missingMods,])) || any(is.na(nulls[!missingMods,,]))) {
    warning(
      "Missing values encountered in the observed test statistics and/or ",
      "in their null distributions. P-values may be biased for these tests.",
      " See 'help(", '"permutationTest"', ")'", 
      immediate. = TRUE
    )

  }
  
  return(p.values)
}

### Exact permutation p-values wrapper
### 
### Wrapper for \code{\link[statmod]{permp}} from the 
### \code{\link[statmod]{statmod}} library, which can crash if FORTRAN 
### libraries are not properly linked.
### 
### @details
### In the case \code{\link[statmod]{permp}} fails, the wrapper will fall back 
### to a slightly more conservative biased estimator: (1+x)/(1+nPerm).
### 
### @param x number of permutations that yielded test statistics at least as 
###  extreme as the observed data. May be a vector or an array of values. 
### @param nperm total number of permutations performed.
### @param ... other arguments to pass to\code{\link[statmod]{permp}}.
### @return 
###  vector or array of p-values, of same dimensions as \code{x}.
###  
### @importFrom statmod permp
### @keywords internal
permp <- function(x, nperm, ...) {
  tryCatch({
    return(statmod::permp(x, nperm, ...))
  }, error=function(e) {
    warning(
      "Error from statmod::permp:", e$message, 
      "\nUsing conservative biased estimator (1+x)/(1+nPerm) instead.",
      immediate. = TRUE
    )
    return(
      (x + 1)/(nperm + 1)  
    )
  })
}

#' How many permutations do I need to test at my desired significance level?
#' 
#' @param alpha desired significance threshold.
#' 
#' @return The minimum number of permutations required to detect any significant
#'  associations at the provided \code{alpha}. The minimum p-value will always
#'  be smaller than \code{alpha}.
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
#' # How many permutations are required to Bonferroni adjust for the 4 modules 
#' # in the example data? 
#' nPerm <- requiredPerms(0.05/4) 
#' 
#' # Note that we recommend running at least 10,000 permutations to make sure 
#' # that the null distributions are representative.
#'
#' preservation <- modulePreservation(
#'  network=network_list, data=data_list, correlation=correlation_list, 
#'  moduleAssignments=labels_list, nPerm=nPerm, discovery="discovery", 
#'  test="test"
#' )
#' 
#' @export
requiredPerms <- function(alpha, alternative="greater") {
  validAlts <- c("two.sided", "less", "greater")
  altMatch <- pmatch(alternative, validAlts)
  if (is.na(altMatch))
    stop("Alternative must be one of ", validAlts)
  
  if (altMatch == 1) {
    1/alpha*2
  } else {
    1/alpha
  }
}
