#' Network Subset Data Summary
#'
#' Network subset eigenvector and proportion of variance explained in C++
#' For a given network subset, returns a vector that summarises the underlying
#' data that network subset was constructed from, along with the proportion of
#' the variance in the underlying data it explains.
#' 
#' @param data a \code{\link[bigmemory]{big.matrix}}: containing the data matrix 
#'  used to construct the network. Rows are the variables (nodes) and columns
#'  are the samples.
#' @param scaledData a \code{\link[bigmemory]{big.matrix}}: containing 
#'  \code{data} where each row (node) has been scaled. See 
#'  \code{\link[=scale-cpp]{scaleBigMatrix}}.
#' @param subsetIndices a vector of integers giving the row indices for the 
#'  network subset of interest in \code{data} and \code{scaledData}
#' 
#' @details
#' This returned vector is the first right singular vector from a 
#' \link[base][=svd]{singular value decomposition}. The vector returned may 
#' differ from that obtained through R's \code{\link[base]{svd}} due to 
#' different C/C++ libraries being used.
#' 
#' Only \code{\link[bigmemory]{big.matrix}} objects of \code{type = "double"}
#' can be used.
#' 
#' The \code{scaledData} is used to orient this vector so that its direction is
#' in alignment with the underlying data, as the SVD will sometimes return a
#' vector which is strongly negatively correlated, rather than positively 
#' correlated, with the data.
#' 
#' @section Warning:
#' A SVD cannot be calculated if there are any \code{NA}, \code{NaN}, or 
#' \code{Inf} in the supplied data. For speed reasons, this function does not
#' check the supplied data. If there are any non finite values in the data, 
#' the function will run for a long time, and fail, print the following 
#' message:
#'  \code{error: svd_econ(): failed to converge}
#' throw a warning, and return a list of \code{NA}s. 
#'
#' @param data
#' @export
dataSummary <- function(data, scaledData, subsetIndices) {
  return(DataSummary(data@address, scaledData@address, subsetIndices))
}
