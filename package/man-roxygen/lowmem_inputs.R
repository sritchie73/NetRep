#' @section Inputs:
#'   The expected object type for the input matrices depends on the value of 
#'   \code{lowmem}. When \code{lowmem} is \code{FALSE}, they are expected to be
#'   objects of type \code{\link[bigmemory]{big.matrix}}. When \code{TRUE}, they
#'   are expected to be file paths pointing to descriptor files for the 
#'   corresponding \code{\link[bigmemory]{big.matrix}} object. 
#'   
#'   The \code{discIndices} and \code{testIndices} arguments are expected to be
#'   integer vectors. If they contain missing values (\code{NA}) then C++
#'   will segfault, crashing your R session and potentionally losing all your
#'   data.
#' 
#' @section Minimising Memory Usage:
#'   When \code{lowmem} is set to \code{TRUE}, these functions will detach the
#'   \code{\link[bigmemory]{big.matrix}} objects after performing the 
#'   calculations, and invoke R's garbage collector, forcibly clearing the 
#'   virtual memory.
#' 
#' @param lowmem logical; See the Minimising Memory Usage
