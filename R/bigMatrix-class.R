# Allows the matrix slot of the big.matrix class to be empty, which
# is desirable if we're forcing modulePreservation to consume minimal memory.
setClassUnion("maybe.big.matrix", c("big.matrix", "NULL"))

# Allows the rownames and colnames slots of the big.matrix class to be empty.
setClassUnion("optional.dimnames", c("character", "NULL"))

#' The 'bigMatrix' class
#' 
#' A \code{'bigMatrix'} is a numeric matrix that is stored in shared memory. 
#' This memory is accessible from multiple parallel R sessions, allowing for 
#' massive parallisation of \code{\link{modulePreservation}} with constant 
#' memory usage. These matrices are stored on disk, and dynamically loaded as 
#' required, allowing instantaneous loading and access of very large matrices
#' in new R sessions.
#' 
#' @details
#' \subsection{Interacting with \code{bigMatrix} objects}{
#' \code{'bigMatrix'} objects can be accessed and manipulated with the same 
#' standard set of functions as \code{'matrix'} objects. These include: 
#' \code{head}, \code{tail}, \code{print}, \code{show}, \code{dimnames},
#' \code{rownames}, \code{colnames}, \code{nrow}, \code{ncol}, \code{dim}, 
#' and element access with \code{[} (see examples). \code{is.bigMatrix} will
#' check whether an object is a \code{bigMatrix}.
#' }
#' \subsection{\code{bigMatrix} types}{
#' The numeric data within a \code{bigMatrix} can be stored either as a 
#' \code{"double"}, \code{"float"}, \code{"integer"}, \code{"short"}, or
#' \code{"char"}. The default is as \code{"double"}: which keeps the values
#' within the \code{'bigMatrix'} at the same numerical precision as when stored
#' in a regular \code{'matrix'}. Storing data as type \code{"float"} results in 
#' a loss of precision, but cuts memory usage in half. Matrices containing whole
#' numbers only can be stored as \code{"integer"}, increasing both speed and 
#' memory efficiency. Memory usage can be reduced further if whole numbers are
#' between -32768 and 32767 by storing the data as type \code{"short"}, and even
#' further if all values are between -128 and 127 by storing values as type 
#' \code{"char"}. \code{typeof} will return the type of data stored within a
#' \code{'bigMatrix'}.
#' }
#' \subsection{\code{bigMatrix} backing files}{
#' \code{bigMatrix} objects are stored on disk at the path provided in the 
#' \code{backingfile} argument. There a four files associated with every
#' \code{bigMatrix}: the file with the ".bin" extension stores the binary 
#' representation of the data. The file with the ".desc" extension stores a 
#' description used by R to load the \code{bigMatrix} object. The files ending
#' with "_rownames.txt" and "_colnames.txt" store the row and column names of 
#' the \code{bigMatrix} that are loaded in by R when calling 
#' \code{load.bigMatrix}.
#' 
#' \strong{Warning:} \code{load.bigMatrix} can also be used to load in 
#' \code{\link[bigmemory]{big.matrix}} objects from the bigmemory package, but 
#' existing row and column information will be stripped out from its descriptor 
#' file. This means if you later load in the data as a \code{big.matrix} instead
#' of a \code{bigMatrix} the row and column names will be missing unless you
#' convert back using \code{\link{as.big.matrix}}.
#' }
#' \subsection{Implementation details}{
#' A \code{bigMatrix} object is simply a wrapper around a 
#' \code{\link[bigmemory]{big.matrix}} object from the 
#' \code{\link[bigmemory]{bigmemory}} package, with a few minor differences:
#' \enumerate{
#'  \item{\code{'bigMatrix'} objects must be backed by a file on disk.}
#'  \item{\code{'bigMatrix'} allows overwriting of these files by the user}
#'  \item{\code{'bigMatrix'} objects are stored as absolute file paths and
#'            are attached only as required.}
#'  \item{The row and column names are stored separately on disk and only
#'            accessible from R}
#' }
#' In particular, keeping the \code{'bigMatrix'} in a detached state makes the 
#' R session more reproducible: the data can be trivially reloaded into an R 
#' session, and R sessions can be reloaed within RStudio without causing the
#' application to crash. Storing the row and column names separately also offers
#' speed improvements for computation on the matrices in C++.
#' }
#' 
#' @slot descriptor path of the descriptor file for the big.matrix.
#' @slot matrix either the big.matrix object, or empty, depending on the value
#'  of \code{atttached}.
#' @slot rownames Optional vector of rownames.
#' @slot colnames Optional vector of colnames.
#' @slot attached Logical; \code{TRUE} when the big.matrix is attached in the
#'  R session, and \code{FALSE} otherwise. 
#' 
#' @seealso 
#'  \code{\link[bigmemory]{big.matrix}}.
#'  
#' @aliases 
#'  [,bigMatrix,ANY,ANY,ANY-method 
#'  [,bigMatrix,ANY,missing,ANY-method
#'  [,bigMatrix,missing,ANY,ANY-method
#'  [,bigMatrix,missing,missing,ANY-method
#'  [<-,bigMatrix,ANY,ANY-method
#'  [<-,bigMatrix,ANY,missing-method
#'  [<-,bigMatrix,missing,ANY-method
#'  [<-,bigMatrix,missing,missing-method
#'  dim,bigMatrix-method
#'  dimnames,bigMatrix-method
#'  dimnames<-,bigMatrix,ANY-method
#'  bigMatrix
#'  
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
#' # 'bigMatrix' objects can be manipulated the same way as regular matrices:
#' head(discovery_data)
#' test_network[1:5, 1:5]
#' discovery_data[,"Node_1"]
#' dim(test_data)
#' rownames(test_data)
#' nrow(test_network)
#' ncol(test_correlation)
#' is.bigMtrix(discovery_data)
#' typeof(discovery_data)
#' 
#' # For matrix algebra the whole matrix must be copied into memory first:
#' discovery_data <- as.matrix(discovery_data)
#' test_data <- test_data[,] # equivalent to 'as.matrix'
#' t(test_network[,])
#' 
#' }
#'
#' @import bigmemory
#' @import methods
#' @name bigMatrix
setClass("bigMatrix",
    slots=list(
      descriptor="character",
      matrix="maybe.big.matrix",
      rownames="optional.dimnames",
      colnames="optional.dimnames",
      attached="logical"
    ),
    validity=function(object) {
      errors <- character()
      if (is.na(object@attached)) {
        msg <- c("Attached must be one of TRUE or FALSE")
        errors <- c(errors, msg)
      } else if (object@attached) {
        if (class(object@matrix) == "NULL") {
          msg <- paste(
            "Invalid state, matrix slot must contain a big.matrix",
            "when attached is TRUE"
          )
          errors <- c(errors, msg)
          return(errors)
        }
        if (!is.null(object@rownames)) {
          if (length(object@rownames) != nrow(object@matrix)) {
            msg <- paste0(
              "The number of rownames provided (", length(object@rownames), 
              ") does not match the number of rows in the big.matrix object (",
              nrow(object@matrix), ")."
            )
            errors <- c(errors, msg)
          }
        } 

        if (!is.null(object@colnames)) {
          if (length(object@colnames) != ncol(object@matrix)) {
            msg <- paste0(
              "The number of colnames provided (", length(object@colnames), 
              ") does not match the number of columns in the big.matrix object (",
              nrow(object@matrix), ")."
            )
            errors <- c(errors, msg)
          }
        } 
        
      } else {
        if (class(object@matrix) == "big.matrix") {
          msg <- c(
            "Invalid state, matrix slot must be empty when attached is FALSE"
          )
          errors <- c(errors, msg)
        }
        tmp <- bigmemory::attach.big.matrix(object@descriptor)
        if (!is.null(object@rownames)) {
          if (length(object@rownames) != nrow(tmp)) {
            msg <- paste0(
              "The number of rownames provided (", length(object@rownames), 
              ") does not match the number of rows in the big.matrix object (",
              nrow(tmp), ")."
            )
            errors <- c(errors, msg)
          }
        } 
        
        if (!is.null(object@colnames)) {
          if (length(object@colnames) != ncol(tmp)) {
            msg <- paste0(
              "The number of colnames provided (", length(object@colnames), 
              ") does not match the number of columns in the big.matrix object (",
              nrow(tmp), ")."
            )
            errors <- c(errors, msg)
          }
        }  
      }
      if (length(errors) == 0) {
        return(TRUE) 
      } else {
        return(errors)
      }
    }
  )
