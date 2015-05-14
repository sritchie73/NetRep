# Allows the matrix slot of the big.matrix class to be empty, which
# is desirable if we're forcing netrep to consume minimal memory.
setClassUnion("maybe.big.matrix", c("big.matrix", "NULL"))

# Allows the rownames and colnames slots of the big.matrix class to be empty.
setClassUnion("optional.dimnames", c("character", "NULL"))

#' The bigMatrix class
#' 
#' A \code{bigMatrix} is an object which points to a data matrix stored on disk,
#' which is accessible via shared memory. This has a number of advantages, 
#' including instananeous load times and shared access to the same object across
#' multiple R sessions, enabling massively parallel operations on a matrix 
#' without the burden of excessive memory consumption.
#' 
#' @details
#' A \code{bigMatrix} object is simply a wrapper around a 
#' \code{\link[bigmemory]{big.matrix}} object from the 
#' \code{\link[bigmemory]{bigmemory}} package. 
#' 
#' When using a \code{big.matrix} object, R will eventually load the full matrix
#' into RAM when performing permutation testing due to the random sampling. The
#' \code{bigMatrix} wrapper allows the \code{big.matrix} to exist in either
#' attached or detached states, allowing \code{netRep} to optimise memory usage,
#' which is advantageous when parallelising on a machine with limited memory.
#' 
#' The authors of the \code{bigmemory} package also explicitly warn against
#' storage of row and column names in a \code{big.matrix} as it can slow down 
#' computation. A \code{bigMatrix} object handles these by keeping the dimension
#' names stored in R.
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
#'  \code{\link[bigmemory]{big.matrix}},
#'  \code{\link{bigMatrix-io}}
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
          if (length(object@colnames) != nrow(object@matrix)) {
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
          if (length(object@colnames) != nrow(tmp)) {
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
