#' The 'disk.matrix' class
#' 
#' A \code{'disk.matrix'} contains a file path to a matrix stored on disk,
#' along with meta data for how to read that file. This allows \pkg{NetRep}
#' to load datasets into RAM only when required, i.e. one at a time. This 
#' significantly reduces the memory usage of R when analysing large datasets.
#' \code{'disk.matrix'} objects may be supplied instead of \code{'matrix'} 
#' objects in the input list arguments \code{'network'}, \code{'data'}, and 
#' \code{'correlation'}, which are common to most of \pkg{NetRep}'s functions.
#' 
#' @param file for \code{attach.disk.matrix} the file name of a matrix on disk. 
#'  For \code{as.disk.matrix} the file name to save the matrix to. For 
#'  \code{serialize.table} the file name of a matrix in table format on disk.
#' @param serialized determines how the matrix will be loaded from disk into R
#'  by \code{as.matrix}. If \code{TRUE}, the \code{readRDS} function 
#'  will be used. If \code{FALSE}, the \code{read.table} function will 
#'  be used. 
#' @param x for \code{as.matrix} a \code{disk.matrix} object to load into R. 
#'  For \code{as.disk.matrix} an object to convert to a \code{disk.matrix}. For 
#'  \code{is.disk.matrix} an object to check if its a \code{disk.matrix}.
#' @param object a \code{'disk.matrix'} object. 
#' @param serialize determines how the matrix is saved to disk by 
#'  \code{as.disk.matrix}. If \code{TRUE} it will be stored as a serialized R 
#'  object using \code{saveRDS}. If \code{FALSE} it will be stored as a 
#'  tab-separated file using \code{write.table}.
#' @param ... arguments to be used by \code{read.table} when reading in matrix 
#'  data from a file in table format.
#' 
#' @details
#' Matrices may either be stored as regular table files that can be read by
#' \code{\link{read.table}}, or as serialized R objects that can be read by
#' \code{\link{readRDS}}. Serialized objects are much faster to load, but 
#' cannot be read by other programs. 
#' 
#' The \code{attach.disk.matrix} function creates a \code{disk.matrix} object
#' from a file path. The \code{as.matrix} function will load the data from disk
#' into the R session as a regular \code{\link{matrix}} object.
#' 
#' The \code{as.disk.matrix} function converts a matrix into a 
#' \code{disk.matrix} by saving its contents to the specified \code{file}. The
#' \code{serialize} argument determines whether the data is stored as a 
#' serialized R object or as a tab-separated file (i.e. \code{sep="\\t"}). We
#' recommend storing the matrix as a serialized R object unless disk space is
#' a concern. More control over the storage format can be obtained by using
#' \code{saveRDS} or \code{write.table} directly.
#' 
#' The \code{serialize.matrix} function converts a file in table format to a
#' serialized R object with the same file name, but with the ".rds" extension.
#' 
#' @section Warning:
#' \code{attach.disk.matrix} does not check whether the specified file can be
#' read into R. \code{as.matrix} will fail and throw an error if this is the
#' case.
#' 
#' @return 
#' A \code{disk.matrix} object (\code{attach.disk.matrix}, \code{as.disk.matrix}),
#' a \code{matrix} (\code{as.matrix}), the file path to a serialized matrix
#' (\code{serialize.table}), or a \code{TRUE} or \code{FALSE} indicating 
#' whether an object is a \code{disk.matrix} (\code{is.disk.matrix}).
#' 
#' @slot file the name of the file where the matrix is saved.
#' @slot read.func either \code{"read.table"} or \code{"readRDS"}.
#' @slot func.args a list of arguments to be supplied to the \code{'read.func'}.
#'
#' @import methods
#' @name disk.matrix
setClass("disk.matrix", 
  slots=list(
    file="character",
    read.func="character",
    func.args="list"
  ),
  validity=function(object) {
    errors <- character()
    if (length(object@file) != 1) {
      msg <- "slot 'file' must be a single value"
      errors <- c(errors, msg)
    }
    if (!file.exists(object@file)) {
      msg <- "slot 'file' must be a file path to an existing file"
      errors <- c(errors, msg)
    }
    if (length(object@read.func) != 1 && object@read.func %nin% c("read.table", "readRDS")) {
      msg <- "slot 'read.func' must be either \"read.table\" or \"readRDS\""
      errors <- c(errors, msg)
    }
    if (length(errors) > 0) { 
      return(errors) 
    } else {
      return(TRUE)
    }
  })

#' @rdname disk.matrix
#' @export
attach.disk.matrix <- function(file, serialized=TRUE, ...) {
  if (is.na(serialized) || length(serialized) != 1) {
    stop("'serialized' must be 'TRUE' or 'FALSE'")
  }
  if (length(file) != 1 || !is.character(file) || !file.exists(file)) {
    stop("'file' must be the name of a file and that file must already exist")
  }
  
  read.func <- ifelse(serialized, "readRDS", "read.table")
  new("disk.matrix", file=normalizePath(file), read.func=read.func, 
      func.args=list(...))
}

#' @rdname disk.matrix
#' @export
serialize.table <- function(file, ...) {
  if (length(file) != 1 || !is.character(file) || !file.exists(file)) {
    stop("'file' must be the name of a file and that file must already exist")
  }
  
  ext <- gsub(".*\\.", "", file)
  serialized.file <- gsub(paste0(ext, "$"), "rds", file)
  
  m <- as.matrix(read.table(file, ...))
  saveRDS(m, serialized.file)
  serialized.file
}

#' @rdname disk.matrix
#' @export
is.disk.matrix <- function(x) {
  "disk.matrix" %in% class(x)
}

#' @rdname disk.matrix
#' @export
setGeneric("as.disk.matrix", function(x, file, serialize=TRUE) {
  standardGeneric("as.disk.matrix")
})

#' @rdname disk.matrix
setMethod("as.disk.matrix", signature(x="disk.matrix"), 
          function(x, file, serialize=TRUE) {
            warning("already a 'disk.matrix'")
            return(x)
          })

#' @rdname disk.matrix
setMethod("as.disk.matrix", signature(x="matrix"), 
          function(x, file, serialize=TRUE) {
            if (is.na(serialize) || length(serialize) != 1) {
              stop("'serialize' must be 'TRUE' or 'FALSE'")
            }
            if (length(file) != 1 || !is.character(file)) {
              stop("'file' must be the name of a file to save the matrix to")
            }
            
            if (serialize) {
              saveRDS(x, file)
              attach.disk.matrix(file)
            } else {
              write.table(x, file, col.names=!is.null(colnames(x)),
                          row.names=!is.null(rownames(x)), sep="\t",
                          quote=FALSE)
              attach.disk.matrix(file, FALSE, header=!is.null(colnames(x)),
                            row.names=ifelse(is.null(rownames(x)), FALSE, 1),
                            sep="\t")
            }
          })

#' @rdname disk.matrix
setMethod("as.disk.matrix", signature(x="ANY"), 
          function(x, file, serialize=TRUE) {
            x <- as.matrix(x)
            as.disk.matrix(x, file, serialize)
          })

#' @rdname disk.matrix 
#' @export
setMethod("as.matrix", signature(x="disk.matrix"), function(x) {
  if (!file.exists(x@file)) {
    stop("file ", prettyPath(x@file), " does not exist")
  }
  as.matrix(do.call(x@read.func, c(file=x@file, x@func.args)))
})

#' @rdname disk.matrix
#' @export
setMethod("show", signature(object="disk.matrix"), function(object) {
  cat("Pointer to matrix stored at", prettyPath(object@file), "\n")
})
