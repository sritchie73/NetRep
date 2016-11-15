#' Load a 'bigMatrix' (deprecated)
#' 
#' The \code{'bigMatrix'} class is no longer implemented in the \code{NetRep}
#' package: the shared memory approach was incompatabile with high performance
#' compute clusters, so the parallel permutation procedure has been translated
#' into C++ code (which is also much faster). The \code{\link{disk.matrix}}
#' class should now be used instead when analysing large datasets.
#' 
#' @details
#' This function will convert \code{'bigMatrix'} data saved by previous 
#' versions of \pkg{NetRep} to a serialized R matrix saved in the same location
#' and return a \code{\link{disk.matrix}} object with the associated file path.
#' If this conversion has taken place already the function will throw a warning.
#' 
#' This function will also convert the \code{'bigMatrix'} descriptor file to a 
#' \code{\link[bigmemory]{big.matrix}} descriptor file to preserve 
#' compatability with functions in the \pkg{bigmemory} package. If this 
#' functionality is not required, the files with the extensions ".bin" and 
#' ".desc" may be removed.
#' 
#' A note for users using multi-node high performance clusters:
#' \code{'big.matrix'} objects are not suitable for general usage. Access
#' to file-backed shared memory segments on multi-node systems is very slow
#' due to consistency checks performed by the operating system. This becomes
#' exponentially worse the more R sessions there are simultaneously accessing
#' the shared memory segment, e.g. through parallel \pkg{foreach} loops.
#' 
#' @param backingfile path to the backingfile for the \code{'bigMatrix'}. The
#'   file extension must be omitted.
#'   
#' @rdname bigMatrix
#' @keywords internal
#' @export
load.bigMatrix <- function(backingfile) {
  if (!pkgReqCheck("bigmemory")) {
    stop("the 'bigmemory' package must be installed")
  }
  warning("The 'bigMatrix' class is deprecated", immediate.=TRUE)
  
  # Get components for interfacing with bigmemory and resolve paths as absolute
  backingname <- basename(backingfile)
  backingpath <- gsub(paste0(backingname, "$"), "", backingfile)
  if (backingpath == "") {
    backingpath <- "."
  }
  backingpath <- normalizePath(backingpath)
  backingfile <- file.path(backingpath, backingname)
  
  # Check for row and column names
  rn <- NULL
  rnFile <- paste0(backingfile, "_rownames.txt")
  if (file.exists(rnFile))
    rn <- as.character(read.table(rnFile, stringsAsFactors=FALSE)[,1])
  
  cn <- NULL
  cnFile <- paste0(backingfile, "_colnames.txt")
  if (file.exists(cnFile))
    cn <- as.character(read.table(cnFile, stringsAsFactors=FALSE)[,1])
  
  # Check if the file is already a 'big.matrix', and handle appropriately
  descFile <- paste0(backingfile, ".desc")
  d1 <- dget(descFile)
  if (!is.null(cn))
    d1@description$colNames <- cn
  if (!is.null(rn))
    d1@description$rowNames <- rn
  
  if (!is.null(cn) || !is.null(rn)) {
    message("Converting descriptor file to the 'big.matrix' format from the",
            " 'bigmemory' package.")
    dput(d1, file=descFile)
    unlink(rnFile, force = TRUE)
    unlink(cnFile, force = TRUE)
  }
  
  serialfile <- paste0(backingfile, ".rds")
  if (!file.exists(serialfile)) {
    message("Serialising matrix to be compatabile with the new 'disk.matrix'",
            " class")
    bm <- bigmemory::attach.big.matrix(descFile)
    dm <- as.disk.matrix(bm, serialfile)
    gc()
    message("Returning 'disk.matrix' object")
    return(dm)
  }
  stop("'bigMatrix' data already converted to 'disk.matrix': use 'attach.matrix'")
}
