#' Load a 'bigMatrix' (deprecated)
#' 
#' The \code{'bigMatrix'} class is no longer implemented in the \code{NetRep}
#' package: the shared memory approach was incompatabile with high performance
#' compute clusters, so the parallel permutation procedure has been translated
#' into C++ code (which is also much faster). For those who had saved their 
#' data in the \code{'bigMatrix'} format, this function will convert these 
#' objects to the similar \code{'\link[bigmemory]{big.matrix}'} class
#' from the \code{\link{bigmemory}} package.
#' 
#' @rdname bigMatrix
#' @importFrom bigmemory attach.big.matrix
#' @export
load.bigMatrix <- function(backingfile) {
  message("The 'bigMatrix' class is deprecated. Converting to file to a ",
          "file-backed 'big.matrix' object (see the bigmemory package).")
  
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
  
  dput(d1, file=descFile)
  unlink(rnFile, force = TRUE)
  unlink(cnFile, force = TRUE)
  
  bigmemory::attach.big.matrix(descFile)
}
