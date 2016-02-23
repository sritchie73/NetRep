#' \code{save.as.bigMatrix} saves a numeric \code{matrix} or \code{data.frame} 
#' in shared memory as a \code{bigMatrix} at the provided \code{backingfile}
#' location. 
#' 
#' @param x \code{save.as.bigMatrix} and \code{as.bigMatrix}: a \code{matrix}
#'  or \code{data.frame} to convert to a \code{bigMatrix}. 
#'  \code{write.bigMatrix} and \code{is.bigMatrix}: a \code{bigMatrix} object. 
#' @param file the name of the file which the (non-bigMatrix) data are to be
#'  read from or written to (see \code{\link[bigmemory]{read.big.matrix}} and
#'  \code{\link[bigmemory]{write.big.matrix}}).
#' @param backingfile location on disk the bigMatrix is, or will be, stored at.
#' @param type the type of the atomic element. \code{"double"} by default,
#'  but can be changed by the user to \code{"float"} , \code{"integer"}, 
#'  \code{"short"}, or \code{"char"} (see details).
#' @param row.names logical; does the first column of the file to be read in
#'  contain row names?
#' @param header logical; does the first line of the file to be read in contain
#'  column names?
#' @param ... additional arguments to pass to 
#'  \code{\link[bigmemory]{read.big.matrix}} or 
#'  \code{\link[bigmemory]{write.big.matrix}}.
#'
#' @return
#'   \code{read.bigMatrix}, \code{load.bigMatrix}, \code{as.bigMatrix}: an 
#'   object of class \code{bigMatrix}.
#'
#' @rdname bigMatrix
#'  
#' @export
#' @import bigmemory
save.as.bigMatrix <- function(
  x, backingfile, type="double"
) {
  if (!(class(x) %in% c("matrix", "data.frame")))
    stop("'x' must be a 'matrix' or 'data.frame'")
  
  if (class(x) == "data.frame") {
    x <- as.matrix(x)
    if (typeof(x) %nin% c("numeric", "integer", "double"))
      stop("could not convert 'x' to a numeric matrix")
  }
  
  if (missing(backingfile)) {
    backingfile <- file.path(tempdir(), getUUID())
    message("no 'backingfile' provided, saving to", backingfile)
  }
    
  # Get components for interfacing with bigmemory and resolve paths as absolute
  backingname <- basename(backingfile)
  backingpath <- gsub(paste0(backingname, "$"), "", backingfile)
  if (backingpath == "") {
    backingpath <- "."
  }
  backingpath <- normalizePath(backingpath)
  backingfile <- file.path(backingpath, backingname)
  
  binFile <- paste0(backingname, ".bin")
  descFile <- paste0(backingname, ".desc")
  
  fullBinFile <- paste0(backingfile, ".bin")
  fullDescFile <- paste0(backingfile, ".bin")
  cnFile <- paste0(backingfile, "_colnames.txt")
  rnFile <- paste0(backingfile,  "_rownames.txt")
  
  # Allow overwriting
  if (file.exists(fullBinFile)) {
    if (!file.remove(fullBinFile)) {
      stop("cannot remove 'bigMatrix' backingfile. Object must first be removed",
           " in R and then the garbage collector ('gc') must be called")
    }
    file.remove(fullDescFile)
    file.remove(cnFile)
    file.remove(rnFile)
    
  }
  
  # Dimension names are saved separately from the big.matrix object.
  if (!is.null(colnames(x))) {
    write.table(
      colnames(x), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE,
      file=cnFile
    )
  }
  if (!is.null(rownames(x))) {
    write.table(
      rownames(x), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE,
      file=rnFile
    ) 
  }
  dimnames(x) <- NULL
  
  # now save the big.matrix to disk
  bigmemory::as.big.matrix(
    x, backingpath=backingpath, type=type,
    backingfile=binFile,
    descriptorfile=descFile
  )
  
  invisible(NULL)
}

#' \code{as.bigMatrix} converts a numeric \code{matrix} or \code{data.frame}
#' into a \code{bigMatrix} object and stores it on disk at the provided 
#' \code{backingfile} location. 
#'  
#' @rdname bigMatrix
#' @export
as.bigMatrix <- function(
  x, backingfile, type="double"
) {
  if (class(x) == "big.matrix") 
    stop(
      "use load.bigMatrix to load in an existing 'big.matrix' as a 'bigMatrix"
    )
  if (!(class(x) %in% c("matrix", "data.frame")))
    stop("'x' must be a 'matrix' or 'data.frame'")
  
  if (class(x) == "data.frame") {
    x <- as.matrix(x)
    if (typeof(x) %nin% c("numeric", "integer", "double"))
      stop("could not convert 'x' to a numeric matrix")
  }
  
  if (missing(backingfile)) {
    backingfile <- file.path(tempdir(), getUUID())
    message("no 'backingfile' provided, saving to", backingfile)
  }
  
  save.as.bigMatrix(x, backingfile, type)
  load.bigMatrix(backingfile)
}

#' \code{load.bigMatrix} loads a \code{bigMatrix} object stored at the provided
#' \code{backingfile} location into R.
#'  
#' @rdname bigMatrix
#' @export
load.bigMatrix <- function(backingfile) {
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
  if (file.exists(descFile)) {
    d1 <- dget(descFile)
    if (!is.null(d1@description$colNames) || 
        !is.null(d1@description$rowNames)) 
    {
      rn <- d1@description$rowNames
      cn <- d1@description$colNames
      d1@description["rowNames"] <- list(NULL)
      d1@description["colNames"] <- list(NULL)
      dput(d1, descFile)
      if (!is.null(rn)) {
        write.table(
          rn, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE,
          file=rnFile
        ) 
      }
      if (!is.null(cn)) {
        write.table(
          cn, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE,
          file=cnFile
        ) 
      }
    }
  }
  
  new("bigMatrix",
    descriptor=descFile,
    rownames=rn,
    colnames=cn,
    attached=FALSE
  )
}

#' \code{read.bigMatrix} reads a file in table format (i.e. 
#' \code{\link[utils]{read.table}}) into R, creates a \code{bigMatrix} from
#' it, and stores it at the provided \code{backingfile} location.
#' 
#' @rdname bigMatrix
#' @export
read.bigMatrix <- function(
  file, backingfile, 
  type="double",
  row.names=TRUE, header=TRUE, ...
) {
  if (missing(backingfile)) {
    backingfile <- file.path(tempdir(), getUUID())
    message("no 'backingfile' provided, saving to", backingfile)
  }
  
  # Get components for interfacing with bigmemory and resolve paths as absolute
  backingname <- basename(backingfile)
  backingpath <- gsub(paste0(backingname, "$"), "", backingfile)
  if (backingpath == "") {
    backingpath <- "."
  }
  backingpath <- normalizePath(backingpath)
  backingfile <- file.path(backingpath, backingname)
  
  binFile <- paste0(backingname, ".bin")
  descFile <- paste0(backingname, ".desc")
  
  fullBinFile <- paste0(backingfile, ".bin")
  fullDescFile <- paste0(backingfile, ".bin")
  cnFile <- paste0(backingfile, "_colnames.txt")
  rnFile <- paste0(backingfile,  "_rownames.txt")
  
  # Allow overwriting
  if (file.exists(fullBinFile)) {
    if (!file.remove(fullBinFile)) {
      stop("cannot remove 'bigMatrix' backingfile. Object must first be removed",
           " in R and then the garbage collector ('gc') must be called")
    }
    file.remove(fullDescFile)
    file.remove(cnFile)
    file.remove(rnFile)
  }
  
  bm <- read.big.matrix(
    filename=file, 
    backingfile=binFile,
    descriptorfile=descFile,
    backingpath=backingpath, type=type,
    has.row.names=row.names, header=header,
    ...
  )
  # For some reason, read.big.matrix doesn't add the row names to the descriptor
  # file when both rownames and colnames are present (in some cases). So we 
  # manually override that.
  dFile <- paste0(backingfile, ".desc")
  desc <- dget(dFile)
  if (!is.null(rownames(bm))) {
    desc@description$rowNames <- rownames(bm)
  } else {
    desc@description["rowNames"] <- list(NULL)
  }
  if (!is.null(colnames(bm))) {
    desc@description$colNames <- colnames(bm)
  } else {
    desc@description["colNames"] <- list(NULL)
  }
  dput(desc, dFile)
  rm(bm)
  gc()
  # load.bigMatrix will throw a warning when converting from a big.matrix 
  # object to a bigMatrix object. Since we're just using this function to avoid
  # code duplication, we can suppress this warning. Any other warnings will be
  # thrown above when using bigmemory::read.big.matrix
  suppressWarnings(load.bigMatrix(backingfile)) 
}

#' \code{write.bigMatrix} writes a \code{bigMatrix} object out as a file in 
#' table format (i.e. \code{\link[utils]{write.table}}).
#'
#' @return
#'   \code{save.as.bigMatrix}, \code{write.bigMatrix}: \code{NULL}, but create
#'   files on disk as a side effect.
#'
#' @rdname bigMatrix
#' @export
write.bigMatrix <- function(x, file, ...) {
  write.table(x=x[,,drop=FALSE], file, ...)
}

#' @return \code{is.bigMatrix}: \code{TRUE} if \code{x} is a \code{bigMatrix}.
#'
#' @rdname bigMatrix
#' @export
is.bigMatrix <- function(x) {
  return(class(x) == "bigMatrix")
}
