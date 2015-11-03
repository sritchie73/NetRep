#' Loading in data as a 'bigMatrix'
#' 
#' Functions for converting to, reading in, and loading a 
#' \code{\link[=bigMatrix-class]{bigMatrix}} object.
#'
#' @param x a regular matrix to convert 
#' @param file the name of the file which the (non-bigMatrix) data are to be
#'  read from (see \code{\link[bigmemory]{read.big.matrix}}).
#' @param backingfile location on disk the bigMatrix is, or will be, stored at.
#' @param type the type of the atomic element 
#'  (\code{options()$bigmemory.default.type} by default - "\code{double}" - 
#'  but can be changed by the user to "\code{integer}", "\code{short}", or 
#'  "\code{char}").
#' @param row.names logical; does the first column of the file to be read in
#'  contain row names?
#' @param header logical; does the first line of the file to be read in contain
#'  column names?
#' @param ... additional arguments to pass to 
#'  \code{\link[bigmemory]{read.big.matrix}}.
#' 
#' @return
#'   \code{read.bigMatrix}, \code{load.bigMatrix}, \code{as.bigMatrix}: an 
#'   object of class \code{\link[=bigMatrix-class]{bigMatrix}}.
#'
#' @details 
#' \code{\link[=bigMatrix-class]{bigMatrix}} objects are stored on disk using
#' across four different files per object. The \code{backingname} designates the
#' basename of these files, and the backingpath designates the directory they
#' are stored in. The matrix data are stored in the file with the ".bin"
#' extension. The file with the ".desc" extension contains a description
#' allowing R to load in this data (see
#' \code{\link[bigmemory]{attach.big.matrix}}). Row and column names are stored
#' separately in the files ending with "_rownames.txt" and "_colnames.txt".
#' 
#' \code{save.as.bigMatrix} takes a standard matrix and writes it out as a 
#' \code{bigMatrix}. This can then be loaded into an R session using 
#' \code{load.bigMatrix}. \code{as.bigMatrix} combines these two steps, taking
#' a standard matrix and returning a \code{bigMatrix} object. Finally, 
#' \code{read.bigMatrix} can be used to read in matrix data from a file (i.e. 
#' data recongisiable by \code{read.table}) as a \code{bigMatrix} object.
#' 
#' \code{load.bigMatrix} can also be used to load in a 
#' \code{\link[bigmemory]{big.matrix}} as a 'bigMarix', but existing row and 
#' column information will be stripped out from the existing backing file. This
#' means if you later load in the data as a 'big.matrix' instead of a
#' 'bigMatrix' the row and column names will be missing unless you convert back
#' using \code{\link{as.big.matrix}}.
#' 
#' @name bigMatrix-get
#' @seealso 
#'  \code{\link[=bigMatrix-class]{bigMatrix}},
#'  \code{\link[bigmemory]{big.matrix}},
#'  \code{\link[bigmemory]{read.big.matrix}}
#' @export
#' @import bigmemory
save.as.bigMatrix <- function(
  x, backingfile, type=options("bigmemory.default.type")[[1]]
) {
  if (!(class(x) == "matrix"))
    stop("'x' must be a 'matrix'")
  
  # Get components for interfacing with bigmemory and resolve paths as absolute
  backingname <- basename(backingfile)
  backingpath <- gsub(paste0(backingname, "$"), "", backingfile)
  if (backingpath == "") {
    backingpath <- "."
  }
  backingpath <- normalizePath(backingpath)
  backingfile <- file.path(backingpath, backingname)
  
  # Dimension names are saved separately from the big.matrix object.
  cnFile <- paste0(backingfile, "_colnames.txt")
  if (!is.null(colnames(x))) {
    write.table(
      colnames(x), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE,
      file=cnFile
    )
  } else {
    if (file.exists(cnFile))
      unlink(cnFile)
  }
  rnFile <- paste0(backingfile,  "_rownames.txt")
  if (!is.null(rownames(x))) {
    write.table(
      rownames(x), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE,
      file=rnFile
    ) 
  } else {
    if (file.exists(rnFile))
      unlink(rnFile)
  }
  
  dimnames(x) <- NULL
  
  # now save the big.matrix to disk
  bigmemory::as.big.matrix(
    x, backingpath=backingpath, type=type,
    backingfile=paste0(backingname, ".bin"),
    descriptorfile=paste0(backingname, ".desc")
  )
  
  invisible(NULL)
}

#' @rdname bigMatrix-get
#' @export
load.bigMatrix <- function(
  backingfile
) {
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
      warning(
        "row and column names will be stripped from existing 'big.matrix'"
      )
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

#' @rdname bigMatrix-get
#' 
#' @export
as.bigMatrix <- function(
  x, backingfile, type=options("bigmemory.default.type")[[1]]
) {
  if (class(x) == "big.matrix") 
    stop(
      "use load.bigMatrix to load in an existing 'big.matrix' as a 'bigMatrix"
    )
  if (class(x) != "matrix")
    stop("Cannot convert from ", class(x), " to 'bigMatrix'")
  
  save.as.bigMatrix(x, backingfile, type)
  load.bigMatrix(backingfile)
}

#' @rdname bigMatrix-get
#' @export
read.bigMatrix <- function(
  file, backingfile, 
  type=options("bigmemory.default.type")[[1]],
  row.names=TRUE, header=TRUE, ...
) {
  # Get components for interfacing with bigmemory and resolve paths as absolute
  backingname <- basename(backingfile)
  backingpath <- gsub(paste0(backingname, "$"), "", backingfile)
  if (backingpath == "") {
    backingpath <- "."
  }
  backingpath <- normalizePath(backingpath)
  backingfile <- file.path(backingpath, backingname)
  
  bm <- read.big.matrix(
    filename=file, 
    backingfile=paste0(backingname, ".bin"),
    descriptorfile=paste0(backingname, ".desc"),
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


#' Writing out data, or converting from, a 'bigMatrix'.
#' 
#' Functions for converting 'bigMatrix' objects to other classes.
#' 
#' @param x an object of type 'bigMatrix'.
#' @param file file to write out the data from \code{x} into.
#' @param ... additional arguments to pass to 
#'  \code{\link[bigmemory]{write.big.matrix}}.
#' 
#' @details
#' \code{as.big.matrix} converts to a \code{\link[bigmemory]{big.matrix}}
#' object.
#' \code{as.matrix} converts to a \code{\link[base]{matrix}} object.
#' \code{write.bigMatrix} will write out the data into a regularly structure
#' file (see \code{\link[utils]{write.table}})
#' 
#' @name bigMatrix-out
#' @export
write.bigMatrix <- function(x, file, ...) {
  write.table(x=x[,], file, ...)
}

#' @rdname bigMatrix-out
setMethod("as.matrix", signature(x="bigMatrix"), function(x) {
  x[,]
})

#' @rdname bigMatrix-out
setMethod(
  "as.big.matrix", signature(x="bigMatrix"), function(x) {
    if (!file.exists(x@descriptor))
      stop("Could not find backing file. Have you changed working directory?")
    rnFile <- gsub(".desc", "_rownames.txt", x@descriptor)
    cnFile <- gsub(".desc", "_colnames.txt", x@descriptor)
    
    desc <- dget(x@descriptor)
    
    if (!is.null(x@colnames))
      desc@description$colNames <- x@colnames
    if (!is.null(x@rownames))
      desc@description$rowNames <- x@rownames
    
    dput(desc, x@descriptor)
    unlink(rnFile)
    unlink(cnFile)
    
    bigmemory::attach.big.matrix(x@descriptor)
  })
