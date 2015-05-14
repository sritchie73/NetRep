#' bigMatrix input and output
#' 
#' Functions for converting to, reading in, and loading a 
#' \code{\link{bigMatrix}} object.
#'
#' @param x matrix to write out as a file-backed `big.matrix`.
#' @param file the name of the file which the (non-bigMatrix) data are to be
#'  read from using \code{\link[utils]{read.table}}.
#' @param backingpath directory the bigMatrix is, or will be, stored in.
#' @param backingname basename (filename without an extension) the
#'  bigMatrix object is, or will be, stored in.
#' @param the type of the atomic element 
#'  (\code{options()$bigmemory.default.type} by default – "\code{double}" – 
#'  but can be changed by the user to "\code{integer}", "\code{short}", or 
#'  "\code{char}").
#' @param ... additional arguments to pass to \code{read.table}
#' 
#' @return
#'   \code{read.bigMatrix}, \code{load.bigMatrix}, \code{as.bigMatrix}: an 
#'   object of class \code{\link{bigMatrix}}.
#'
#' @details 
#' \code{\link{bigMatrix}} objects are stored on disk using across four 
#' different files per object. The \code{backingname} designates the basename 
#' of these files, and the backingpath designates the directory they are stored 
#' in. The matrix data are stored in the file with the ".bin" extension. The 
#' file with the ".desc" extension contains a description allowing R to load in
#' this data (see \code{\link[bigmemory]{attach.big.matrix}}). Row and 
#' column names are stored separately in the files ending with "_rownames.txt"
#' and "_colnames.txt". 
#' 
#' \code{write.bigMatrix} takes a standard matrix and writes it out as a 
#' \code{bigMatrix}. This can then be loaded into an R session using 
#' \code{load.bigMatrix}. \code{as.bigMatrix} combines these two steps, taking
#' a standard matrix and returning a \code{bigMatrix} object. Finally, 
#' \code{read.bigMatrix} can be used to read in matrix data from a file (i.e. 
#' using \code{read.table}) as a \code{bigMatrix} object (note that the 
#' appropriate backing files will also be created on disk).
#' 
#' @name bigMatrix-io
#' @seealso 
#'  \code{\link{bigMatrix}},
#'  \code{\link[bigmemory]{big.matrix}},
#'  \code{\link[utils]{read.table}}
#' @export
#' @import bigmemory
write.bigMatrix <- function(
  x, backingname, backingpath, type=options("bigmemory.default.type")
) {
  # Take a matrix x, 
  if (class(x) != "matrix")
    stop("'x' must be a matrix.")
  
  # Dimension names are saved separately from the big.matrix object.
  if (!is.null(colnames(x))) {
    write.table(
      rownames(x), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE,
      file=file.path(backingpath, paste(backingname, "_rownames.txt"))
    )
  }
  if (!is.null(rownames(x))) {
    write.table(
      colnames(x), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE,
      file=file.path(backingpath, paste(backingname, "_colnames.txt"))
    ) 
  }
  
  # Remove dimension names before saving the big.matrix
  dimnames(x) <- NULL
  
  # now save the big.matrix to disk
  as.big.matrix(
    x, backingpath=backingpath, type=type,
    backingfile=paste0(backingname, ".bin"),
    descriptorfile=paste0(backingname, ".desc")
  )
  invisible(NULL)
}

#' @rdname bigMatrix-io
#' @export
load.bigMatrix <- function(
  backingname, backingpath
) {
  
  rn <- NULL
  rnFile <- file.path(backingpath, paste0(backingname, "_rownames.txt"))
  if (file.exists(rnFile))
    rn <- read.table(rnFile)
  
  cn <- NULL
  cnFile <- file.path(backingpath, paste0(backingname, "_colnames.txt"))
  if (file.exists(cnFile))
    rn <- read.table(cnFile)

  new("bigMatrix",
    descriptor=file.path(backingpath, paste0(backingname, ".desc")),
    rownames=rn,
    colnames=cn,
    attached=FALSE
  )
}

#' @rdname bigMatrix-io
#' @export
as.bigMatrix <- function(
  x, backingname, backingpath, type=options("bigmemory.default.type")
) {
  write.bigMatrix(x, backingname, backingpath, type)
  load.bigMatrix(backingname, backingpath)
}

#' @rdname bigMatrix-io
#' @export
read.bigMatrix <- function(
  file, backingname, backingpath, type=options("bigmemory.default.type"), ...
) {
  as.bigMatrix(read.table(file, ...), backingname, backingpath, type)
}
