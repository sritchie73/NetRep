#' Attaching and detaching a bigMatrix
#' 
#' In the attached state, R is aware of the \code{\link[bigmemory]{big.matrix}}
#' pointer object. In the detached state, R is only aware of the 
#' location of the 'big.matrix' descriptor file, and the big.matrix is attached
#' and removed every time its elements are requested.
#' 
#' @param x the bigMatrix object to attach or detach
#' 
#' @return
#'  The bigMatrix object with its state updated.
#'  
#' @name bigMatrix-state
attach.bigMatrix <- function(x) {
  if (class(x) != 'bigMatrix')
    stop("object is not of class 'bigMatrix'")
  if (x@attached)
    return(x)
  if (!file.exists(x@descriptor))
    stop("Could not find backing file. Have you changed working directory?")
  x@matrix <- bigmemory::attach.big.matrix(x@descriptor)
  poke(x@matrix)
  x@attached <- TRUE
  x
}

#' @rdname bigMatrix-state
detach.bigMatrix <- function(x) {
  if (class(x) != 'bigMatrix')
    stop("object is not of class 'bigMatrix'")
  if (!x@attached)
    return(x)
  x@matrix <- NULL
  x@attached <- FALSE
  gc()
  x
}

#' Checks if an object is a bigMatrix
#' 
#' @param x object to check
#' 
#' @return
#'  \code{TRUE} if 'x' is a 'bigMatrix', \code{FALSE} otherwise.
#'  
#' @seealso \code{\link[=bigMatrix-class]{bigMatrix}}
#' @export
is.bigMatrix <- function(x) {
  return(class(x) == "bigMatrix")
}

#-------------------------------------------------------------------------------
# Methods for standard functions, e.g. dim, nrow, ncol, etc.
# 
# With the exception of 'show', these simply attach the big.matrix object and
# call the appropriate function on the matrix object.
#-------------------------------------------------------------------------------
setMethod("show", signature(object = "bigMatrix"), function(object) {
  cat(
    '"bigMatrix" of type "', typeof(object), '" with ', nrow(object), 
    " rows and ", ncol(object), ' columns stored at "',  
    gsub(".desc", ".bin", object@descriptor), '"\n', sep=""
  )
})

setMethod("print", signature(x = "bigMatrix"), function(x) {
  show(x)
})

setMethod("dim", signature(x = "bigMatrix"), function(x) {
  # Attach the big.matrix object if not attached yet
  is.attached <- x@attached
  if (!is.attached)
    x <- attach.bigMatrix(x)
  
  res <- dim(x@matrix)
  
  # detach big.matrix object if it was detached to begin with
  if (!is.attached)
    x <- detach.bigMatrix(x)
  
  res
})

setMethod("nrow", signature(x = "bigMatrix"), function(x) {
  # Attach the big.matrix object if not attached yet
  is.attached <- x@attached
  if (!is.attached)
    x <- attach.bigMatrix(x)
  
  res <- nrow(x@matrix)
  
  # detach big.matrix object if it was detached to begin with
  if (!is.attached)
    x <- detach.bigMatrix(x)
  
  res
})

setMethod("ncol", signature(x = "bigMatrix"), function(x) {
  # Attach the big.matrix object if not attached yet
  is.attached <- x@attached
  if (!is.attached)
    x <- attach.bigMatrix(x)
  
  res <- ncol(x@matrix)
  
  # detach big.matrix object if it was detached to begin with
  if (!is.attached)
    x <- detach.bigMatrix(x)
  
  res
})

setMethod("typeof", signature(x = "bigMatrix"), function(x) {
  # Attach the big.matrix object if not attached yet
  is.attached <- x@attached
  if (!is.attached)
    x <- attach.bigMatrix(x)
  
  res <- typeof(x@matrix)
  
  # detach big.matrix object if it was detached to begin with
  if (!is.attached)
    x <- detach.bigMatrix(x)
  
  res
})

setMethod("head", signature(x = "bigMatrix"), function(x, n=6){
  n <- min(as.integer(n), nrow(x))
  if (n < 1 | n > nrow(x)) 
    stop("n must be between 1 and nrow(x)")
  x[1:n,]
})

setMethod("tail", signature(x = "bigMatrix"), function(x, n=6){
  n <- min(as.integer(n), nrow(x))
  if (n < 1 | n > nrow(x)) 
    stop("n must be between 1 and nrow(x)")
  x[(nrow(x) - n + 1):nrow(x),]
})

setMethod("dimnames", signature(x = "bigMatrix"), function(x) {
  list(x@rownames, x@colnames)
})

setMethod("dimnames<-", 
  signature(x = "bigMatrix", value="ANY"), function(x, value) {
    if (!is.null(value) | !is.list(value))
      stop("'dimnames' must be a list")
    if (is.null(value)) {
      dn <- list(NULL, NULL)
    } else if (length(value) > 2) {
      stop(
        "length of 'dimnames' [", length(value), 
        "] must match that of 'dims' [2]"
      )
    } else if (length(value) == 0) {
      dn <- list(NULL, NULL)
    } else if (length(value) == 1) {
      dn <- list(as.character(unlist(value[[1L]])), NULL)
    } else {
      dn <- list(
        as.character(unlist(value[[1L]])), 
        as.character(unlist(value[[2L]]))
      )
    }
    if (length(dn[[1L]]) == 0) {
      dn[1L] <- list(NULL)
    }
    if (length(dn[[2L]]) == 0) {
      dn[2L] <- list(NULL)
    }
    
    if (!is.null(dn[[1L]]) & length(dn[[1L]]) != nrow(x)) {
      stop(
        "length of 'dimnames' [", length(dn[[1L]]), 
        "] not equal to array extent"
      )
    }
    if (!is.null(dn[[2L]]) & length(dn[[2L]]) != ncol(x)) {
      stop(
        "length of 'dimnames' [", length(dn[[2L]]), 
        "] not equal to array extent"
      )
    }
    
    # We need to update the backing files. Make sure this is successful before
    # updating the matrix wrapper!
    if (!file.exists(x@descriptor))
      stop("Could not find backing file. Have you changed working directory?")
    rnFile <- gsub(".desc", "_rownames.txt", x@descriptor)
    cnFile <- gsub(".desc", "_colnames.txt", x@descriptor)
    
    is.attached <- x@attached

    if (is.null(dn[[1L]])) {
      unlink(rnFile)
      x@rownames <- NULL
    } else {
      write.table(
        dn[[1L]], quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE,
        file=rnFile
      ) 
      x@rownames <- dn[[1L]]
    }
    if (is.null(dn[[2L]])) {
      unlink(cnFile)
      x@colnames <- NULL
    } else {
      write.table(
        dn[[2L]], quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE,
        file=cnFile
      ) 
      x@colnames <- dn[[2L]]
    }
    x
})


#-------------------------------------------------------------------------------
# Matrix subsetting methods
#-------------------------------------------------------------------------------
setMethod(
  "[", signature(x = "bigMatrix", i="ANY", j="ANY", drop="ANY"), 
  function(x, i, j, drop){

    # Check if subscript type is valid
    if (!(class(i) %in% c("character", "integer", "numeric", "logical")))
      stop("invalid subscript type '", class(i), "'")
    if (!(class(j) %in% c("character", "integer", "numeric", "logical")))
      stop("invalid subscript type '", class(j), "'")   
    
    # Handle indexing by row or column names
    if (class(i) == "character")
      i <- match(i, x@rownames)
    if (class(j) == "character")
      j <- match(j, x@colnames)
    if(any(is.na(i)) || any(is.na(j)))
      stop("subscript out of bounds")
    
    # The default behaviour is for drop to be TRUE
    if (is.null(drop) || is.na(drop) || class(drop) != "logical")
      drop <- TRUE
    
    # Attach the big.matrix object if not attached yet
    is.attached <- x@attached
    if (!is.attached)
      x <- attach.bigMatrix(x)
     
    # Fetch the results
    ret <- x@matrix[i, j, drop=drop]
    
    # detach big.matrix object if it was detached to begin with
    if (!is.attached)
      x <- detach.bigMatrix(x)
    
    # Give the results the approriate row and column names
    if (is.null(dim(ret))) {
      # Case: length(ret) == 1, we take EITHER the row or column
      # names if only one exists. Otherwise, no names
      if (length(ret) == 1) {
        if (!is.null(x@colnames) & is.null(x@rownames)) {
          names(ret) <- x@colnames[j]
        } else if (is.null(x@colnames) & !is.null(x@rownames)) {
          names(ret) <- x@rownames[i]
        }
      }
      # Case: length(i) == 1: number returned rows == 1.
      else if (length(i) == 1) {
        names(ret) <- x@colnames[j]
      }
      # Case: length(j) > 1, i is logical and j is not: number of returned 
      # rows == 1 because sum(i) == 1
      else if (is.logical(i) & !is.logical(j) & length(j) > 1) {
        names(ret) <- x@colnames[j]
      }
      # Case: length(j) == 1: number returned columns == 1.
      else if (length(j) == 1) {
        names(ret) <- x@rownames[i]
      }
      # Case: length(i) > 1, j is logical and i is not: number of returned 
      # columns == 1 because sum(j) == 1
      else if (is.logical(j) & !is.logical(i) & length(i) > 1) {
        names(ret) <- x@colnames[j]
      }
      else if (is.logical(j) & is.logical(i)) {
        nr <- sum(i)
        nc <- sum(j)
        if (nr == 1) {
          names(ret) <- x@colnames[j]
        } else if (nc == 1) {
          names(ret) <- x@rownames[i]
        }
      }
    } else {
      if (!is.null(x@rownames) & !is.null(x@colnames)) {
        dimnames(ret) <- list(x@rownames[i], x@colnames[j])
      } else if (!is.null(x@rownames)) {
        rownames(ret) <- x@rownames[i]
      } else if (!is.null(x@colnames)) {
        colnames(ret) <- x@colnames[j]
      }
    }
    ret
})

setMethod(
  "[", signature(x = "bigMatrix", i="ANY", j="missing", drop="ANY"), 
  function(x, i, j, drop){
    
    # Check if subscript type is valid
    if (!(class(i) %in% c("character", "integer", "numeric", "logical")))
      stop("invalid subscript type '", class(i), "'")

    # Handle indexing by row or column names
    if (class(i) == "character")
      i <- match(i, x@rownames)
    if(any(is.na(i)))
      stop("subscript out of bounds")
    
    # The default behaviour is for drop to be TRUE
    if (is.null(drop) || is.na(drop) || class(drop) != "logical")
      drop <- TRUE
    
    # Attach the big.matrix object if not attached yet
    is.attached <- x@attached
    if (!is.attached)
      x <- attach.bigMatrix(x)
    
    # Fetch the results
    ret <- x@matrix[i, , drop=drop]
    
    # detach big.matrix object if it was detached to begin with
    if (!is.attached)
      x <- detach.bigMatrix(x)
    
    # Give the results the approriate row and column names
    if (is.null(dim(ret)) & length(ret) != 0) {
      # Case: 1 column, length(ret) == 1, we take EITHER the row or column
      # names if only one exists. Otherwise, no names
      if (ncol(x) == 1 & length(ret) == 1) {
        if (!is.null(x@colnames) & is.null(x@rownames)) {
          names(ret) <- x@colnames
        } else if (is.null(x@colnames) & !is.null(x@rownames)) {
          names(ret) <- x@rownames[i]
        }
      }
      # Case: 1 column, multiple rows returned. Just take the rownames
      else if (ncol(x) == 1) {
        names(ret) <- x@rownames[i]
      }
      # Case: nrow, ncol > 1 we can only get a vector if length(i) == 1, in 
      # which case we take the column names
      else {
        names(ret) <- x@colnames
      }
    } else if (!is.null(x@rownames) & !is.null(x@colnames)) {
      dimnames(ret) <- list(x@rownames[i], x@colnames)
    } else if (!is.null(x@rownames)) {
      rownames(ret) <- x@rownames[i]
    } else if (!is.null(x@colnames)) {
      colnames(ret) <- x@colnames
    }
    ret
})

setMethod(
  "[", signature(x = "bigMatrix", i="missing", j="ANY", drop="ANY"), 
  function(x, i, j, drop){
    
    # Check if subscript type is valid
    if (!(class(j) %in% c("character", "integer", "numeric", "logical")))
      stop("invalid subscript type '", class(j), "'")   
    
    # Handle indexing by row or column names
    if (class(j) == "character")
      j <- match(j, x@colnames)
    if(any(is.na(j)))
      stop("subscript out of bounds")
    
    # The default behaviour is for drop to be TRUE
    if (is.null(drop) || is.na(drop) || class(drop) != "logical")
      drop <- TRUE
    
    # Attach the big.matrix object if not attached yet
    is.attached <- x@attached
    if (!is.attached)
      x <- attach.bigMatrix(x)
    
    # Fetch the results
    ret <- x@matrix[, j, drop=drop]
    
    # detach big.matrix object if it was detached to begin with
    if (!is.attached)
      x <- detach.bigMatrix(x)
    
    # handle the row and column names
    if (is.null(dim(ret)) & length(ret) != 0) {
      # Case: 1 row, length(ret) == 1, we take EITHER the row or column
      # names if only one exists. Otherwise, no names
      if (nrow(x) == 1 & length(ret) == 1) {
        if (!is.null(x@colnames) & is.null(x@rownames)) {
          names(ret) <- x@colnames[j]
        } else if (is.null(x@colnames) & !is.null(x@rownames)) {
          names(ret) <- x@rownames
        }
      }
      # Case: 1 column, multiple rows returned. Just take the rownames
      else if (nrow(x) == 1) {
        names(ret) <- x@colnames[j]
      }
      # Case: nrow, ncol > 1 we can only get a vector if length(j) == 1, in 
      # which case we take the column names
      else {
        names(ret) <- x@rownames
      }
    } else if (!is.null(x@rownames) & !is.null(x@colnames)) {
      dimnames(ret) <- list(x@rownames, x@colnames[j])
    } else if (!is.null(x@rownames)) {
      rownames(ret) <- x@rownames
    } else if (!is.null(x@colnames)) {
      colnames(ret) <- x@colnames[j]
    }
    ret
})

setMethod(
  "[", signature(x = "bigMatrix", i="missing", j="missing", drop="ANY"), 
  function(x, i, j, drop){
    
    # The default behaviour is for drop to be TRUE
    if (is.null(drop) || is.na(drop) || class(drop) != "logical")
      drop <- TRUE
    
    # Attach the big.matrix object if not attached yet
    is.attached <- x@attached
    if (!is.attached)
      x <- attach.bigMatrix(x)
    
    # Fetch the results
    ret <- x@matrix[,,drop=drop]
    
    # detach big.matrix object if it was detached to begin with
    if (!is.attached)
      x <- detach.bigMatrix(x)
    
    # Give the results the approriate row and column names
    if (is.null(dim(ret)) & length(ret) != 0) {
      if (nrow(x) == 1) {
        names(ret) <- x@colnames
      } else {
        names(ret) <- x@rownames
      }
    } else if (!is.null(x@rownames) & !is.null(x@colnames)) {
      dimnames(ret) <- list(x@rownames, x@colnames)
    } else if (!is.null(x@rownames)) {
      rownames(ret) <- x@rownames
    } else if (!is.null(x@colnames)) {
      colnames(ret) <- x@colnames
    }
    ret
})

#-------------------------------------------------------------------------------
# Matrix assignment methods
#-------------------------------------------------------------------------------
setMethod(
  "[<-", signature(x = "bigMatrix", i="ANY", j="ANY"), 
  function(x, i, j, value){
    
    # Check if subscript type is valid
    if (!(class(i) %in% c("character", "integer", "numeric", "logical")))
      stop("invalid subscript type '", class(i), "'")
    if (!(class(j) %in% c("character", "integer", "numeric", "logical")))
      stop("invalid subscript type '", class(j), "'")
    
    # Handle indexing by row or column names
    if (class(i) == "character")
      i <- match(i, x@rownames)
    if (class(j) == "character")
      j <- match(j, x@colnames)
    if(any(is.na(i)) || any(is.na(j)))
      stop("subscript out of bounds")
    
    # Attach the big.matrix object if not attached yet
    is.attached <- x@attached
    if (!is.attached)
      x <- attach.bigMatrix(x)
    
    # assign the value to the appropriate location
    x@matrix[i, j] <- value
    
    # detach big.matrix object if it was detached to begin with
    if (!is.attached)
      x <- detach.bigMatrix(x)
    
    x
})

setMethod(
  "[<-", signature(x = "bigMatrix", i="ANY", j="missing"), 
  function(x, i, j, value){
    
    # Check if subscript type is valid
    if (!(class(i) %in% c("character", "integer", "numeric", "logical")))
      stop("invalid subscript type '", class(i), "'")
    
    # Handle indexing by row or column names
    if (class(i) == "character")
      i <- match(i, x@rownames)
    if(any(is.na(i)))
      stop("subscript out of bounds")
    
    # Attach the big.matrix object if not attached yet
    is.attached <- x@attached
    if (!is.attached)
      x <- attach.bigMatrix(x)
    
    # assign the value to the appropriate location
    x@matrix[i,]<-value
    
    # detach big.matrix object if it was detached to begin with
    if (!is.attached)
      x <- detach.bigMatrix(x)
    
    x
})

setMethod(
  "[<-", signature(x = "bigMatrix", i="missing", j="ANY"), 
  function(x, i, j, value){
    
    # Check if subscript type is valid
    if (!(class(j) %in% c("character", "integer", "numeric", "logical")))
      stop("invalid subscript type '", class(j), "'")   
    
    # Handle indexing by row or column names
    if (class(j) == "character")
      j <- match(j, x@colnames)
    if(any(is.na(j)))
      stop("subscript out of bounds")
    
    # Attach the big.matrix object if not attached yet
    is.attached <- x@attached
    if (!is.attached)
      x <- attach.bigMatrix(x)
    
    # assign the value to the appropriate location
    x@matrix[,j]<-value
    
    # detach big.matrix object if it was detached to begin with
    if (!is.attached)
      x <- detach.bigMatrix(x)
    
    x
})

setMethod(
  "[<-", signature(x = "bigMatrix", i="missing", j="missing"), 
  function(x, i, j, value){
    
    # Attach the big.matrix object if not attached yet
    is.attached <- x@attached
    if (!is.attached)
      x <- attach.bigMatrix(x)
    
    # assign the value to the appropriate location
    x@matrix[,j]<-value
    
    # detach big.matrix object if it was detached to begin with
    if (!is.attached)
      x <- detach.bigMatrix(x)
    
    x
})
