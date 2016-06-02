#' Value Matching and Subsetting
#' 
#' @description
#' This set of functions provides shortcuts for value matching and subsetting,
#' on top of the functionality provided by \code{\link[base]{\%in\%}}.
#' 
#' \code{\%nin\%} returns a logical vector indicating if elements of \code{x} 
#' are not in \code{table}, This is the opposite of \code{\%in\%}.
#' 
#' \code{\%sub_in\%} returns the elements \code{x} that are \code{\%in\%} 
#' \code{table} rather than a logical vector.
#' 
#' \code{\%sub_nin\%} returns the elements \code{x} that are \code{\%nin\%} 
#' \code{table} rather than a logical vector.
#' 
#' @param x vector or \code{NULL}: the values to be matched. \link[base]{Long vectors} 
#'   are supported.
#' @param table vector or \code{NULL}: the values to be matched against. 
#'   \link[base]{Long vectors} are not supported.
#' @name matchsub
NULL

#' @rdname matchsub
`%nin%` <- function(x, table) !(x %in% table)

#' @rdname matchsub
`%sub_in%` <- function(x, table) x[x %in% table]

#' @rdname matchsub
`%sub_nin%` <- function(x, table) x[x %nin% table]

#' Verbose Concatenate and Print with Indentation
#' 
#' Concatenate and output the objects only if the \code{verbose} flag is set 
#' to \code{TRUE}. Allows for indentation, adding a series of spaces to the 
#' beginning of each line, 2 for every increment in \code{ind}. 
#' 
#' @details
#'  \code{vCat} is slightly more intelligent than regulat \code{cat} in the way 
#'  it formats the output and breaks it into lines. As a result, \code{fill} is 
#'  set to \code{TRUE} by default. Another notable difference from \code{cat} is
#'  the way newline objects are handled. For example, the call: 
#'  \code{cat("hello", "world", "\\n", "foo", "bar")} won't wrap the newline 
#'  character with spaces. This avoids the need to set \code{sep} to \code{""} 
#'  and embed multiple \code{paste} calls. Finally, a newline character is 
#'  appended to the end of the whole message, avoiding the need to manually 
#'  specify this when calling \code{vCat}.
#'  
#' @seealso \code{\link[base]{cat}}
#' @param verbose logical. If \code{TRUE}, passes the rest of the arguments to
#'   \code{\link{cat}}
#' @param ind an integer corresponding to the level of indentation. Each
#'   indentation level corresponds to two spaces.
#' @param ... Arguments to pass to \code{\link[base]{cat}}
#' @param sep a character vector of strings to append after each element.
#' @param fill a logical or (positive) numeric controlling how the output is 
#'   broken into successive lines. If \code{FALSE}, only newlines created 
#'   explicitly by "\\n" are printed. Otherwise, the output is broken into lines
#'   with print width equal to the option width if fill is \code{TRUE} 
#'   (default), or the value of fill if this is numeric. Non-positive fill 
#'   values are ignored, with a warning.
#' @param labels character vector of labels for the lines printed. Ignored if 
#'   fill is \code{FALSE}.
vCat <- function(verbose, ind=0,  ..., sep=" ", fill=TRUE, labels=NULL) {
  if (!(is.vector(verbose) && !is.list(verbose) && is.logical(verbose) &&
        length(verbose) == 1 && !is.na(verbose))) {
    stop("'verbose' must be one of 'TRUE' or 'FALSE'")
  }
  
  if(verbose) {
    # Put a timestamp at the start of the message:
    if (is.null(labels))
      labels <- paste0("[", format(Sys.time(), usetz=TRUE), "] ")
    
    # We need to format each line with the indendation level
    if (ind > 0) {
      indent <- paste(rep("  ", ind), collapse="")
    } else {
      indent = ""
    }
    
    args <- list(...)
    if (is.null(names(args))) {
      str <- paste(args, collapse=sep)
      named <- NULL
    } else {
      str <- paste(args[names(args) == ""], collapse=sep)
      named <- args[names(args) != ""]
    }
    str <- gsub(" \n ", "\n", str) # make it easier to insert newlines
    lines <- strsplit(str, "\n")[[1]]
    # Handle automatic line wrapping
    if (fill) {
      if (is.logical(fill)) {
        fillWidth <- options("width")
      } else if (!is.numeric(fill)) {
        stop("invalid 'fill' argument")
      } else if (fill < 1) {
        warning("non-positive 'fill' argument will be ignored")
        fillWidth <- options("width")
      } else {
        fillWidth <- fill
      }
      words <- strsplit(lines, " ")
      # Create new lines by accumulating words in each line until the addition
      # of the next word would exceed the fillWidth. 
      formatted <- lapply(words, function(lw) {
        newlines <- c("")
        curnl <- 1
        for (w in lw) {
          if (newlines[curnl] == "") {
            newlines[curnl] <- paste0(labels, indent, w)
          } else if(nchar(newlines[curnl]) + 1 + nchar(w) < fillWidth) {
            newlines[curnl] <- paste(newlines[curnl], w)
          } else {
            curnl <- curnl + 1
            labels <- paste(rep(" ", nchar(labels)), collapse="")
            newlines[curnl] <- paste0(labels, indent, w)
          }
        }
        paste(newlines, collapse="\n")
      })
      lines <- strsplit(paste(formatted, collapse="\n"), "\n")[[1]]
    }
    str <- paste0(lines, "\n", collapse="")
    if (is.null(named)) {
      cat(str)
    } else {
      # build expression from remaining named arguments
      args = sapply(seq_along(named), function(n) {
        if (is.character(named[[n]])) {
          paste0(names(named)[n], "=", "'", named[[n]], "'")
        } else {
          paste0(names(named)[n], "=", named[[n]])
        }
      })
      eval(parse(text=paste0(paste0(c("cat(str", args), collapse=", "), ")")))
    }
  }
}

#' Poke attached big matrix objects to initialise them
#' 
#' This is a magical speed hack. See details
#' 
#' @details
#' For some reason, one of the datasets i tried running on took an incredibly
#' long time for some very small network subsets. After playing around it turns
#' out that the first time any of the Rcpp functions are called, regardless of 
#' subset size, it was taking forever. This can be solved simply by accessing
#' the big.matrix objects beforehand.
#' 
#' @param ... a number of big.matrix objects
poke <- function(...) {
  for (o in list(...)) {
    if (!is.null(o)) {
      pCols <- min(ncol(o), 5)
      pRows <- min(nrow(o), 5)
      o[1:pRows, 1:pCols]
    }
  }
}

#' Get a universally unique identifier
#' 
#' Thanks to thelatemail's answer on stackoverflow:
#' http://stackoverflow.com/questions/10492817/how-can-i-generate-a-guid-in-r
#' 
getUUID <- function() {
  baseuuid <- paste(sample(c(letters[1:6],0:9),30,replace=TRUE),collapse="")
  
  paste(
    substr(baseuuid,1,8),
    "-",
    substr(baseuuid,9,12),
    "-",
    "4",
    substr(baseuuid,13,15),
    "-",
    sample(c("8","9","a","b"),1),
    substr(baseuuid,16,18),
    "-",
    substr(baseuuid,19,30),
    sep="",
    collapse=""
  )
}

#' Insert NAs into a vector at specified positions
#' 
#' Useful for inserting NAs into the correct positions when examining module 
#' probes that do not exist in the test dataset.
#' 
#' @param vec vector to insert NAs to.
#' @param na.indices indices the NAs should be located at in the final vector.
#' 
#' @return
#' The vector with NAs inserted in the correct positions.
insert.nas <- function(vec, na.indices) {
  res <- vector(typeof(vec), length(vec) + length(na.indices))
  res[na.indices] <- NA
  res[!is.na(res)] <- vec
  res
}

#' Order the module vector numerically
#' 
#' The module assingments may be numeric, but coded as characters.
#' 
#' @param vec module vector to order
#' 
#' @return the order of the vector
orderAsNumeric <- function(vec) {
  tryCatch({
    order(as.integer(vec))
  }, warning=function(w) {
    order(vec)
  })
}

#' Remove unnecessary list structure 
#'
#' Removes entries that are \code{NULL} are extracts element if the length is 1.
#'
#' @param l a nested list
#' @param depth depth to traverse to
#' 
#' @return a list
simplifyList <- function(l, depth) {
  # Recursively traverse until we hit the depth requested or we cant go deeper
  stopifnot(is.numeric(depth) && depth > 0)
  if (depth == 1) {
    if (is.null(l)) {
      return(NULL)
    }
    # Delete empty leaf nodes
    for (ii in rev(seq_along(l))) {
      if (is.null(l[[ii]]))
        l[[ii]] <- NULL
    }
    if (length(l) == 0) {
      return(NULL)
    } else if (length(l) == 1) {
      return(l[[1]])
    } else {
      return(l) 
    }
  } else {
    for (ii in rev(seq_along(l))) {
      l[[ii]] <- simplifyList(l[[ii]], depth-1)
    }
    if (length(l) == 0) {
      return(NULL)
    } else if (length(l) == 1) {
      return(l[[1]])
    } else {
      return(l) 
    }
  }
}

#' Get a sorted list of module names
#' 
#' If module labels are numeric, sorts numerically, otherwise sorts 
#' alphabetically.
#' 
#' @param modules a vector of module labels to sort
#' @return a sorted vector   
sortModuleNames <- function(modules) {
  tryCatch({
    modules[order(as.numeric(modules))]
  }, warning=function(w) {
    sort(modules)
  })
}

#' Silently check and load a package into the namespace
#' 
#' @param pkg name of the package to check
#' 
#' @return logical; \code{TRUE} if the package is installed and can be loaded.
pkgReqCheck <- function(pkg) {
  suppressMessages(suppressWarnings(requireNamespace(pkg, quietly=TRUE)))
}
