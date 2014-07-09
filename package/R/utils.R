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
#' @export
`%nin%` <- function(x, table) !(x %in% table)

#' @rdname matchsub
#' @export
`%sub_in%` <- function(x, table) x[x %in% table]

#' @rdname matchsub
#' @export
`%sub_nin%` <- function(x, table) x[x %nin% table]

#' Verbose Concatenate and Print with Indentation
#' 
#' Concatenate and output the objects only if the \code{verbose} flag is set 
#' to \code{TRUE}. Allows for indentation, adding a series of spaces to the 
#' beginning of each line, 2 for every increment in \code{indLevel}. 
#' 
#' @details
#'  \code{vCat} is slightly more intelligent than regulat \code{cat} in the way
#'  it formats the output and breaks it into lines. As a result, \code{fill} is
#'  set to \code{TRUE} by default. The other notable difference from \code{cat}
#'  is the way newline objects are handled. For example, the call:
#'  \code{cat("hello", "world", "\n", "foo", "bar")} won't wrap the newline 
#'  character with spaces. This avoids the need to set \code{sep} to \code{""} 
#'  and embed multiple \code{paste} calls.
#' 
#' @seealso \code{\link[base]{cat}}
#' @param verbose logical. If \code{TRUE}, passes the rest of the arguments to
#'   \code{\link{cat}}
#' @param indLevel an integer corresponding to the level of indentation. Each
#'   indentation level corresponds to two spaces.
#' @param ... Arguments to pass to \code{\link[base]{cat}}
#' @param sep a character vector of strings to append after each element.
#' @param fill a logical or (positive) numeric controlling how the output is 
#'   broken into successive lines. If \code{FALSE}, only newlines created 
#'   explicitly by "\n" are printed. Otherwise, the output is broken into lines
#'   with print width equal to the option width if fill is \code{TRUE} 
#'   (default), or the value of fill if this is numeric. Non-positive fill 
#'   values are ignored, with a warning.
#' @param labels character vector of labels for the lines printed. Ignored if 
#'   fill is \code{FALSE}.
vCat <- function(verbose, ..., indLevel=0, sep=" ", fill=TRUE, labels=NULL) {
  if(verbose) {
    # We need to format each line with the indendation level
    if (indLevel > 0) {
      indent <- paste(rep("  ", indLevel), collapse="")
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
            newlines[curnl] <- paste0(labels, " ", indent, w)
          } else if(nchar(newlines[curnl]) + 1 + nchar(w) < fillWidth) {
            newlines[curnl] <- paste(newlines[curnl], w)
          } else {
            curnl <- curnl + 1
            newlines[curnl] <- paste0(labels, " ", indent, w)
          }
        }
        paste(newlines, collapse="\n")
      })
      lines <- strsplit(paste(formatted, collapse="\n"), "\n")[[1]]
    }
    str <- paste0(lines, collapse="\n")
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

#' Binds two dimensional arrays along the third dimension.
#' @importFrom abind abind
abind3 <- function(...) abind(..., along=3)
