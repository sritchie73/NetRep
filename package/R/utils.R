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
  if(verbose) {
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

#' Combine two-dimensional arrays
#' 
#' Binds two dimensional arrays along the third dimension.
#' 
#' @param ... any number of two dimensional objects, to be bound together along
#'   the third dimension.
#' @importFrom abind abind
abind3 <- function(...) abind(..., along=3)

#' Log Progress to File for Parallel Processes
#' 
#' Update progress of a parallel \code{\link[foreach]{foreach}} process to
#' files.
#' 
#' @details
#'   When running code in parallel using \code{foreach}, new R
#'   instances will be spawned, one for each core as setup by the parallel
#'   backend. As a result no output will be sent back to the main R process
#'   until the \code{foreach} loop has finished running.
#' 
#' @seealso \code{\link[utils]{txtProgressBar}}
#' @name parProgress
NULL


#' @param chunk The chunk of indices for the current parallel instance of the 
#'   \code{foreach} loop.
#' @param nChunks The total number of chunks.
#' @return
#'   an object of class "\code{txtProgressBar}"
#' @rdname parProgress
#' @importFrom utils txtProgressBar
setupParProgressLogs <- function(chunk, nChunks) {
  chunkNum <- ceiling(chunk[1]/length(chunk))
  min <- chunk[1] - 1
  max <- chunk[length(chunk)]
  logFile <- file.path("run-progress", paste0("chunk", chunkNum, ".log"))
  file.create(logFile)
  # In our monitoring code, we will rotate through each chunk, printing out 
  # something along the lines of: 
  #  Chunk N:  |========                 | 30%
  progWidth = options(width) - nchar("Chunks ") - nchar(nChunks) - 1
  txtProgressBar(
      min, max, min, char="=", width=progWidth, style=3, 
      file=file(logFile, open="wt")
    )
}

#' @param pb an object of class "\code{txtProgressBar}"
#' @param i new value for the progress bar.
#' @importFrom utils setTxtProgressBar
#' @rdname parProgress
updateParProgress <- function(pb, i) {
  setTxtProgressBar(pb, i)
}

#' Monitor Parallel Progress
#' 
#' Monitor the progress of null distribution calculation. 
#' 
#' @details
#'   This must be run in a new R session.
#' 
#' @note
#'  If you see the message \emph{"Waiting for parallel code to start..."} and 
#'  \code{\link{netRep}} is already calculating the null distributions, your R
#'  session has been initialized in the wrong directory.
#'  
#' @param updateFreq number or seconds to wait between chunk switches.
#' 
#' @export
monitorProgress <- function(updateFreq=2) {
  if(!file.exists("run-progress")) {
    vCat(TRUE, 0, "Waiting for parallel code to start...")
  }
  files <- list.files("run-progress")
  if (nChunks == 0) {
    nChunks <- length(files)
  }
  while(TRUE) {
    if (!file.exists("run-progress")) {
      break
    }
    for (file in files) {
      num <- gsub("Chunk|.log", "", file)
      progress <- readLines(file, warn=FALSE)
      nSpaces <- nchar(nChunks) - nchar(num)
      cat(sep="", "\rChunk ", rep(" ", nSpaces), num, ": ", prog)
      Sys.sleep(updateFreq)
    }
  }
  vCat(TRUE, 0, "\nAll Done!")
  invisible() # Nothing to return
}