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

#' Iterator with master worker
#' 
#' Create an iterator of chunks of indices from 1 to \code{n}, along with an 
#' initial element that designates the master worker. You can specify either the
#' number of pieces, using the \code{chunks} argument, or the maximum size of
#' the pieces, using the \code{chunkSize} argument.
#' 
#' @details
#'  If \code{verbose} is \code{FALSE}, the returned iterator is exactly the 
#'  same as if calling \code{\link[itertools]{isplitIndices}}. If \code{TRUE},
#'  then the iterator is prepended with a -1. This allows for a 
#'  \code{\link[foreach]{foreach}} loop to easily operate on chunks of tasks, 
#'  while also designating one worker thread as the master of reporting 
#'  progress (see \code{\link[=parProgress]{monitorProgress}}).
#' 
#' @seealso 
#'  \code{\link[itertools]{isplitIndices}} \code{\link[iterators]{idiv}}
#'  \code{\link[=parProgress]{monitorProgress}}
#' @param verbose logical. Controls the type of iterator returned, see details.
#' @param n Maximum index to generate.
#' @param cores the number of cores to divide \code{n} across. If \code{verbose}
#'  is \code{TRUE}, \code{n} is distributed over \code{cores - 1}, while 1 core
#'  is reserved as the task monitor.
#' @return
#'  An iterator that returns -1 (for the master worker), and vectors of indices 
#'  from 1 to \code{n} for the other worker threads.
#' @importFrom iterators idiv
#' @importFrom iterators nextElem
#' @importFrom itertools isplitIndices
ichunkTasks <- function(verbose, n, cores) {
  if (verbose & (cores > 1)) {
    it <- idiv(n, chunks=cores-1)
    i <- 1L
    first = TRUE
    nextEl <- function() {
      if (first) {
        first <<- FALSE
        -1L
      } else {
        m <- as.integer(nextElem(it))
        j <- i
        i <<- i + m
        seq(j, length=m)
      }
    }
    object <- list(nextElem = nextEl)
    class(object) <- c("abstractiter", "iter")
    object
  } else {
    isplitIndices(n, chunks=cores)
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

#' Floating Point Comparison
#'
#' Tests elements in a vector for equality to a specified floating point value.
#' 
#' @param vector a numeric vector of doubles
#' @param value a double
#' @return logical; \code{TRUE} where an element is equal to \code{value}, 
#'   \code{NA} where a comparison is made with an \code{NA}, and \code{FALSE}
#'   otherwise.
is.equal <- function(vector, value) {
  sapply(vector, function(element) {
    isTRUE(all.equal(element, value))
  })
}

#' Unify the datastructure to be a list of things
#' @param x object to convert
unifyDS <- function(x) {
  if (!is.list(x))
    x <- list(x)
  x
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

#' Combine null distribution calculations
#' 
#' @param ... any number of three dimensional arrays holding the null 
#'   distributions calculated by \code{\link{modulePreservation}}.
#'
#' @return
#'  A three dimensional array combined from the input. The rows correspond to
#'  the modules, columns to the module preservation statistics.
#'
#' @importFrom abind abind
#' @export
combineNulls <- function(...) {
  abind(..., along=3)
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

#' Get variables composing a module
#' 
#' Get the variables composing the user specified modules in the order they 
#' appear in the 'moduleAssignments' vector for the discovery dataset.
#' 
#' @param moduleAssignments a vector of the module assignment for each variable
#'  in the discovery dataset. If there are multiple discovery datasets 
#' then this argument should be a list of such vectors.  
#' @param modules a vector of modules to apply the function to.
#' @param discovery name or index denoting which dataset the module of
#'  interest was discovered in.
#'  
#' @rdname getModuleVarsUnsorted
getModuleVarsUnsorted <- function(
  moduleAssignments, modules, discovery=1
) {
  foreach(mi = modules, .combine=c) %do% {
    names(moduleAssignments[[discovery]] %sub_in% modules)
  }
}

#' Set up a parallel backend
#'
#' Set up a backend with the requested number of cores, or use existing backend
#' if the user has set one up already.
#'
#' @param nCores number of cores to use.
setupParallel <- function(nCores=NULL, verbose) {
  # First, check whether the user has already set up a parallel backend. In this
  # case, we can ignore the `nCores` argument.
  if (getDoParWorkers() > 1) {
    vCat(
      verbose, 0, "Ignoring 'nCores': parallel backend detected.", 
      "Reserving 1 core for progress reporting.",
      getDoParWorkers() - 1, "cores will be used for computation"
    )
    nCores <- getDoParWorkers()
  } 
  
  # If the user is on a Windows machine, we have to use the `doParallel` package 
  else if (.Platform$OS.type == "windows" & nCores > 1) {
    # Quietly load parallel backend packages. Throw our own warning and 
    # continue
    if(suppressWarnings(suppressMessages(requireNamespace("doParallel")))) {
      # we need an additional thread to monitor and report progress
      if (verbose)  
        nCores <- nCores + 1
      cl <- parallel::makeCluster(nCores)
      doParallel::registerDoParallel(cl)
      on.exit({
        parallel::stopCluster(cl)
      }, add=TRUE)
      vCat(verbose, 0, "Running on", nCores - 1, "cores.")
      if ((nCores - 1) > parallel::detectCores()) {
        stop(
          "Requested number of threads (", nCores - 1, ") is higher than the ",
          "number of available cores (", parallel::detectCores(), 
          "). Using too many threads may cause the machine to thrash/freeze."
        )
      }
    } else {
      nCores <- 1
      # We want to immediately print a warning for the user, not at the end 
      # once the analysis has finished.
      vCat(
        TRUE, 0, file=stderr(),
        "Warning: unable to find 'doParallel' package, running on 1 core." 
      )
    }
  } else if (.Platform$OS.type == "unix" & nCores > 1) {
    # Quietly load parallel backend packages. Throw our own warning and 
    # continue
    if(suppressWarnings(suppressMessages(requireNamespace("doMC")))) {
      # we need an additional thread to monitor and report progress
      if (verbose) 
        nCores <- nCores + 1
      doMC::registerDoMC(nCores)
      vCat(verbose, 0, "Running on", nCores - 1, "cores.")
      if ((nCores - 1) > parallel::detectCores()) {
        stop(
          "Requested number of threads (", nCores - 1, ") is higher than the ",
          "number of available cores (", parallel::detectCores(), 
          "). Using too many threads may cause the machine to thrash/freeze."
        )
      }
    } else {
      nCores <- 1
      # We want to immediately print a warning for the user, not at the end 
      # once the analysis has finished.
      vCat(
        TRUE, 0, file=stderr(),
        "Unable to find 'doMC' package, running on 1 core."
      )
    }
  } else {
    vCat(verbose, 0, "Running on 1 cores.")
  }
  
  # Suppress annoying foreach warning generated when using %dopar% and running 
  # in serial
  if (nCores == 1) {
    suppressWarnings({
      ii <- 0 # suppress R CMD check note
      foreach(ii = 1:2) %dopar% { ii }
    })
  }

  # Since we expect the user to explicitly handle the number of parallel threads,
  # we will disable the potential implicit parallelism on systems where R has
  # been compiled against a multithreaded BLAS, e.g. OpenBLAS. 
  omp_set_num_threads(1)
  blas_set_num_threads(1)
}

#' Remove unnecessary list structure at depth = 2
#'
#' Removes entries that are \code{NULL} are extracts element if the length is 1.
#'
#' @param l a nested list
#' 
#' @return a list
simplifyList2 <- function(l) {
  filterNulls <- function(l) {
    l <- l[!sapply(l, is.null)]
    if (length(l) == 0)
      return(list(NULL))
    return(l)
  }
  collapse <- function(l) {
    if (length(l) == 1)
      return(l[[1]])
    return(l)
  }
  # At depth 2
  for (i1 in rev(seq_along(l))) {
    l[[i1]] <- filterNulls(l[[i1]])
    l[[i1]] <- collapse(l[[i1]])
  }
  # At depth 1
  l <- filterNulls(l)
  l <- collapse(l)
  return(l)
}


#' Remove temporarily created bigMatrix objects
#' 
#' Prevents the temporary directory from growing excessively large if using many
#' package functions in the same R session
#'  
cleanTempDir <- function() {
  for (pat in c("*.bin", "*.desc", "*_rownames.txt", "*_colnames.txt")) {
    path <- file.path(tempdir(), list.files(tempdir(), pat))
    unlink(path)
  }
}
