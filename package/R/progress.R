#' Monitor Parallel Process Progress
#' 
#' Update progress of a parallel \code{\link[foreach]{foreach}} process to
#' files.
#' 
#' @details
#'   When running code in parallel using \code{foreach}, new R
#'   instances will be spawned, one for each core as setup by the parallel
#'   backend. As a result no output will be sent back to the main R process
#'   until the \code{foreach} loop has finished running.
#' @param ind an integer corresponding to the level of indentation. Each
#'   indentation level corresponds to two spaces.
#' @param nChunks The total number of chunks.
#' @seealso \code{\link[utils]{txtProgressBar}}
#' @name parProgress
NULL

#' @param chunk The chunk of indices for the current parallel instance of the 
#'   \code{foreach} loop.
#' @return
#'   \code{setupParProgressLogs}: an object of class "\code{txtProgressBar}"
#' @rdname parProgress
#' @importFrom utils txtProgressBar
setupParProgressLogs <- function(chunk, nChunks, ind) {
  chunkNum <- ceiling(chunk[1]/length(chunk))
  min <- chunk[1] - 1
  max <- chunk[length(chunk)]
  filename <- file.path("run-progress", paste0("chunk", chunkNum, ".log"))
  file.create(filename)
  logfile <- file(filename, open="wt")
  # In our monitoring code, we will rotate through each chunk, printing out 
  # something along the lines of: 
  #  Chunk N: |========                 | 30%
  # The 2 corresponds to ": "
  width <- options("width")[[1]]
  if (nChunks > 1) {
    progWidth = width - ind*2 - nchar("Chunks ") - nchar(nChunks) - 2
  } else {
    # If only one worker/core, no need to prepend with "Chunk N: "
    progWidth = width - ind*2 
  }
  txtProgressBar(min, max, min, width=progWidth, style=3)
}

#' @param pb an object of class "\code{txtProgressBar}"
#' @param i new value for the progress bar.
#' @importFrom utils setTxtProgressBar
#' @rdname parProgress
updateParProgress <- function(pb, i) {
  setTxtProgressBar(pb, i)
}

#' @description 
#'  \code{monitorProgress}: Monitor the progress of parallel workers.
#' @rdname parProgress
monitorProgress <- function(nChunks, ind) {
  init <- FALSE
  while(TRUE) {
    files <- list.files("run-progress")
    if (init & length(files) == 0) {
      break;
    } else {
      init <- TRUE
    }
    for (file in files) {
      num <- gsub("Chunk|.log", "", file)
      progress <- readLines(file.path("run-progress", file), warn=FALSE)
      if (nChunks > 1) {
        nSpaces <- nchar(nChunks) - nchar(num)
        cat(sep="", "\r", rep("  ", ind), "Chunk ", rep(" ", nSpaces), num, 
            ": ", progress, file=stdout())       
      } else {
        cat(sep="", "\r", rep("  ", ind), progress, file=stdout())
      }
      Sys.sleep(2)
    }
  }
  cat("\n")
  invisible() # Nothing to return
}

#' @description 
#'  \code{reportProgress}: Report the progress of a sequential loop.
#' @rdname parProgress
reportProgress <- function(ind) {
  progress <- readLines(
      file.path("run-progress", file="Chunk1.log"), warn=FALSE
    )
  cat(sep="", "\r", rep("  ", ind), progress, file=stdout())
}