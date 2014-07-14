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
#' @param ind an integer corresponding to the level of indentation. Each
#'   indentation level corresponds to two spaces.
#' @param nChunks The total number of chunks.
#' @seealso \code{\link[utils]{txtProgressBar}}
#' @name parProgress
NULL

#' @param chunk The chunk of indices for the current parallel instance of the 
#'   \code{foreach} loop.
#' @return
#'   an object of class "\code{txtProgressBar}"
#' @rdname parProgress
#' @importFrom utils txtProgressBar
setupParProgressLogs <- function(chunk, nChunks, ind) {
  chunkNum <- ceiling(chunk[1]/length(chunk))
  min <- chunk[1] - 1
  max <- chunk[length(chunk)]
  logFile <- file.path("run-progress", paste0("chunk", chunkNum, ".log"))
  file.create(logFile)
  # In our monitoring code, we will rotate through each chunk, printing out 
  # something along the lines of: 
  #  Chunk N:  |========                 | 30%
  # The 2 corresponds to ": "
  progWidth = options(width) - ind*2 - nchar("Chunks ") - nchar(nChunks) - 2
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
      progress <- readLines(file, warn=FALSE)
      nSpaces <- nchar(nChunks) - nchar(num)
      cat(sep="", "\r", rep("  ", ind), "Chunk ", rep(" ", nSpaces), num, ": ", 
          prog, file=stdout())
      Sys.sleep(2)
    }
  }
  cat("\n")
  invisible() # Nothing to return
}
