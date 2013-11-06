#' Prints the message only if the verbose flag is set.
#' 
#' @param ... message to print
#' @param verbose Logical; prints message if true
VerbosePrint <- function(..., verbose) {
  if (verbose) {
    print(paste(...))
  }
}

#' Prints the progress of a looping function each time the loop progresses by
#' at least 1\%.
#'
#' @param this the current iteration number.
#' @param total the total number of iterations
PrintProgress <- function(this, total) {  
  # Error Checking
  CheckArgument(this, class=c("numeric", "integer"))
  CheckArgument(total, class=c("numeric", "integer"))
  
  last <- floor((this-1)/total*100)  # % complete of previous iteration
  cur <- floor(this/total*100)  # % completed of the current iteration
  if (cur > last) {
    print(paste(cur, "% completed", sep=""))
  }
  return()
}

#' Log the progress of a looping function to file.
#' 
#' Outputs the progress of a looping function to a file. This can then be 
#' monitored using the unix command \code{tail -f log.txt}. Useful for parallel
#' loops, where \code{\link{PrintProgress}} won't be able to output to the R terminal
#'
#' @param this the current iteration number.
#' @param total the total number of iterations
#' @param file The log file to output to.
LogProgress <- function(this, total, file="log.txt") {  
  # Error Checking
  CheckArgument(this, class=c("numeric", "integer"))
  CheckArgument(total, class=c("numeric", "integer"))
  
  if(!file.exists(file)) { 
    cat("Starting Loop!\n", file=file)
  }
  
  last <- floor((this-1)/total*100)  # % complete of previous iteration
  cur <- floor(this/total*100)  # % completed of the current iteration
  if (cur > last) {
    cat(paste(cur, "% completed\n", sep=""), file=file, append=TRUE)
  }
  return()
}
