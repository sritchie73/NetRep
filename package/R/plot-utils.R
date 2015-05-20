#' Map a variable to a color gradient
#'
#' For a color gradient mapped onto a scale of values (\code{vlim}), retrieve
#' the color for any given \code{value}
#' 
#' @param values a vector of values to retrieve the colors for.
#' @param palette a vector of colors
#' @param vlim limits of the values
#' 
#' @return
#'  The color for a given value. 
getColFromPalette <- function(values, palette, vlim) {
  if (missing(vlim)) {
    vlim <- range(values)                                                                                                                                                                                                                     
  }
  vpal <- seq(vlim[1], vlim[2], length=length(palette)+1)
  
  sapply(values, function(vv) {
    if (vv >= vlim[2]) {
      tail(palette, n=1)
    } else if (vv <= vlim[1]) {
      palette[1]
    } else {
      palette[sum(vv > vpal)]
    }   
  })  
}

#' Get the plot limits to set for the desired plot window
#'
#' \code{plot} manually adds an extra region to the plot on top of the given
#' 'xlim' and 'ylim', amounting to 4% of the plot region in either direction.
#' This function tells you what limits to set so that the boundaries of the plot
#' are precisely the min and max you want.
#' 
#' @param dlim the limits you want to set
#' @return
#'  the min and max to set for the desired plot limits  
forceLim <- function(dlim) {
  A = matrix(c(1.04, -0.04, -0.04, 1.04), nrow=2, ncol=2)
  B = matrix(dlim, nrow=2, ncol=1)
  as.vector(solve(A, B))
}

#' Create an empty plot window with the specified limits
#' 
#' @param ... extra arguments to pass to \code{\link[graphics]{plot}}.
#' @param xlim limits for the x axis.
#' @param ylim limits for the y axis.
#' @param hardlim logical; if \code{TRUE}, the plot window is exactly 
#'  constrained to the specified \code{xlim} and \code{ylim}.
#'  
emptyPlot <- function(..., xlim, ylim, hardlim=TRUE) {
  if (!missing(xlim) & !missing(ylim)) {
    if(hardlim) {
      xlim <- forceLim(xlim)
      ylim <- forceLim(ylim)
    }
    plot(
      ..., 0, type='n', xaxt='n', yaxt='n', xlim=xlim, ylim=ylim, xlab="", ylab=""
    )
  } else if (!missing(ylim)) {
    if (hardlim) {
      ylim <- forceLim(ylim)
    }
    plot(
      ..., ylim=ylim, type='n', xaxt='n', yaxt='n', xlab="", ylab=""
    )
  } else if (!missing(xlim)) {
    if (hardlim) {
      xlim <- forceLim(xlim)
    }
    plot(
      ..., xlim=xlim, type='n', xaxt='n', yaxt='n', xlab="", ylab=""
    )
  } else {
    plot(
      ..., type='n', xaxt='n', yaxt='n', xlab="", ylab=""
    )
  }
}

#' Get the hexidecimal value for a named color
#' 
#' @param name the color name
#' @return the hex value for that color
col2hex <- function(name) { 
  paste0("#", paste(as.hexmode(col2rgb(name)), collapse="")) 
}

#' Adds alpha transparency to any color 
#' 
#' @param col color
#' @param alpha alpha transparency, should be between 0 and 1
#' @return the hexidecimal value of the color with the added alpha channel
addAlpha <- function(col, alpha) {
  if (!grepl("^#", col)) {
    col <- col2hex(col)
  }
  paste0(col, as.hexmode(floor(255*alpha)))
}
