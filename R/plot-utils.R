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
#' @param xlab label for the x axis.
#' @param ylab label for the y axis.
#' @param hardlim logical; if \code{TRUE}, the plot window is exactly
#'  constrained to the specified \code{xlim} and \code{ylim}.
#'
emptyPlot <- function(..., xlim, ylim, xlab="", ylab="", hardlim=TRUE) {
  if (!missing(xlim) && !missing(ylim)) {
    if(hardlim) {
      xlim <- forceLim(xlim)
      ylim <- forceLim(ylim)
    }
    plot(
      ..., 0, type='n', xaxt='n', yaxt='n', xlim=xlim, ylim=ylim, xlab=xlab,
      ylab=ylab
    )
  } else if (!missing(ylim)) {
    if (hardlim) {
      ylim <- forceLim(ylim)
    }
    plot(
      ..., ylim=ylim, type='n', xaxt='n', yaxt='n', xlab=xlab, ylab=ylab
    )
  } else if (!missing(xlim)) {
    if (hardlim) {
      xlim <- forceLim(xlim)
    }
    plot(
      ..., xlim=xlim, type='n', xaxt='n', yaxt='n', xlab=xlab, ylab=ylab
    )
  } else {
    plot(
      ..., type='n', xaxt='n', yaxt='n', xlab=xlab, ylab=ylab
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

#' Color palette for correlation heatmaps
#'
#' RColorBrewer palette "RdYlBu" with the middle color replaced with white.
#' This gives a nicer contrast than the "RdBu" palette
#' 
#' @import RColorBrewer
correlation.palette <- function() {
  cols <- rev(brewer.pal(11, "RdYlBu"))
  cols[6] <- "#FFFFFF"
  cols
}

#' Color palette for network heatmaps
#'
#' RColorBrewer palette "RdYlBu" with the middle color replaced with white.
network.palette <- function() {
  correlation.palette()[6:11]
}

#' Get module break points on the x-axis
#'
#' @param mas ordered subset of the moduleAssignments vector
#'
#' @return
#'  a vector of positions on the x-axis where one module begins and another ends
getModuleBreaks <- function(mas) {
  sizes <- rle(mas)
  breaks = numeric(length(sizes$lengths) + 1)
  breaks[1] <- 0.5
  for (mi in seq_along(sizes$lengths)) {
    breaks[mi + 1] <- breaks[mi] + sizes$lengths[mi]
  }
  breaks
}

#' Get module mid-points on the x-axis
#'
#' @param mas ordered subset of the moduleAssignments vector
#'
#' @return
#'  a vector of positions on the x-axis indicating the centre of a module
getModuleMidPoints <- function(mas) {
  breaks <- getModuleBreaks(mas)
  mids <- numeric(length(breaks) - 1)
  for (bi in seq_along(breaks)[-1]) {
    mids[bi - 1] <- (breaks[bi] - breaks[bi - 1])/2 + breaks[bi - 1]
  }
  mids
}

#' Check if a character vector of colors is valid
#' 
#' Courtesy of Josh O'Brien's stackoverflow answer at
#' \url{http://stackoverflow.com/a/13290832/2341679}
#' 
#' @param colvec a character vectors of colors (hex or name) to validate.
areColors <- function(colvec) {
  sapply(colvec, function(col) {
    tryCatch(is.matrix(col2rgb(col)), error = function(e) FALSE)
  })
}
