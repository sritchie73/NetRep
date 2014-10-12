## Authors:
##   Scott Ritchie <sritchie73@gmail.com>
##
## File Purpose:
##   To provide a tool for plotting a gradient on a bar with arbitrary break
##   points. Useful for plotting bars with variable width bins, as well as
##   plotting legends.
##

#' @title Plot a Gradient along a Bar
#'
#' @description
#' \code{gradient.bar} plots a gradient along a bar in either the x or y
#' directions, using user specified break points for each step in the gradient.
#' This works by plotting the gradient over the entire dimensions of a plot
#' window.
#'
#' @param range The range of the plot in the direction of the gradient (e.g.
#'   xlim if the gradient is plotted over the x axis). Useful if the
#'   \code{break.points} specified by the user have been determined by some
#'   other function over some range.
#' @param break.points A vector of locations where each step of the gradient
#'   transitions into the next.
#' @param nBins Splits the gradient bar into \code{nBins} evenly spaced bins.
#'   Ignored if \code{break.points} are given.
#' @param col A vector of colors used to color each bin in the gradient bar.
#' @param col.gradient Colors to generate a gradient between. For example,
#'   c("white", "red"), or, c("red", "white", "blue"). See
#'   \code{help(colorRampPalette)} for more details. Ignored if \code{col} is
#'   provided.
#' @param border The color of the border around the entire bar, or NA for no
#'   border.
#' @param lines The color(s) of the lines to split each gradient step, or NA
#'   to not draw lines.
#' @param direction The direction to render the gradient, must be one of 'x' or
#'   'y'. To reverse the direction the gradient is draw on an axis, simply
#'   reverse the color vectors provided to either \code{code} or
#'   \code{col.gradient}.
#' @param main Text to render above the gradient bar.
#' @param bin.lab Text labels for each of the gradient bins. Will either plot 
#'    the labels next to the bins, or the breaks between bins, depending on how
#'    many labels are provided. Labels are rendered either on the left, or
#'    bottom of the graph depending on the direction of the gradient. Useful for
#'    labelling legends
#'    
#' @keywords gradient legend
#' @examples
#' # Intervals of size 2^x
#' gradient.bar(range=c(0, 1024), break.points=2^(1:10),
#'   col.gradient=c("red", "white"))
#'
#' # Plot a legend next to a bar plot
#' layout(matrix(c(1, 2), nrow=1), widths=c(0.8, 0.2))
#' barplot(1:10, col=1:10)
#' gradient.bar(nBins=10, col=1:10, direction='y', main="legend",
#'              bin.lab=letters[1:10], lines="black")
#'
#' # Rainbow box
#' gradient.bar(nBins=256, col=rainbow(256))
#'
#' @export
gradient.bar <- function(range=c(0,1), break.points=NULL, nBins=NULL, col=NULL,
                         col.gradient=NULL, border=par('fg'), lines=NA,
                         direction="x", main="", bin.lab=NULL) {
  # Error Check
  if (!(direction %in% c("x", "y"))) {
    stop("Error: Direction must be either 'x' or 'y'.")
  }
  
  # Helper function - makes sure break points are within range and match up
  # with boundaries
  fix.points <- function(points, range) {
    if (is.null(points)) {
      return(range)
    }
    points <- sort(points)
    points <- points[which(points > range[1])]
    points <- points[which(points < range[2])]
    return(c(range[1], points, range[2]))
  }
  
  # Handle mutually exclusive arguments `break.points` and `nBins`.
  # `break.points` overrides `nBins`.
  if (!is.null(break.points)) {
    if(!is.null(nBins)) {
      warning("Ignoring 'nBins' arguments, using 'break.points' instead.")
      nBins <- NULL
    }
    break.points <- fix.points(break.points, range)
  } else if (!is.null(nBins)) {  # Use the number of nBins to determine break.points
    break.points <- seq(range[1], range[2], l=(nBins + 1))
    break.points <- fix.points(break.points, range)
  } else {
    break.points <- c()
  }  
  
  # Handle mutually exlusive arguments `col` and `col.gradient`.
  # `col` overrides `col.gradient`
  if (!is.null(col)) {
    if ((length(break.points) - 1) %% length(col)) {
      stop("Error: number of colors in 'col' is not a replacement length of the",
           " number of bins in the gradient.")
    }
    if (!is.null(col.gradient)) {
      warning("Ignoring 'col.gradient' argument, using supplied colors in 'col'",
              "instead.")
      col.gradient <- NULL
    }
  } else if (is.null(col.gradient)) {
    col <- par('bg')
  }
  
  # Create blank template plot
  # Set the limits of the graph according to the range and direction
  
  if (direction == 'x') {
    plot(NA, type='n', ann=FALSE, xlim=range, ylim=c(0,1), xaxt='n', yaxt='n',
         bty='n')
  } else {
    plot(NA, type='n', ann=FALSE, xlim=c(0,1), ylim=range, xaxt='n', yaxt='n',
         bty='n')
  }
  
  # Create gradient
  if (!is.null(col.gradient)) {
    col <- colorRampPalette(col.gradient)(length(break.points)-1)
  }
  
  # Plot gradient
  if (direction == "x") {
    rect(head(break.points, -1), 0, tail(break.points, -1), 1,
         col=col, border=lines)
  } else  {
    break.points <- rev(break.points)  # So it is rendered top down
    rect(0, head(break.points, -1), 1, tail(break.points, -1),
         col=col, border=lines)
  }
  
  # Render border over the top if specified
  if (!is.null(border)) {
    par(new=TRUE)
    if (direction == 'x') {
      rect(range[1], 0, range[2], 1, col="#00000000", border=border)
    } else {
      rect(0, range[1], 1, range[2], col="#00000000", border=border)
    }
  }
  
  # Handle text to label each bin with. Useful for legend plotting
  if (!is.null(bin.lab)) {
    side <- switch(direction, "x" = 1, "y" = 2)
    if (length(bin.lab) == length(break.points)) {
      at <- break.points
    } else if (length(bin.lab) == (length(break.points) - 1)) {
      # Get the midpoint for each interval
      offsets <- sapply(head(seq_along(break.points),-1), function(i) {
        (break.points[i+1] - break.points[i])/2
      })
      at <- head(break.points, -1) + offsets
    } else {
      stop("Invalid number of labels provided for bin.lab!")
    }
    mtext(bin.lab, side=side, at=at, 0.5, las=2)
  }
  if (direction  == "x") {
    at <- (range[2] - range[1])/2
  } else {
    at <- 0.5
  }
  mtext(main, side=3, at=at, 0, las=1, cex=1)
}
