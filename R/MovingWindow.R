#' MovingWindow returns a moving average of a vector.
#' @title Return a moving average of a vector.
#' @author Marc Girondot
#' @return A vector
#' @param x The vector to analyze
#' @param window The window size
#' @param fill TRUE or FALSE, should the vector return NA
#' @param hole Should the returned vector have the same length than x
#' @param FUN Function to apply to the window
#' @description Return a moving average of a vector./cr
#' hole parameter can be none, bothL, bothR, both, begin, end.
#' @examples
#' MovingWindow(1:10, window = 4, fill = TRUE, hole="bothL")
#' MovingWindow(1:10, window = 4, fill = TRUE, hole="bothR")
#' MovingWindow(1:10, window = 4, fill = TRUE, hole="both")
#' MovingWindow(1:10, window = 4, fill = TRUE, hole="none")
#' MovingWindow(1:10, window = 4, fill = TRUE, hole="begin")
#' MovingWindow(1:10, window = 4, fill = TRUE, hole="end")
#' MovingWindow(1:10, window = 4, fill = TRUE, hole="end", FUN=sd)
#' @export

MovingWindow <- function(x, window, hole="begin", fill=TRUE, FUN=mean) {
  fill <- tolower(fill)
  if (!is.na(fill))  hole <- match.arg(hole, choices = c("begin", "end", "bothL", "bothR", "both", "none"))
  y <- rep(NA, length(x))
  for (i in window:length(x)) {
    y[i] <- FUN(x[(i-window+1):i])
  }
  if (hole=="none") y <- y[is.na(y)]
  if (hole=="begin") {
    if (fill) {
      y[1:(window-1)] <- y[window]
    }
  }
  if (hole=="end") {
    if (fill) {
      y <- c(y[window:length(y)], rep(y[length(y)], window-1))
    } else {
      y <- c(y[window:length(y)], rep(NA, window-1))
    }
  }
  if (hole=="bothL") {
    w_l <- floor((window-1)/2)
    w_r <- window-1-w_l
  }
  if (hole=="bothR") {
    w_l <- floor((window-1)/2+0.5)
    w_r <- window-1-w_l
  }
  if (hole=="both") {
    w_l <- floor((window-1)/2+sample(c(0, 0.5), 1))
    w_r <- window-1-w_l
  }
  if (substr(hole, 1, 4)=="both") {
    if (fill) {
      y <- c(rep(y[window], w_l), y[window:length(y)], rep(y[length(y)], w_r))
    } else {
      y <- c(rep(NA, w_l), y[window:length(y)], rep(NA, w_r))
    }
  }
  return(y)
}
