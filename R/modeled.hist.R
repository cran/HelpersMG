#' modeled.hist returns the theoretical value for the histogram bar based on a model of distribution.
#' @title Return the theoretical value for the histogram bar
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with x (the center of the bar) and y components
#' @param breaks Vector with the breaks; it can be obtained directly from hist()
#' @param FUN Function to be used to integrate the density, ex. pnorm
#' @param sum Total numbers in the histogram; 1 for emperical frequencies
#' @param ... Parameters to be used by FUN
#' @description Return the theoretical value for the histogram bar based on
#' a model of distribution.
#' @examples
#' \dontrun{
#' n <- rnorm(100, mean=10, sd=2)
#' breaks <- 0:20
#' hist(n, breaks=breaks)
#' 
#' s <- modeled.hist(breaks=breaks, FUN=pnorm, mean=10, sd=2, sum=100)
#' 
#' points(s$x, s$y, pch=19)
#' lines(s$x, s$y)
#' 
#' n <- rlnorm(100, meanlog=2, sdlog=0.4)
#' b <- hist(n, ylim=c(0, 70))
#' 
#' s <- modeled.hist(breaks=b$breaks, FUN=plnorm, meanlog=2, sdlog=0.4, sum=100)
#' 
#' points(s$x, s$y, pch=19)
#' lines(s$x, s$y)
#' }
#' @export

modeled.hist <- function(breaks, FUN, ..., sum=1) {
  # breaks is a vector with the breaks; it can be obtained directly from hist()
  # FUN is the function to integrate the density, ex. pnorm
  # ... are the parameters of FUN
  # sum is the total numbers in the histogram; 1 for emperical frequencies
  by <- breaks[2]-breaks[1]
  xp <- c(breaks, tail(breaks, 1)+by)
  par.fun <- list(...)
  par.fun <- modifyList(list(log.p=FALSE, q=xp, lower.tail=TRUE), 
                        par.fun)
  s <- do.call(FUN, par.fun)
  # s <- pfun(xp, par.fun[1], par.fun[2], log=FALSE)
  s <- utils::head(c(s[-1], 0)-s, length(breaks))
  return(list(x=breaks+by/2, y=s*sum))
}

