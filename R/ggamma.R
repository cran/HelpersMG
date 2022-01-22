#' Generalized gamma distribution
#' 
#' \code{pggamma}, \code{qggamma}, \code{dggamma}, and \code{rggamma} are used to model the 
#' generalized gamma distribution.
#' 
#' The code is modified from \url{https://rpubs.com/FJRubio/GG}.
#' 
# #' @name ggamma
#' @title Generalized gamma distribution. 
#' @author Marc Girondot
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param lower.tail logical; if TRUE (default), probabilities are \emph{P[X <= x]} otherwise, \emph{P[X > x]}.
#' @param theta scale parameter.
#' @param kappa shape parameter.
#' @param delta shape parameter.
#' @param log,log.p logical; if TRUE, probabilities \emph{p} are given as \emph{log(p)}.
#' @return \code{dggamma} gives the density, \code{pggamma} gives the distribution function, 
#' \code{qggamma} gives the quantile function, and \code{rggamma} generates random deviates.
#' @family Distributions
#' @examples
#' # To reproduce the wikipedia page graphic
#' x <- seq(from=0, to=8, by=0.1)
#' plot(x, dggamma(x, theta=2, kappa=0.5, delta=0.5), lty=1, col="blue", 
#'      type="l", lwd=2, xlab="x", ylab="PDF")
#' lines(x, dggamma(x, theta=1, kappa=1, delta=0.5), lty=1, col="green", lwd=2)
#' lines(x, dggamma(x, theta=2, kappa=1, delta=2), lty=1, col="red", lwd=2)
#' lines(x, dggamma(x, theta=5, kappa=1, delta=5), lty=1, col="yellow", lwd=2)
#' lines(x, dggamma(x, theta=7, kappa=1, delta=7), lty=1, col="grey", lwd=2)
#' legend("topright", legend=c("a=2, d=0.5, p=0.5", "a=1, d=1, p=0.5", 
#'                             "a=2, d=1, p=2", "a=5, d=1, p=5", "a=7, d=1, p=7"), 
#'                             col=c("blue", "green", "red", "yellow", "grey"), 
#'                             lty=1, lwd=2, bty="n")
#' @export


#' @section More details here:
#' The generalized gamma is described here \url{https://en.wikipedia.org/wiki/Generalized_gamma_distribution}.\cr
#' With \code{a} being \code{theta}, \code{b} being \code{kappa}, and \code{p} being \code{delta}.\cr
#' \code{theta}, \code{kappa} and \code{delta} must be all > 0.

#' @describeIn ggamma Density of the generalized gamma.
# #' @section Another section after function section:
dggamma <- function(x, theta, kappa, delta, log = FALSE){
  val <- log(delta) - kappa*log(theta) - lgamma(kappa/delta) + (kappa - 1)*log(x) -
    (x/theta)^delta
  if(log) return(val) else return(exp(val))
}
#' @export
#' @describeIn ggamma Distribution function of the generalized gamma.
# #' @section Another section after function section:
pggamma <- function(q, theta, kappa, delta, lower.tail = TRUE, log.p = FALSE){
  val <- pgamma( q^delta, shape = kappa/delta, scale = theta^delta, log.p = TRUE) 
  if (!lower.tail) val <- log(1 - exp(val))
  if(log.p) return(val) else return(exp(val))
}
#' @export

#' @describeIn ggamma Quantile of the generalized gamma.
# #' @section Another section after function section:
qggamma <- function(p, theta, kappa, delta, lower.tail = TRUE,
                    log.p = FALSE){
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  out <- qgamma(p, shape = kappa/delta, scale = theta^delta)^(1/delta)
  return(out)
}
#' @export

#' @describeIn ggamma Random of the generalized gamma.
# #' @section Another section after function section:
rggamma <- function(n, theta, kappa, delta){
  p <- runif(n)
  out <- qgamma(p, shape = kappa/delta, scale = theta^delta)^(1/delta)
  return(as.vector(out))
}
#' @export
