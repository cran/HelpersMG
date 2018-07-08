#' newdbeta returns the density for the Beta distributions
#' @title Density for the Beta distributions. 
#' @author Marc Girondot
#' @return newdbeta gives the density for the Beta distributions
#' @param x vector of quantiles.
#' @param mu mean of the Beta distribution.
#' @param v variance of the Beta distribution.
#' @param shape1 non-negative parameters of the Beta distribution.
#' @param shape2 non-negative parameters of the Beta distribution.
#' @param ncp non-centrality parameter.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @description Density for the Beta distribution with parameters mu and v or shape1 and 
#' shape2 (and optional non-centrality parameter ncp).
#' @details The Beta distribution with parameters shape1 = a and shape2 = b has density \cr
#' gamma(a+b)/(gamma(a)gamma(b))x^(a-1)(1-x)^(b-1)\cr
#' for a > 0, b > 0 and 0 <= x <= 1 where the boundary values at x=0 or x=1 are defined as by 
#' continuity (as limits).\cr
#' The mean is a/(a+b) and the variance is ab/((a+b)^2 (a+b+1)). These moments and all 
#' distributional properties can be defined as limits.
#' @examples
#' pi <- rbeta(100, shape1=0.48, shape2=0.12)
#' hist(pi, freq=FALSE, breaks=seq(from=0, to=1, by=0.1), ylim=c(0, 8), las=1)
#' library("HelpersMG")
#' mx <- ScalePreviousPlot()$ylim["end"]/
#'       max(newdbeta(seq(from=0.01, to=0.99, by=0.01), mu = 0.8, v=0.1))
#' curve(newdbeta(x, mu = 0.8, v=0.1)*mx, add=TRUE, col="red")
#' @export

newdbeta <- function(x, mu=NULL, v=NULL, shape1, shape2, ncp = 0, log = FALSE) {
	if (!is.null(mu) & !is.null(v)) {
  		shape1 <- mu*(((mu*(1-mu))/v)-1)
  		shape2 <- shape1*(1-mu)/mu
  	}
  dbeta(x=x, shape1=shape1, shape2 = shape2, ncp=ncp, log=log)
}
