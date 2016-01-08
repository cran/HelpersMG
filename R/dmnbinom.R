#' dmnbinom returns the density for the sum of random variable with negative binomial distributions
#' @title Density for the sum of random variable with negative binomial distributions. 
#' @author Marc Girondot
#' @references Furman, E., 2007. On the convolution of the negative binomial random variables. Statistics & Probability Letters 77, 169-172.
#' @return dmnbinom gives the density
#' @param x vector of (non-negative integer) quantiles.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param mu alternative parametrization via mean.
#' @param log	logical; if TRUE, probabilities p are given as log(p).
#' @param infinite Number of maximal iterations; check different values to determine the error in estimation
#' @description Density for the sum of random variable with negative binomial distributions.
#' @examples
#' alpha <- c(1, 2, 5, 1, 2)
#' p <- c(0.1, 0.12, 0.13, 0.14, 0.14)
#' # Test with lower iterations: 2 or 50 rather than 10 [default]; precision is very good still with 10
#' dmnbinom(20, size=alpha, prob=p, infinite=50)
#' dmnbinom(20, size=alpha, prob=p, infinite=10)
#' dmnbinom(20, size=alpha, prob=p, infinite=2)
#' # However it is not always the case; It depends on the parametrization (see Furman 2007)
#' dmnbinom(20, size=2, mu=c(0.01, 0.02, 0.03), infinite=1000)
#' dmnbinom(20, size=2, mu=c(0.01, 0.02, 0.03), infinite=100)
#' dmnbinom(20, size=2, mu=c(0.01, 0.02, 0.03), infinite=50)
#' dmnbinom(20, size=2, mu=c(0.01, 0.02, 0.03), infinite=10)
#' dmnbinom(20, size=2, mu=c(0.01, 0.02, 0.03), infinite=2)
#' # Test with a single distribution
#' dmnbinom(20, size=1, mu=20)
#' # when only one distribution is available, it is the same as dnbinom()
#' dnbinom(20, size=1, mu=20)
#' # If a parameter is supplied as only one value, it is supposed to be constant
#' dmnbinom(20, size=1, mu=c(14, 15, 10))
#' # The function is vectorized:
#' plot(0:200, dmnbinom(0:200, size=alpha, prob=p), bty="n", type="h", xlab="x", ylab="Density")
#' @export

dmnbinom <- function(x=stop("You must provide a x value"), 
                     size=stop("size parameter is mandatory"), 
                     prob=NULL, mu=NULL, log = FALSE, infinite=50) {
  
  # prob=NULL; mu=NULL; log = FALSE; infinite=10
  
  if (!is.null(mu)) prob <- size/(size+mu)
  
  alpha <- size
  p <- prob
  
  q <- 1-p
  p1 <- max(p)
  q1 <- 1-p1
  
  R <- prod(((q*p1)/(q1*p))^(-alpha))
  
  xi <- rep(NA, infinite)
  delta <- c(1, xi)
 
 for(i in 1:infinite) {
   xi[i] <- sum((alpha*(1-((q1*p)/(q*p1)))^i)/i)
   Ks <- 1:i
   delta[i+1] <- (1/i)*sum(Ks*xi[Ks]*delta[i-Ks+1])
}
  
  Pr <- sapply(x, function(S) {
  PrS <- R*sum(delta*dnbinom(S, size=sum(alpha)+seq_along(delta), prob=p1))
  if (log) PrS <- log(PrS)
  return(PrS)
  }
  )
  return(Pr)
}