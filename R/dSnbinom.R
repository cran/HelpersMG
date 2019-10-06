#' dSnbinom returns the density for the sum of random variable with negative binomial distributions
#' @title Density for the sum of random variable with negative binomial distributions. 
#' @author Marc Girondot
#' @references Furman, E., 2007. On the convolution of the negative binomial random variables. Statistics & Probability Letters 77, 169-172.
#' @return dSnbinom gives the density
#' @param x vector of (non-negative integer) quantiles.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param mu alternative parametrization via mean.
#' @param log	logical; if TRUE, probabilities p are given as log(p).
#' @param tol Tolerance for recurrence
#' @param infinite No more used - Maintained for retrocompatibility
#' @description Density for the sum of random variable with negative binomial distributions.\cr
#' If all prob values are the same, infinite is automatically set to 0.
#' @family Distribution of sum of random variable with negative binomial distributions
#' @examples
#' \dontrun{
#' alpha <- c(1, 2, 5, 1, 2)
#' p <- c(0.1, 0.12, 0.13, 0.14, 0.14)
#' dSnbinom(20, size=alpha, prob=p)
#' dSnbinom(20, size=2, mu=c(0.01, 0.02, 0.03))
#' # Test with a single distribution
#' dSnbinom(20, size=1, mu=20)
#' # when only one distribution is available, it is the same as dnbinom()
#' dnbinom(20, size=1, mu=20)
#' # If a parameter is supplied as only one value, it is supposed to be constant
#' dSnbinom(20, size=1, mu=c(14, 15, 10))
#' # The function is vectorized:
#' plot(0:200, dSnbinom(0:200, size=alpha, prob=p), bty="n", type="h", xlab="x", ylab="Density")
#' # Comparison with simulated distribution using rep replicates
#' alpha <- c(2.1, 2.05, 2)
#' mu <- c(10, 30, 20)
#' rep <- 100000
#' distEmpirique <- rSnbinom(rep, size=alpha, mu=mu)
#' tabledistEmpirique <- rep(0, 301)
#' names(tabledistEmpirique) <- as.character(0:300)
#' tabledistEmpirique[names(table(distEmpirique))] <- table(distEmpirique)/rep
#' 
#' plot(0:300, dSnbinom(0:300, size=alpha, mu=mu), type="h", bty="n", 
#'    xlab="x", ylab="Density", ylim=c(0,0.02))
#' plot_add(0:300, tabledistEmpirique, type="l", col="red")
#' legend(x=200, y=0.02, legend=c("Empirical", "Theoretical"), 
#'    text.col=c("red", "black"), bty="n")
#' 
#' # Example with the approximation mu=mean(mu)
#' plot(0:300, dSnbinom(0:300, size=alpha, mu=mu), type="h", bty="n", 
#'    xlab="x", ylab="Density", ylim=c(0,0.02))
#' plot_add(0:300, tabledistEmpirique, type="l", col="red")
#' legend(x=200, y=0.02, legend=c("Empirical", "Theoretical"), 
#'    text.col=c("red", "black"), bty="n")
#'    
#' # example to fit the distribution
#' data <- rnbinom(1000, size=1, mu=10)
#' hist(data)
#' ag <- rep(1:100, 10)
#' r <- aggregate(data, by=list(ag), FUN=sum)
#' hist(r[,2])
#' 
#' parx <- c(size=1, mu=10)
#' 
#' dSnbinomx <- function(x, par) {
#'   -sum(dSnbinom(x=x[,2], mu=rep(par["mu"], 10), size=par["size"], log=TRUE))
#' }
#' 
#' fit_mu_size <- optim(par = parx, fn=dSnbinomx, x=r, method="BFGS", control=c(trace=TRUE))
#' fit_mu_size$par
#' }
#' @export

dSnbinom <- function(x=stop("You must provide a x value"), 
                     size=NULL, 
                     prob=NULL, mu=NULL, log = FALSE,  tol=1E-6,
                     infinite=1000) {
  
  infinite <- 1000
  # prob=NULL; mu=NULL; log = FALSE; infinite=10
  if (is.null(mu) + is.null(size) + is.null(prob) != 1) stop("Two values among mu, size and prob must be provided")
  
  m <- max(c(length(size), length(prob), length(mu)))
  if (!is.null(mu)) mu <- rep(mu, m)[1:m]
  if (!is.null(size)) size <- rep(size, m)[1:m]
  if (!is.null(prob)) prob <- rep(prob, m)[1:m]
  
  if (is.null(prob)) {
    prob <- size/(size+mu)
    prob <- ifelse(prob>1-(1e-9), 1-1e-6, prob)
  }
  if (is.null(mu)) mu <- size/prob - size
  if (is.null(size))  size  <- (prob * mu) / (1 - prob)
  
  if ((infinite==0) | (all(prob==prob[1]))) {
    
#    if (length(prob)<length(size)) prob <- rep(prob, length(size))[1:length(size)]
#    if (length(size)<length(prob)) size <- rep(size, length(prob))[1:length(prob)]
    
    return(dnbinom(x, size=sum(size), prob=mean(prob), log=log))
    
  } else {
    
    
    alpha <- size
    p <- prob
    
    q <- 1-p
    p1 <- max(p)
    q1 <- 1-p1
    
    # R <- prod(((q*p1)/(q1*p))^(-alpha))
    R <- exp(sum(log(((q*p1)/(q1*p))^(-alpha))))
    if (R == 0) warning("Probability too low to be estimated")
    
    xi <- rep(NA, infinite)
    delta <- c(1, xi)
    
    for(i in 1:infinite) {
      xi[i] <- sum((alpha*(1-((q1*p)/(q*p1)))^i)/i)
      Ks <- 1:i
      delta[i+1] <- (1/i)*sum(Ks*xi[Ks]*delta[i-Ks+1])
      if (abs((delta[i+1] - delta[i])) < tol) {
        delta <- delta[!is.na(delta)]
        break
      }
    }
    
    
    
    Pr <- sapply(x, function(S) {
      PrS <- R*sum(delta*dnbinom(S, size=sum(alpha)+seq_along(delta)-1, prob=p1))
      if (log) PrS <- log(PrS)
      return(PrS)
    }
    )
    return(Pr)
  }
}
