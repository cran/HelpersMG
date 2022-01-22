#' rcutter returns random values based on fitted distribution with cut.
#' @title Random values of unobserved values of cut distribution. 
#' @author Marc Girondot
#' @return A vector with the random numbers.
#' @param cutter The fitted model obtained with cutter()
#' @param lower_detection_limit The lower detection limit
#' @param upper_detection_limit The upper detection limit
#' @param method_cut What method is used to cut the distribution: "censored", "truncated"?
#' @param n number of random numbers
#' @param observed_detection_limit If TRUE, will use the pattern of detection limit as in observations 
#' @param random_method How to get parameters; it can be "ML", "medianMCMC", or "MCMC"
#' @param index_mcmc For MCMC random_method, the index of data to be used.
#' @description Return n random numbers.\cr
#' It can be used to get the posterior predictive distribution; see example.\cr
#' If random_method is "ML", the parameter values obtained using maximum likelihood are used.\cr
#' If random_method is "medianMCMC", the parameter values obtained using median of posterior distribution are used.\cr
#' If random_method is "MCMC", the parameter values are one sample of the MCMC posterior distribution.\cr
#' if observed_detection_limit is set to TRUE, the number of random number is equal to the number of observations; n is not used.\cr
#' rcutter is the abbreviation for random-cutter.
#' @family Distributions
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' # _______________________________________________________________
#' # right censored distribution with gamma distribution
#' # _______________________________________________________________
#' # Detection limit
#' DL <- 100
#' # Generate 100 random data from a gamma distribution
#' obc <- rgamma(100, scale=20, shape=2)
#' # remove the data below the detection limit
#' obc[obc>DL] <- +Inf
#' # search for the parameters the best fit these censored data
#' result <- cutter(observations=obc, upper_detection_limit=DL, 
#'                            cut_method="censored")
#' result
#' # Posterior predictive distribution
#' r <- rcutter(cutter=result, upper_detection_limit=DL, n=100)
#' hist(r)
#' # _______________________________________________________________
#' # left censored distribution with gamma distribution
#' # _______________________________________________________________
#' # Detection limit
#' DL <- 10
#' # Generate 100 random data from a gamma distribution
#' obc <- rgamma(100, scale=20, shape=2)
#' # remove the data below the detection limit
#' obc[obc<DL] <- -Inf
#' # search for the parameters the best fit these truncated data
#' result <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored")
#' result
#' plot(result, breaks=seq(from=0, to=200, by=10))
#' r <- rcutter(cutter=result, n=100)
#' hist(r, breaks=seq(from=0, to=200, by=10))
#' r <- rcutter(cutter=result, lower_detection_limit=DL, n=100)
#' hist(r, breaks=seq(from=0, to=250, by=10))
#' # With censored method, some values are replaced with +Inf or -Inf
#' any(is.infinite(r))
#' r <- rcutter(cutter=result, upper_detection_limit=DL, n=100, 
#'              method_cut="truncated")
#' # With truncated method, the values below LDL or upper UDL are not present
#' any(is.infinite(r))
#' hist(r, breaks=seq(from=0, to=10, by=0.25))
#' r <- rcutter(cutter=result, observed_detection_limit=TRUE)
#' hist(r, breaks=seq(from=0, to=300, by=10))
#' }
#' @export

rcutter <- function(cutter=stop("A result of cutter() must be provided"), 
                    n=1, 
                    lower_detection_limit = NULL, 
                    upper_detection_limit = NULL, 
                    method_cut = c("censored", "truncated"), 
                    observed_detection_limit = FALSE, 
                    random_method = c("medianMCMC", "MCMC", "ML"), 
                    index_mcmc=NULL) {
  
  # cutter=result
  # n=1
  # lower_detection_limit = NULL
  # upper_detection_limit = NULL
  # method_cut = c("censored", "truncated")
  # observed_detection_limit = FALSE
  # random_method = c("ML", "medianMCMC", "MCMC")
  # index_mcmc=NULL
  
  n.mixture <- cutter$n.mixture
  
  getparcutter <- getFromNamespace(".getparcutter", ns="HelpersMG")
  
  random_method <- tolower(random_method)[1]
  random_method <- match.arg(random_method, choices = c("ml", "medianmcmc", "mcmc"))
  method_cut <- tolower(method_cut)[1]
  method_cut  <- match.arg(method_cut, choices = c("truncated", "censored"))
  distribution <- cutter$distribution
  
  if (random_method == "ml") {
    par <- cutter$par
  } else {
    if (random_method == "medianmcmc") {
      par <- cutter$par_mcmc[2, ]
    } else {
      if (is.null(index_mcmc)) {
        index_mcmc <- sample(x = 1:nrow(cutter$mcmc$resultMCMC[[1]]), size=1)
      }
      par <- cutter$mcmc$resultMCMC[[1]][index_mcmc,]
    }
  }
  
  obc <- cutter$observations[!is.na(cutter$observations[, "Observations"]), ]
  
  if (observed_detection_limit) {
    n <- nrow(obc)
    lower_detection_limit <- cutter$observations[!is.na(cutter$observations[, "Observations"]), "LDL"]
    upper_detection_limit <- cutter$observations[!is.na(cutter$observations[, "Observations"]), "UDL"]
    method_cut <- cutter$observations[!is.na(cutter$observations[, "Observations"]), "Cut"]
  } else {
    if (is.null(lower_detection_limit)) {
      lower_detection_limit <- rep(x=NA, n)
    } else {
      lower_detection_limit <- rep(x=lower_detection_limit, n)[1:n]
    }
    if (is.null(upper_detection_limit)) {
      upper_detection_limit <- rep(x=NA, n)
    } else {
      upper_detection_limit <- rep(x=upper_detection_limit, n)[1:n]
    }
    method_cut <- rep(x=method_cut, n)[1:n]
  }
  
  pparX <- getparcutter(par, set=NULL)
  
  rdistr <- switch(distribution, 
                   gamma = rgamma, 
                   lognormal = rlnorm, 
                   normal = rnorm, 
                   weibull = rweibull, 
                   generalized.gamma=rggamma)
  
  x2 <- NULL
  cpt <- 0
  # Je dois continuer si length(x2) < n et si cpt != 10*n
  while((length(x2) < n) & (cpt != (100*n))) {
    x1 <- NULL

    while((length(x1) < 1) & (cpt != (100*n))) {
      
      x1 <- NULL
      for (p in 1:n.mixture)
      x1 <- c(x1, do.call(rdistr, args = c(list(n = 1), as.list(getparcutter(par, set=p)))))
      
      
      x1 <- x1[which(runif(n=1) < cumsum(pparX))[1]]
      
      if (!is.na(lower_detection_limit[length(x2)+1])) {
        if (x1 < lower_detection_limit[length(x2)+1]) {
          if (method_cut[length(x2)+1] == "censored") {
            x1 <- -Inf
          } else {
            x1 <- NULL
          }
        }
      }
      if (!is.null(x1))
        if (!is.na(upper_detection_limit[length(x2)+1]) & is.finite(x1)) {
          if (x1 > upper_detection_limit[length(x2)+1]) {
            if (method_cut[length(x2)+1] == "censored") {
              x1 <- +Inf
            } else {
              x1 <- NULL
            }
          }
        }
      cpt <- cpt + 1
    }
    
    x2 <- c(x2, x1)
    
    if (cpt == n*100) {
      if (length(x2) == 0) {
        x2 <- NA
      }
      
        x2 <- rep(x2, n)[1:n]
      
    }
  }
  

  
  attributes(x2) <- c(attributes(x2), list(par=par))  
  
  return(x2)
}


