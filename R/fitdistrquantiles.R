#' fitdistrquantiles returns the parameters of beta, normal or gamma distribution
#' @title Parameters of beta, normal or gamma distribution based on quantiles. 
#' @author Marc Girondot
#' @return Parameters of beta, normal or gamma distribution based on quantiles.
#' @param quantiles Vector of quantiles.
#' @param probs Numeric vector of probabilities with values in [0,1].
#' @param scaled Used scaled least-square.
#' @param distribution Distribution to be fitted: beta, normal, or gamma.
#' @description Return the parameters of beta or gamm that fits the best the
#' quantiles. The vector of probabilities can be obtained from names of 
#' quantiles.
#' @examples
#' rd <- rbeta(100000, shape1 = 0.7, shape2 = 6.2, ncp=0)
#' (q <- quantile(rd, probs=c(0.025, 0.5, 0.975)))
#' 
#' (best <- fitdistrquantiles(quantiles = q, probs = c(0.025, 0.5, 0.975), 
#'                            scaled=FALSE, distribution = "beta"))
#' rd10000 <- rbeta(10000, shape1 = best["shape1"], shape2 = best["shape2"], ncp=best["ncp"])
#' quantile(rd10000, probs=c(0.025, 0.5, 0.975))
#' 
#' # Here the probabilities are obtained from names of quantiles
#' (best <- fitdistrquantiles(quantiles = q, scaled=FALSE, distribution = "beta"))
#' rd10000 <- rbeta(10000, shape1 = best["shape1"], shape2 = best["shape2"], ncp=best["ncp"])
#' quantile(rd10000, probs=c(0.025, 0.5, 0.975))
#' 
#' # If only two quantiles are provided, ncp cannot be fitted
#' (q2 <- quantile(rd, probs=c(0.025, 0.975)))
#' (best <- fitdistrquantiles(quantiles = q2, scaled=FALSE, distribution = "beta"))
#' rd10000 <- rbeta(10000, shape1 = best["shape1"], shape2 = best["shape2"])
#' quantile(rd10000, probs=c(0.025, 0.975))
#' x <- seq(from=0.00, to=1, by=0.001)
#' plot(x=x, y=pbeta(x, shape1 = best["shape1"], shape2 = best["shape2"]), 
#'      las=1, bty="n", type="l", ylim=c(0, 1))
#' segments(x0=q2[1], x1=q2[1], y0=0, y1=1, lty=2)
#' segments(x0=q2[2], x1=q2[2], y0=0, y1=1, lty=2)
#' 
#' (best <- fitdistrquantiles(quantiles = q, probs = c(0.025, 0.5, 0.975), 
#'                            scaled=FALSE, distribution = "gamma"))
#' rd10000 <- rgamma(10000, shape = best["shape"], scale = best["scale"])
#' quantile(rd10000, probs=c(0.025, 0.5, 0.975))
#' 
#' (best <- fitdistrquantiles(quantiles = c(10, 20, 30), probs = c(0.025, 0.5, 0.975), 
#'                            scaled=FALSE, distribution = "normal"))
#' rd10000 <- rnorm(10000, mean = best["mean"], sd = best["sd"])
#' quantile(rd10000, probs=c(0.025, 0.5, 0.975))
#' @export

fitdistrquantiles <- function(quantiles=stop("At least two quantiles must be provided"), 
                              probs=NULL, 
                              scaled=FALSE, distribution="beta") {
  distribution <- match.arg(distribution, choices = c("beta", "gamma", "normal"))
  
  optimfitdistrquantiles <- function(x, quantiles, probs, scaled=FALSE, distribution="beta") {
    if (distribution == "beta") {
      shape1 <- x["shape1"]
      shape2 <- x["shape2"]
      ncp <- x["ncp"]
      if (is.na(ncp)) {
        out_distr <- suppressWarnings(qbeta(p=probs, shape1 = shape1, shape2 = shape2))
        
      } else {
        out_distr <- suppressWarnings(qbeta(p=probs, shape1 = shape1, shape2 = shape2, ncp=ncp))
      }
    }
    if (distribution =="gamma") {
      rate <- x["rate"]
      shape <- x["shape"]
      scale <- x["scale"]
      # if (is.na(rate)) {
      out_distr <- suppressWarnings(qgamma(p=probs, shape=shape, scale=scale))
      # } else {
      #   out_distr <- qgamma(p=probs, rate=rate, shape=shape)
      # }
    }
    if (distribution =="normal") {
      mean <- x["mean"]
      sd <- x["sd"]
      out_distr <- suppressWarnings(qnorm(p=probs, mean=mean, sd=sd))
    }
    
    if (scaled) {
      c <- sum(((out_distr-quantiles)/out_distr)^2)
    } else {
      c <- sum((out_distr-quantiles)^2)
    }
    c <- ifelse(is.na(c), +Inf, c)
    return(c)
  }
  
  if (is.null(probs) & is.null(names(quantiles))) {
    stop("Or probs or named quantiles must be provided.")
  }
  if (is.null(probs)) {
    probs <- as.numeric(gsub("%", "", names(quantiles)))/ifelse(grepl("%", names(quantiles)), 100, 1)
  }
  quantile_median <- unname(quantiles[which.min(abs(probs-0.5))][1])
  if (distribution == "beta") {
    # The mean is shape1/(shape1+shape2) and the variance is shape1*shape2/((shape1+shape2)^2 (shape1+shape2+1)). 
    if (length(quantiles)>2) {
      par <- c(shape1=0.5, shape2=(0.5-quantile_median*0.5)/quantile_median, ncp=0)
    } else {
      par <- c(shape1=0.5, shape2=(0.5-quantile_median*0.5)/quantile_median)
    }
    lower <- rep(1e-7, length(par))
  }
  if (distribution == "gamma") {
    # The mean and variance are E(X) = shape*scale and Var(X) = shape*scale^2.
    par <- c(shape=1, scale=1/quantile_median)
    lower <- rep(1e-7, length(par))
  }
  if (distribution == "normal") {
    par <- c(mean=quantile_median, sd=1)
    lower <- c(-Inf, 1e-7)
  }
  
  
  best <- optim(par=par, 
                lower=lower, upper=rep(Inf, length(par)), 
                method="L-BFGS-B", 
                fn=optimfitdistrquantiles, 
                quantiles=quantiles, probs=probs, scaled=scaled, distribution=distribution)
  return(best$par)
  
}

