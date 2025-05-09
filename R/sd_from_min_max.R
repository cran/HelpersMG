#' sd_from_min_max returns standard deviation and/or mean from minimum and maximum
#' @title Distribution from minimum and maximum
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return sd_from_min_max returns a mcmcComposite object
#' @param observedMinimum The observed minimum
#' @param observedMaximum The observed maximum
#' @param observedMean The observed mean (can be omited)
#' @param priors The priors (can be omited)
#' @param n Number of observations used to estimate observedMinimum and observedMaximum
#' @param fittedparameters Can be "sd" or "mean-sd"
#' @param n.iter Number of iterations for each chain
#' @param n.chains Number of chains
#' @param n.adapt Number of iteration to stabilize likelihood
#' @param thin Interval for thinning likelihoods
#' @param adaptive Should an adaptive process for SDProp be used
#' @param trace Or FALSE or period to show progress
#' @description Bayesian estimate of underlying distribution.\cr
#' The distribution of extrema is expected to be a Gaussian distribution; we could do better using 
#' generalized extreme value (GEV) distribution.\cr
#' Weisstein, Eric W. "Extreme value distribution". mathworld.wolfram.com
#' @examples
#' \dontrun{
#' minobs <- 45
#' maxobs <- 53
#' n <- 4
#' # To estimate only the sd of the distribution
#' out_sd_mcmc <- sd_from_min_max(n=n, observedMinimum=minobs, 
#'                            observedMaximum=maxobs, 
#'                            fittedparameters="sd")
#' plot(out_sd_mcmc, what="MarkovChain", parameters="sd")
#' plot(out_sd_mcmc, what="posterior", parameters="sd")
#' as.parameters(out_sd_mcmc, index="quantile")
#' # To be compared with the rule of thumb:
#' print(paste0("sd = ", as.character((maxobs - minobs) / 4))) # SD Clearly biased
#' 
#' # To estimate both the sd and mean of the distribution
#' 
#' out_sd_mean_mcmc <- sd_from_min_max(n=n, observedMinimum=minobs, 
#'                            observedMaximum=maxobs, 
#'                            fittedparameters="mean-sd")
#' 
#' plot(out_sd_mean_mcmc, what="MarkovChain", parameters="sd")
#' plot(out_sd_mean_mcmc, what="MarkovChain", parameters="mean")
#' plot(out_sd_mean_mcmc, what="posterior", parameters="mean", xlim=c(0, 100), 
#'      breaks=seq(from=0, to=100, by=5))
#' as.parameters(out_sd_mean_mcmc, index="quantile")
#' # To be compared with the rule of thumb:
#' print(paste0("mean = ", as.character((maxobs + minobs) / 2))) # Mean Not so bad
#' print(paste0("sd = ", as.character((maxobs - minobs) / 4))) # SD Clearly biased
#' 
#' # Covariation of sd and mean is nearly NULL
#' cor(x=as.parameters(out_sd_mean_mcmc, index="all")[, "mean"], 
#'     y=as.parameters(out_sd_mean_mcmc, index="all")[, "sd"])^2
#' plot(x=as.parameters(out_sd_mean_mcmc, index="all")[, "mean"], 
#'      y=as.parameters(out_sd_mean_mcmc, index="all")[, "sd"], 
#'      xlab="mean", ylab="sd")
#'
#' 
#' }
#' @export

sd_from_min_max <- function(n                                    , 
                            observedMean=NULL                    , 
                            observedMinimum                      , 
                            observedMaximum                      , 
                            priors=NULL                          , 
                            fittedparameters=c("mean-sd")        , 
                            n.iter=10000                         ,
                            n.chains = 1                         ,
                            n.adapt = 100                        ,
                            thin=30                              ,
                            adaptive = FALSE                     , 
                            trace = 100) {
  
  # https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution
  
  fittedparameters <- match.arg(fittedparameters, choices = c("sd", "mean-sd"))
  
  L_sd_min_max <- getFromNamespace(".L_sd_min_max", ns="HelpersMG")
  
  if (is.null(priors)) {
    if (fittedparameters == "sd") {
      priors <- setPriors(par=c(sd=(observedMaximum - observedMinimum) / 4), 
                          se=c(sd=1), 
                          density = "dunif", 
                          rules = data.frame(Name="sd", Min=1E-6, Max=10*(observedMaximum - observedMinimum) / 4))
    } else {
      priors <- setPriors(par=c(sd=(observedMaximum - observedMinimum) / 4, 
                                mean=(observedMaximum + observedMinimum) / 2), 
                          se=c(sd=1, mean=1), 
                          density = "dunif", 
                          rules = rbind(data.frame(Name="sd", Min=1E-6, Max=10*(observedMaximum - observedMinimum) / 4), 
                                        data.frame(Name="mean", Min=0, Max=10*(observedMaximum + observedMinimum) / 2)))
      
    }
    
  }
  
  out_sd_mean_mcmc <- MHalgoGen(likelihood = L_sd_min_max                                 , 
                                parameters = priors                                       , 
                                meanobs=observedMean                                      , 
                                n=n                                                       , 
                                minobs=observedMinimum                                    , 
                                maxobs=observedMaximum                                    , 
                                n.iter=n.iter                                               , 
                                n.chains = n.chains                                               , 
                                n.adapt = n.adapt                                              , 
                                thin=thin                                                    , 
                                trace=trace                                               )
  
  print(as.parameters(out_sd_mean_mcmc, index="quantile"))
  
  return(invisible(out_sd_mean_mcmc))
  
}

.L_sd_min_max <- function(x, n, meanobs=NULL, minobs, maxobs, rep=10000, dist_extrema=dnorm) {
  if (is.null(meanobs)) {
    meanobs <- (minobs + maxobs) / 2
  }
  sd <- x["sd"]
  if (any(names(x) == "mean")) meanobs <- x["mean"]
  dsit <- sapply(1:rep, FUN=function(i) {k <- (rnorm(n=n, mean=meanobs, sd=sd)); return(c(min(k), max(k)))})
  kk <- matrix(dsit, ncol=2, byrow = TRUE)
  mean_min <- mean(kk[, 1]); sd_min <- sd(kk[, 1])
  mean_max <- mean(kk[, 2]); sd_max <- sd(kk[, 2])
  # Make the hypothesis that distribution of maximum and minimum are Gaussian
  L <- - dist_extrema(minobs, mean_min, sd_min, log=TRUE) - dist_extrema(maxobs, mean_max, sd_max, log=TRUE)
  return(L)
}