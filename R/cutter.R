#' cutter returns the fitted distribution without cut
#' @title Distribution of the fitted distribution without cut. 
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The parameters of distribution of values below or above the detection limit.
#' @param observations The observations; see description
#' @param observations.colname If observations is a data.frame, the name of column with observations
#' @param lower_detection_limit.colname If observations is a data.frame, the name of column with lower detection limit
#' @param upper_detection_limit.colname If observations is a data.frame, the name of column with upper detection limit
#' @param cut_method.colname If observations is a data.frame, the name of column with cut method, being "censored" or "truncated"
#' @param par Initial values for parameters of distribution 
#' @param lower_detection_limit Value for lower detection limit 
#' @param upper_detection_limit Value for upper detection limit
#' @param cut_method Value for cut method, being "censored" or "truncated"
#' @param distribution Can be gamma, normal, weibull, lognormal, or generalized.gamma
#' @param n.iter Number of iteration for Bayesian MCMC and to estimate the goodness-of-fit
#' @param n.mixture Number of distributions
#' @param debug If TRUE, show some information
#' @param progress.bar If TRUE, show a progress bar for MCMC
#' @param session The session of a shiny process
#' @description If observations is a data.frame, it can have 4 columns:\cr
#' A column for the measurements;\cr
#' A column for the lower detection limit;\cr
#' A column for the upper detection limit;\cr
#' A column for the truncated of censored nature of the data.\cr
#' The names of the different columns are in the observations.colname, lower_detection_limit.colname, 
#' upper_detection_limit.colname and cut_method.colname.\cr
#' If lower_detection_limit.colname is NULL or if the column does not exist, 
#' the data are supposed to not be left-cut 
#' and if upper_detection_limit.colname is NULL or if the column does not exist, 
#' the data are supposed to not be right-cut.\cr
#' If observations is a vector, then the parameters lower_detection_limit and/or upper_detection_limit 
#' must be given. Then cut_method must be also provided.\cr
#' In abservations, -Inf must be used to indicate a value below the lower detection limit and +Inf must be used 
#' for a value above the upper detection limit.\cr
#' Be careful: NA is used to represent a missing data and not a value below of above the detection limit.\cr
#' If lower_detection_limit, upper_detection_limit or cut_method are only one value, they are supposed to 
#' be used for all the observations.\cr
#' Definitions for censored or truncated distribution vary, and the two terms are sometimes used interchangeably. Let the following data set:\cr
#' 1 1.25 2 4 5\cr
#' \cr
#' Censoring: some observations will be censored, meaning that we only know that they are below (or above) some bound. This can for instance occur if we measure the concentration of a chemical in a water sample. If the concentration is too low, the laboratory equipment cannot detect the presence of the chemical. It may still be present though, so we only know that the concentration is below the laboratory's detection limit.\cr
#' \cr
#' If the detection limit is 1.5, so that observations that fall below this limit is censored, our example data set would become:\cr
#' <1.5 <1.5 2 4 5; that is, we don't know the actual values of the first two observations, but only that they are smaller than 1.5.\cr
#' \cr
#' Truncation: the process generating the data is such that it only is possible to observe outcomes above (or below) the truncation limit. This can for instance occur if measurements are taken using a detector which only is activated if the signals it detects are above a certain limit. There may be lots of weak incoming signals, but we can never tell using this detector.\cr
#' \cr
#' If the truncation limit is 1.5, our example data set would become: 2 4 5; and we would not know that there in fact were two signals which were not recorded.\cr
#' If n.iter is NULL, no Bayesian MCMC is performed but credible interval will not be available.
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
#' plot(result, xlim=c(0, 150), breaks=seq(from=0, to=150, by=10), col.mcmc=NULL)
#' plot(result, xlim=c(0, 150), breaks=seq(from=0, to=150, by=10))
#' # _______________________________________________________________
#' # The same data seen as truncated data with gamma distribution
#' # _______________________________________________________________
#' obc <- obc[is.finite(obc)]
#' # search for the parameters the best fit these truncated data
#' result <- cutter(observations=obc, upper_detection_limit=DL, 
#'                            cut_method="truncated")
#' result
#' plot(result, xlim=c(0, 150), breaks=seq(from=0, to=150, by=10))
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
#' plot(result, xlim=c(0, 200), breaks=seq(from=0, to=300, by=10))
#' # _______________________________________________________________
#' # left censored distribution with mixture of gamma distribution
#' # _______________________________________________________________
#' #' # Detection limit
#' library(HelpersMG)
#' # Generate 200 random data from a gamma distribution
#' set.seed(1234)
#' obc <- c(rgamma(100, scale=10, shape=5), rgamma(100, scale=20, shape=10))
#' LDL <- 20
#' l <- seq(from=0, to=LDL, length.out=1001)
#' p <- pgamma(l, scale=10, shape=5)*0.5+pgamma(l, scale=20, shape=10)
#' deltal <- l[2]-l[1]
#' expected_LDL <- sum((l[-1]-deltal/2)*(p[-1]-p[-length(p)]))/sum((p[-1]-p[-length(p)]))
#' # remove the data below the detection limit
#' obc[obc<LDL] <- -Inf
#' 
#' UDL <- 300
#' l <- seq(from=UDL, to=1000, length.out=1001)
#' p <- pgamma(l, scale=10, shape=5)*0.5+pgamma(l, scale=20, shape=10)
#' deltal <- l[2]-l[1]
#' expected_UDL <- sum((l[-1]-deltal/2)*(p[-1]-p[-length(p)]))/sum((p[-1]-p[-length(p)]))
#' obc[obc>UDL] <- +Inf
#' 
#' # search for the parameters the best fit these truncated data
#' result1_gamma <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                           upper_detection_limit = UDL, 
#'                           distribution="gamma", 
#'                           cut_method="censored", n.iter=5000, debug=0)
#' result1_normal <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                           upper_detection_limit = UDL, 
#'                           distribution="normal", 
#'                           cut_method="censored", n.iter=5000)
#' result1_lognormal <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                           upper_detection_limit = UDL, 
#'                           distribution="lognormal", 
#'                           cut_method="censored", n.iter=5000)
#' result1_Weibull <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                           upper_detection_limit = UDL,  
#'                           distribution="Weibull", 
#'                           cut_method="censored", n.iter=5000)
#' result1_generalized.gamma <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                           upper_detection_limit = UDL, 
#'                           distribution="generalized.gamma", 
#'                           cut_method="censored", n.iter=5000)
#' result2_gamma <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                           upper_detection_limit = UDL, 
#'                           distribution="gamma", 
#'                           n.mixture=2, 
#'                           cut_method="censored", n.iter=5000)
#' result2_normal <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                           upper_detection_limit = UDL, 
#'                           distribution="normal", 
#'                           n.mixture=2, 
#'                           cut_method="censored", n.iter=5000)
#' result2_lognormal <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                           upper_detection_limit = UDL, 
#'                           distribution="lognormal", 
#'                           n.mixture=2, 
#'                           cut_method="censored", n.iter=5000)
#' result2_Weibull <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                           upper_detection_limit = UDL, 
#'                           distribution="Weibull", 
#'                           n.mixture=2, 
#'                           cut_method="censored", n.iter=5000)
#' result2_generalized.gamma <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                           upper_detection_limit = UDL, 
#'                           distribution="generalized.gamma", 
#'                           n.mixture=2, 
#'                           cut_method="censored", n.iter=5000)
#'                           
#' compare_AIC(nomixture.gamma=result1_gamma, 
#'              nomixture.normal=result1_normal, 
#'              nomixture.lognormal=result1_lognormal, 
#'              nomixture.Weibull=result1_Weibull, 
#'              nomixture.generalized.gamma=result1_generalized.gamma, 
#'              mixture.gamma=result2_gamma, 
#'              mixture.normal=result2_normal, 
#'              mixture.lognormal=result2_lognormal, 
#'              mixture.Weibull=result2_Weibull, 
#'              mixture.generalized.gamma=result2_generalized.gamma)
#'              
#' plot(result2_gamma, xlim=c(0, 600), breaks=seq(from=0, to=600, by=10))
#' plot(result2_generalized.gamma, xlim=c(0, 600), breaks=seq(from=0, to=600, by=10))
#' 
#' # _______________________________________________________________
#' # left and right censored distribution
#' # _______________________________________________________________
#' # Generate 100 random data from a gamma distribution
#' obc <- rgamma(100, scale=20, shape=2)
#' # Detection limit
#' LDL <- 10
#' # remove the data below the detection limit
#' obc[obc<LDL] <- -Inf
#' # Detection limit
#' UDL <- 100
#' # remove the data below the detection limit
#' obc[obc>UDL] <- +Inf
#' # search for the parameters the best fit these censored data
#' result <- cutter(observations=obc, lower_detection_limit=LDL, 
#'                            upper_detection_limit=UDL, 
#'                           cut_method="censored")
#' result
#' plot(result, xlim=c(0, 150), col.DL=c("black", "grey"), 
#'                              col.unobserved=c("green", "blue"), 
#'      breaks=seq(from=0, to=150, by=10))
#' # _______________________________________________________________
#' # Example with two values for lower detection limits
#' # corresponding at two different methods of detection for example
#' # with gamma distribution
#' # _______________________________________________________________
#' obc <- rgamma(50, scale=20, shape=2)
#' # Detection limit for sample 1 to 50
#' LDL1 <- 10
#' # remove the data below the detection limit
#' obc[obc<LDL1] <- -Inf
#' obc2 <- rgamma(50, scale=20, shape=2)
#' # Detection limit for sample 1 to 50
#' LDL2 <- 20
#' # remove the data below the detection limit
#' obc2[obc2<LDL2] <- -Inf
#' obc <- c(obc, obc2)
#' # search for the parameters the best fit these censored data
#' result <- cutter(observations=obc, 
#'                            lower_detection_limit=c(rep(LDL1, 50), rep(LDL2, 50)), 
#'                           cut_method="censored")
#' result
#' # It is difficult to choose the best set of colors
#' plot(result, xlim=c(0, 150), col.dist="red", 
#'      col.unobserved=c(rgb(red=1, green=0, blue=0, alpha=0.1), 
#'                       rgb(red=1, green=0, blue=0, alpha=0.2)), 
#'      col.DL=c(rgb(red=0, green=0, blue=1, alpha=0.5), 
#'                       rgb(red=0, green=0, blue=1, alpha=0.9)), 
#'      breaks=seq(from=0, to=200, by=10))
#'      
#' # ________________________________________________________________________
#' # left censored distribution comparison of normal, lognormal, 
#' # weibull, generalized gamma, and gamma without Bayesian MCMC
#' # Comparison with Akaike Information Criterion
#' # ________________________________________________________________________
#' # Detection limit
#' DL <- 10
#' # Generate 100 random data from a gamma distribution
#' obc <- rgamma(100, scale=20, shape=2)
#' # remove the data below the detection limit
#' obc[obc<DL] <- -Inf
#' 
#' result_gamma <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored", distribution="gamma", 
#'                           n.iter=NULL)
#' plot(result_gamma)
#' 
#' result_lognormal <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored", distribution="lognormal", 
#'                           n.iter=NULL)
#' plot(result_lognormal)
#' 
#' result_weibull <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored", distribution="weibull", 
#'                           n.iter=NULL)
#' plot(result_weibull)
#' 
#' result_normal <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored", distribution="normal", 
#'                           n.iter=NULL)
#' plot(result_normal)
#' 
#' result_generalized.gamma <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored", distribution="generalized.gamma", 
#'                           n.iter=NULL)
#' plot(result_generalized.gamma)
#' 
#' compare_AIC(gamma=result_gamma, 
#'             lognormal=result_lognormal, 
#'             normal=result_normal, 
#'             Weibull=result_weibull, 
#'             Generalized.gamma=result_generalized.gamma)
#' 
#' # ______________________________________________________________________________
#' # left censored distribution comparison of normal, lognormal, 
#' # weibull, generalized gamma, and gamma
#' # ______________________________________________________________________________
#' # Detection limit
#' DL <- 10
#' # Generate 100 random data from a gamma distribution
#' obc <- rgamma(100, scale=20, shape=2)
#' # remove the data below the detection limit
#' obc[obc<DL] <- -Inf
#' # search for the parameters the best fit these truncated data
#' result_gamma <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored", distribution="gamma")
#' result_gamma
#' plot(result_gamma, xlim=c(0, 250), breaks=seq(from=0, to=250, by=10))
#' 
#' result_lognormal <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored", distribution="lognormal")
#' result_lognormal
#' plot(result_lognormal, xlim=c(0, 250), breaks=seq(from=0, to=250, by=10))
#' 
#' result_weibull <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored", distribution="weibull")
#' result_weibull
#' plot(result_weibull, xlim=c(0, 250), breaks=seq(from=0, to=250, by=10))
#' 
#' result_normal <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored", distribution="normal")
#' result_normal
#' plot(result_normal, xlim=c(0, 250), breaks=seq(from=0, to=250, by=10))
#' 
#' result_generalized.gamma <- cutter(observations=obc, lower_detection_limit=DL, 
#'                           cut_method="censored", distribution="generalized.gamma")
#' result_generalized.gamma
#' plot(result_generalized.gamma, xlim=c(0, 250), breaks=seq(from=0, to=250, by=10))
#' 
#' # ___________________________________________________________________
#' # Test for similarity in gamma left censored distribution between two
#' # datasets
#' # ___________________________________________________________________
#' obc1 <- rgamma(100, scale=20, shape=2)
#' # Detection limit for sample 1 to 50
#' LDL <- 10
#' # remove the data below the detection limit
#' obc1[obc1<LDL] <- -Inf
#' obc2 <- rgamma(100, scale=10, shape=2)
#' # remove the data below the detection limit
#' obc2[obc2<LDL] <- -Inf
#' # search for the parameters the best fit these censored data
#' result1 <- cutter(observations=obc1, 
#'                   distribution="gamma", 
#'                   lower_detection_limit=LDL, 
#'                   cut_method="censored", n.iter=NULL)
#' logLik(result1)
#' plot(result1, xlim=c(0, 200), 
#'      breaks=seq(from=0, to=200, by=10))
#' result2 <- cutter(observations=obc2, 
#'                   distribution="gamma", 
#'                   lower_detection_limit=LDL, 
#'                   cut_method="censored", n.iter=NULL)
#' logLik(result2)
#' plot(result2, xlim=c(0, 200), 
#'      breaks=seq(from=0, to=200, by=10))
#' result_totl <- cutter(observations=c(obc1, obc2), 
#'                       distribution="gamma", 
#'                       lower_detection_limit=LDL, 
#'                       cut_method="censored", n.iter=NULL)
#' logLik(result_totl)
#' plot(result_totl, xlim=c(0, 200), 
#'      breaks=seq(from=0, to=200, by=10))
#'      
#' compare_AIC(Separate=list(result1, result2), 
#'             Common=result_totl, factor.value=1)
#' compare_BIC(Separate=list(result1, result2), 
#'             Common=result_totl, factor.value=1)           
#' }
#' @export

cutter <- function(observations=stop("Observations must be provided"), 
                   observations.colname="Observations", 
                   lower_detection_limit.colname="LDL", 
                   upper_detection_limit.colname="UDL", 
                   cut_method.colname="Cut", 
                   par=NULL, 
                   lower_detection_limit=NULL, 
                   upper_detection_limit=NULL, 
                   cut_method="censored", 
                   distribution="gamma", 
                   n.mixture=1, 
                   n.iter=5000, debug=FALSE, progress.bar=TRUE, 
                   session=NULL) {
  
  # observations=NULL 
  # observations.colname="Observations" 
  # lower_detection_limit.colname="LDL" 
  # upper_detection_limit.colname="UDL" 
  # cut_method.colname="Cut" 
  # par=NULL 
  # lower_detection_limit=NULL 
  # upper_detection_limit=NULL 
  # cut_method="censored" 
  # n.iter=5000
  # n.mixture=1
  # distribution="gamma"
  # debug <- FALSE
  # progress.bar=TRUE
  # n.iter=5000
  # debug=FALSE
  # progress.bar=TRUE
  # session=NULL
  
  # observations=obc
  # lower_detection_limit=LDL 
  # upper_detection_limit = UDL 
  # distribution="generalized.gamma" 
  # n.mixture=2 
  # cut_method="truncated"
  # n.iter=5000
  # debug = 1
  
  # observations = Ce
  # lower_detection_limit = 0.001
  # cut_method = "censored"
  # n.iter = 5000
  # n.mixture = 2
  # distribution = "gamma"
  
  fitn <- function(par, observations=NULL, distribution="gamma", 
                   n.mixture=NULL, debug=FALSE, limits.lower=NULL, 
                   limits.upper=NULL, log = TRUE) {
    - dcutter(par=par, observations=observations, distribution=distribution, 
              n.mixture=n.mixture, debug=debug, limits.lower=limits.lower, 
              limits.upper=limits.upper, log = log)
  }
  
  getparcutter <- getFromNamespace(".getparcutter", ns="HelpersMG")
  
  distribution <- tolower(distribution)
  distribution <- match.arg(distribution, choices = c("gamma", "normal", "lognormal", "weibull", "generalized.gamma"))
  
  
  if (is.data.frame(observations)) {
    if (all(colnames(observations) != observations.colname)) {
      stop("No observations are found in data.frame")
    }
    if (all(colnames(observations) != lower_detection_limit.colname)) {
      observations <- cbind(observations, LDL=NA)
      colnames(observations)[ncol(observations)] <- lower_detection_limit.colname
    }
    if (all(colnames(observations) != upper_detection_limit.colname)) {
      observations <- cbind(observations, UDL=NA)
      colnames(observations)[ncol(observations)] <- upper_detection_limit.colname
    }
    if (all(colnames(observations) != cut_method.colname)) {
      observations <- cbind(observations, Cut="censored")
      colnames(observations)[ncol(observations)] <- cut_method.colname
    }
    observations <- observations[, c(observations.colname, lower_detection_limit.colname, 
                                     upper_detection_limit.colname, cut_method.colname)]
    colnames(observations) <- c("Observations", "LDL", "UDL", "Cut")
  } else {
    observations <- data.frame(Observations=observations, 
                               LDL=NA, UDL=NA, Cut=NA)
    if (!is.null(lower_detection_limit)) observations$LDL <- lower_detection_limit
    if (!is.null(upper_detection_limit)) observations$UDL <- upper_detection_limit
    if (!is.null(cut_method)) observations$Cut <- cut_method
  }
  
  observations$Cut <- tolower(observations$Cut)
  
  if (any((observations$Cut == "truncated") & (is.infinite(observations$Observations)))) {
    stop("When data are truncated, all values must be finite")
  }
  
  obc <- observations$Observations[is.finite(observations$Observations)]
  qdistr <- switch(distribution, 
                   gamma=qgamma, 
                   lognormal=qlnorm, 
                   normal=qnorm, 
                   weibull=qweibull, 
                   generalized.gamma=qggamma)
  
  pdistr <- switch(EXPR = distribution, 
                   gamma=pgamma, 
                   lognormal=plnorm, 
                   normal=pnorm, 
                   weibull=pweibull, 
                   generalized.gamma=pggamma)
  
  ddistr <- switch(EXPR = distribution, 
                   gamma=dgamma, 
                   lognormal=dlnorm, 
                   normal=dnorm, 
                   weibull=dweibull, 
                   generalized.gamma=dggamma)
  
  if (is.null(par)) {
    
    mean <- mean(obc)
    var <- var(obc)
    if (distribution == "gamma") {
      s <- var/mean
      a <-  mean/s
      par <- c(shape=a, scale=s)
    }
    if (distribution == "lognormal") {
      mean <- mean(log(obc))
      var <- var(log(obc))
      par <- c(meanlog=mean, sdlog=sqrt(var))
    }
    if (distribution == "normal") {
      par <- c(mean=mean, sd=sqrt(var))
    }
    if (distribution == "weibull") {
      par <- c(shape=1, scale=10)
    }
    if (distribution == "generalized.gamma") {
      s <- var/mean
      a <-  mean/s
      par <- c(theta=s,  kappa=a, delta=1)
    }
    previouspar <- par
    names(par) <- paste0(names(par), "1")
    
    if (n.mixture > 1) {
      
      for (p in 2:n.mixture) {
        parint <- previouspar*0.95
        names(parint) <- paste0(names(parint), as.character(p))
        par <- c(par, parint)
      }
      parint <- rep(0.5, (n.mixture-1))
      names(parint) <- paste0("p", as.character(1:(n.mixture-1)))
      par <- c(par, parint)
    }
  }
  
  
  lower <- NULL
  upper <- NULL
  for (i in 1:length(par)) {
    parint <- par[i]
    nparint <- names(par[i])
    if (grepl("^p[0-9]", nparint)) {
      lowerint <- 1e-4
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- 100
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^shape[0-9]", nparint)) {
      lowerint <- parint/10
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- parint*10
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^scale[0-9]", nparint)) {
      lowerint <- parint/10
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- parint*10
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^delta[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^meanlog[0-9]", nparint)) {
      lowerint <- -abs(parint*10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^sdlog[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^mean[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^sd[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^theta[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^kappa[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
  }
  
  # fitn(par, observations=observations, distribution=distribution, n.mixture=n.mixture, debug=2)
  
  # first analysis by ML
  result2 <- optim(par, fitn, observations=observations, distribution=distribution, 
                   method="L-BFGS-B", 
                   limits.lower=lower, limits.upper=upper, 
                   lower=lower, 
                   upper=upper, 
                   hessian=FALSE, n.mixture=n.mixture, debug=debug)
  
  par <- result2$par
  
  if (!is.null(n.iter)) {
    
    par[par == lower] <- par[par == lower] * 2
    par[par == upper] <- par[par == upper] / 2
    
    lower <- NULL
    upper <- NULL
    for (i in 1:length(par)) {
      parint <- par[i]
      nparint <- names(par[i])
      if (grepl("^p[0-9]", nparint)) {
        lowerint <- 1e-4
        names(lowerint) <- nparint
        lower <- c(lower, lowerint)
        upperint <- 100
        names(upperint) <- nparint
        upper <- c(upper, upperint)
      }
      if (grepl("^shape[0-9]", nparint)) {
        lowerint <- parint/10
        names(lowerint) <- nparint
        lower <- c(lower, lowerint)
        upperint <- parint*10
        names(upperint) <- nparint
        upper <- c(upper, upperint)
      }
      if (grepl("^scale[0-9]", nparint)) {
        lowerint <- parint/10
        names(lowerint) <- nparint
        lower <- c(lower, lowerint)
        upperint <- parint*10
        names(upperint) <- nparint
        upper <- c(upper, upperint)
      }
      if (grepl("^delta[0-9]", nparint)) {
        lowerint <- abs(parint/10)
        names(lowerint) <- nparint
        lower <- c(lower, lowerint)
        upperint <- abs(parint*10)
        names(upperint) <- nparint
        upper <- c(upper, upperint)
      }
      if (grepl("^meanlog[0-9]", nparint)) {
        lowerint <- -abs(parint*10)
        names(lowerint) <- nparint
        lower <- c(lower, lowerint)
        upperint <- abs(parint*10)
        names(upperint) <- nparint
        upper <- c(upper, upperint)
      }
      if (grepl("^sdlog[0-9]", nparint)) {
        lowerint <- abs(parint/10)
        names(lowerint) <- nparint
        lower <- c(lower, lowerint)
        upperint <- abs(parint*10)
        names(upperint) <- nparint
        upper <- c(upper, upperint)
      }
      if (grepl("^mean[0-9]", nparint)) {
        lowerint <- abs(parint/10)
        names(lowerint) <- nparint
        lower <- c(lower, lowerint)
        upperint <- abs(parint*10)
        names(upperint) <- nparint
        upper <- c(upper, upperint)
      }
      if (grepl("^sd[0-9]", nparint)) {
        lowerint <- abs(parint/10)
        names(lowerint) <- nparint
        lower <- c(lower, lowerint)
        upperint <- abs(parint*10)
        names(upperint) <- nparint
        upper <- c(upper, upperint)
      }
      if (grepl("^theta[0-9]", nparint)) {
        lowerint <- abs(parint/10)
        names(lowerint) <- nparint
        lower <- c(lower, lowerint)
        upperint <- abs(parint*10)
        names(upperint) <- nparint
        upper <- c(upper, upperint)
      }
      if (grepl("^kappa[0-9]", nparint)) {
        lowerint <- abs(parint/10)
        names(lowerint) <- nparint
        lower <- c(lower, lowerint)
        upperint <- abs(parint*10)
        names(upperint) <- nparint
        upper <- c(upper, upperint)
      }
    }
    
    # prior <- data.frame(Density=rep('dnorm', length(par)), 
    #                     Prior1=par, 
    #                     Prior2=abs(par/4), 
    #                     SDProp=abs(log(abs(par))), 
    #                     Min=lower, 
    #                     Max=upper, 
    #                     Init=par, 
    #                     stringsAsFactors = FALSE, 
    #                     row.names=names(par))
    
    prior <- data.frame(Density=rep('dunif', length(par)), 
                        Prior1=lower, 
                        Prior2=upper, 
                        SDProp=abs(log(abs(par))), 
                        Min=lower, 
                        Max=upper, 
                        Init=par, 
                        stringsAsFactors = FALSE, 
                        row.names=names(par))
    
    if (!is.null(session)) {
      # Je suis en shiny
      mcmc_run <- MHalgoGen(n.iter=n.iter, parameters=prior, observations=observations, 
                            distribution=distribution, 
                            parameters_name="par", adaptive = TRUE, 
                            n.mixture=n.mixture, debug=debug, 
                            likelihood=fitn, n.chains=1, n.adapt=100, thin=1, trace=FALSE, 
                            progress.bar.ini=NULL, 
                            session=session, 
                            progress.bar=function(iter, session) {
                              getFromNamespace("updateProgressBar", ns="shinyWidgets")(
                                session = session,
                                id = "pb2",
                                value = iter, total = 10000,
                                title = "MCMC algorithm"
                              )
                            })
      
      
    } else {
      if (progress.bar) {
        
        pb <- txtProgressBar(min=0, max=n.iter*2, style=3)
        mcmc_run <- MHalgoGen(n.iter=n.iter, parameters=prior, observations=observations, 
                              distribution=distribution, 
                              parameters_name="par", adaptive = TRUE, 
                              n.mixture=n.mixture, debug=debug, 
                              likelihood=fitn, n.chains=1, n.adapt=100, thin=1, trace=FALSE, 
                              session=NULL, 
                              progress.bar.ini=NULL, 
                              progress.bar=function(iter, session=NULL) {setTxtProgressBar(get("pb", envir = parent.env(environment())), iter)})
        
      } else {
        mcmc_run <- MHalgoGen(n.iter=n.iter, parameters=prior, observations=observations, 
                              distribution=distribution, 
                              parameters_name="par", adaptive = TRUE, 
                              n.mixture=n.mixture, debug=debug, 
                              likelihood=fitn, n.chains=1, n.adapt=100, thin=1, trace=FALSE)
      }
    }
    par_mcmc <- as.quantiles(mcmc_run, namepar=names(par)[1], probs = c(0.025, 0.5, 0.975))
    
    for (i in seq_along(par)[-1])
      par_mcmc <- cbind(par_mcmc, 
                        as.quantiles(mcmc_run, namepar=names(par)[i], probs = c(0.025, 0.5, 0.975)))
    
    colnames(par_mcmc) <- names(par)
    
    result2 <- c(result2, list(mcmc=mcmc_run), list(par_mcmc=par_mcmc))
    
    parX <- suppressMessages(as.parameters(mcmc_run))
  } else {
    # result2 <- c(result2, list(mcmc=NA), list(par_mcmc=NA))
    parX <- par
  }
  
  par <- parX
  
  lower <- NULL
  upper <- NULL
  for (i in 1:length(par)) {
    parint <- par[i]
    nparint <- names(par[i])
    if (grepl("^p[0-9]", nparint)) {
      lowerint <- 1e-4
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- 100
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^shape[0-9]", nparint)) {
      lowerint <- parint/10
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- parint*10
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^scale[0-9]", nparint)) {
      lowerint <- parint/10
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- parint*10
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^delta[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^meanlog[0-9]", nparint)) {
      lowerint <- -abs(parint*10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^sdlog[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^mean[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^sd[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^theta[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
    if (grepl("^kappa[0-9]", nparint)) {
      lowerint <- abs(parint/10)
      names(lowerint) <- nparint
      lower <- c(lower, lowerint)
      upperint <- abs(parint*10)
      names(upperint) <- nparint
      upper <- c(upper, upperint)
    }
  }
  
  result3 <- optim(par, fitn, observations=observations, distribution=distribution, 
                   method="L-BFGS-B", 
                   limits.lower=lower, limits.upper=upper, 
                   lower=lower, 
                   upper=upper, n.mixture=n.mixture, 
                   hessian=TRUE, debug=debug)
  
  result2 <- modifyList(result2, result3)
  
  parX <- result2$par
  
  
  # Inverse the hessian matrix to get SE for each parameters
  mathessian <- result2$hessian
  inversemathessian <- try(solve(mathessian), silent = TRUE)
  if (inherits(inversemathessian, "try-error")) {
    res <- rep(NA, length(lower))
    names(res) <- colnames(mathessian)
  } else {
    if (any(diag(inversemathessian) < 0)) {
      res <- rep(NA, length(lower))
      names(res) <- colnames(mathessian)
    } else {
      res <- sqrt(diag(inversemathessian))
    }
  }
  result2 <- c(result2, list(distribution=distribution))
  result2 <- c(result2, list(observations=observations))
  result2 <- c(result2, list(n.mixture=n.mixture))
  result2 <- c(result2, list(n.LDL=sum(observations[, "Observations"] == c(-Inf))))
  result2 <- c(result2, list(n.UDL=sum(observations[, "Observations"] == c(Inf))))
  result2 <- c(result2, list(n.quantified=sum(is.finite(observations[, "Observations"]))))
  
  result2 <- c(result2, list(SE=res), list(AIC=2*result2$value+2*length(par)), 
               list(AICc=2*result2$value+4+(2*length(par)*(length(par)+1))/(nrow(observations)-length(par)-1)), 
               list(BIC=2*result2$value+log(nrow(observations))*length(par)))
  
  
  nlength <- 100
  
  mx <- max(obc[is.finite(obc)], na.rm = FALSE)
  for (p in 1:n.mixture) {
    mx_ec <- do.call(qdistr, args=c(list(p = 0.99999, lower.tail = TRUE, log.p = FALSE), as.list(getparcutter(parX, set=p))))
    if (is.finite(mx_ec)) {
      mx <- max(c(mx, mx_ec))
    }
  }
  
  seqx <- seq(from=0, to=mx, length.out=nlength)
  pparX <- getparcutter(parX, set=NULL)
  seqpg <- rep(0, length(nlength))
  for (p in 1:n.mixture) {
    seqpg <- seqpg + do.call(pdistr, args = modifyList(list(q=seqx, lower.tail = TRUE,
                                                            log.p = FALSE), as.list(getparcutter(parX, set=p))))*pparX[p]
  }
  
  meanx <- sum((seqpg[2:nlength]-seqpg[1:(nlength-1)])*(seqx[1:(nlength-1)]+(seqx[2]-seqx[1])/2))/(sum((seqpg[2:nlength]-seqpg[1:(nlength-1)])))
  sdx <- sqrt(sum(((seqx[1:(nlength-1)]+(seqx[2]-seqx[1])/2) - meanx)^2*(seqpg[2:nlength]-seqpg[1:(nlength-1)])))
  
  result2 <- c(result2, list(total.distribution=c(mean=unname(meanx), sd=unname(sdx))))
  
  attributes(result2) <- modifyList(attributes(result2), 
                                    list(nall=nrow(result2$observations), 
                                         df=length(result2$par), 
                                         nobs=nrow(result2$observations)))
  
  if (!is.null(n.iter)) {
    LDL_x <- unique(observations[, "LDL"])
    
    
    if (all(!is.na(LDL_x))) {
      prob_mcmc_LDL_tot <- NULL
      LDL_totp <- NULL
      
      for (LDL in LDL_x) {
        mcmc_xp <- NULL
        prob_mcmc_LDL <- NULL
        
        for (i in 1:n.iter) {
          nlength <- 100
          seqx <- seq(from=0, to=LDL, length.out=nlength)
          seqpg <- rep(0, length(seqx))
          pparX_mcmc <- getparcutter(mcmc_run$resultMCMC[[1]][i, ], set=NULL)
          for (p in 1:n.mixture)
            seqpg <- seqpg + do.call(pdistr, args = modifyList(list(q=seqx, lower.tail = TRUE,
                                                                    log.p = FALSE), as.list(getparcutter(mcmc_run$resultMCMC[[1]][i, ], set=p)))) * pparX_mcmc[p]
          
          mcmc_xp <- c(mcmc_xp, sum((seqpg[2:nlength]-seqpg[1:(nlength-1)])*(seqx[1:(nlength-1)]+(seqx[2]-seqx[1])/2))/(sum((seqpg[2:nlength]-seqpg[1:(nlength-1)]))))
          prob_mcmc_LDL <- c(prob_mcmc_LDL, do.call(pdistr, args = modifyList(list(q=LDL, lower.tail = TRUE,
                                                                                   log.p = FALSE), as.list(getparcutter(mcmc_run$resultMCMC[[1]][i, ], set=p)))) * pparX_mcmc[p])
        }
        LDL_totp <- c(LDL_totp, quantile(mcmc_xp, probs=c(0.025, 0.5, 0.975), na.rm=TRUE))
        prob_mcmc_LDL_tot <- c(prob_mcmc_LDL_tot, quantile(prob_mcmc_LDL, probs=c(0.025, 0.5, 0.975), na.rm=TRUE))
      }
      LDL_totp <- matrix(LDL_totp, nrow=3, byrow=FALSE)
      colnames(LDL_totp) <- as.character(LDL_x)
      rownames(LDL_totp) <- c("2.5%", "50%", "97.5%")
      result2 <- c(result2, list(LDL=LDL_totp))
      
      prob_mcmc_LDL_tot <- matrix(prob_mcmc_LDL_tot, nrow=3, byrow=FALSE)
      colnames(prob_mcmc_LDL_tot) <- as.character(LDL_x)
      rownames(prob_mcmc_LDL_tot) <- c("2.5%", "50%", "97.5%")
      result2 <- c(result2, list(prob.LDL=prob_mcmc_LDL_tot))
    }
    
    
    
    UDL_x <- unique(observations[, "UDL"])
    
    if (all(!is.na(UDL_x))) {
      prob_mcmc_UDL_tot <- NULL
      UDL_totp <- NULL
      
      for (UDL in UDL_x) {
        # mcmc_x <- NULL
        mcmc_xp <- NULL
        prob_mcmc_UDL <- NULL
        for (i in 1:n.iter) {
          # Alternative plus correcte
          nlength <- 100
          seqx <- seq(from=UDL, to=mx, length.out=nlength)
          seqpg <- rep(0, length(seqx))
          pparX_mcmc <- getparcutter(mcmc_run$resultMCMC[[1]][i, ], set=NULL)
          for (p in 1:n.mixture)
            seqpg <- seqpg + do.call(pdistr, args = modifyList(list(q=seqx, lower.tail = TRUE,
                                                                    log.p = FALSE), as.list(getparcutter(mcmc_run$resultMCMC[[1]][i, ], set=p)))) * pparX_mcmc[p]
          
          mcmc_xp <- c(mcmc_xp, sum((seqpg[2:nlength]-seqpg[1:(nlength-1)])*(seqx[1:(nlength-1)]+(seqx[2]-seqx[1])/2))/(sum((seqpg[2:nlength]-seqpg[1:(nlength-1)]))))
          
          prob_mcmc_UDL <- c(prob_mcmc_UDL, do.call(pdistr, args = modifyList(list(q=UDL, lower.tail = FALSE,
                                                                                   log.p = FALSE), as.list(getparcutter(mcmc_run$resultMCMC[[1]][i, ], set=p)))) * pparX_mcmc[p])
        }
        UDL_totp <- c(UDL_totp, quantile(mcmc_xp, probs=c(0.025, 0.5, 0.975), na.rm=TRUE))
        prob_mcmc_UDL_tot <- c(prob_mcmc_UDL_tot, quantile(prob_mcmc_UDL, probs=c(0.025, 0.5, 0.975), na.rm=TRUE))
        
      }
      UDL_totp <- matrix(UDL_totp, nrow=3, byrow=FALSE)
      colnames(UDL_totp) <- as.character(UDL_x)
      rownames(UDL_totp) <- c("2.5%", "50%", "97.5%")
      result2 <- c(result2, list(UDL=UDL_totp))
      
      prob_mcmc_UDL_tot <- matrix(prob_mcmc_UDL_tot, nrow=3, byrow=FALSE)
      colnames(prob_mcmc_UDL_tot) <- as.character(UDL_x)
      rownames(prob_mcmc_UDL_tot) <- c("2.5%", "50%", "97.5%")
      result2 <- c(result2, list(prob.UDL=prob_mcmc_UDL_tot))
    }
    
    
    
    L_MCMC <- do.call(fitn, args = list(observations=observations, 
                                        distribution=distribution, 
                                        par=result2$par_mcmc[2, ]))
    
    # Goodness of fit using likelihood of posterior predictive distribution
    
    L_PPD <- NULL
    obc <- observations
    for (i in 1:n.iter) {
      
      if (!is.null(session)) {
        getFromNamespace(x="updateProgressBar", ns="shinyWidgets")(
          session = session,
          id = "pb2",
          value = n.iter+i, total = 10000,
          title = "Posterior predictive distribution")
        
      } else {
        if (progress.bar) setTxtProgressBar(pb, n.iter+i)
      }
      # print(i)
      r <- rcutter(cutter=result2, 
                   observed_detection_limit=TRUE, 
                   random_method = "MCMC", 
                   index_mcmc = i)
      if (all(!is.na(r)) & any(is.finite(r))) {
        obc[!is.na(obc[, "Observations"]), "Observations"] <- r
        
        L <- do.call(fitn, args = list(observations=obc, 
                                       distribution=distribution, 
                                       par=result2$par_mcmc[2, ]))
        
        
        L_PPD <- c(L_PPD, L)
      }
    }
    
    result2 <- c(result2, list(GoF=sum(L_PPD>L_MCMC)/length(L_PPD)))
    result2 <- c(result2, list(n.GoF=length(L_PPD)))
    
  } else {
    # result2 <- c(result2, list(UDL=NA, LDL=NA))
    # result2 <- c(result2, list(GoF=NA))
    
    LDL_x <- unique(observations[, "LDL"])
    
    if (all(!is.na(LDL_x))) {
      LDL_totp <- NULL
      prob_mcmc_LDL_tot <- NULL
      for (LDL in LDL_x) {
        nlength <- 100
        seqx <- seq(from=0, to=LDL, length.out=nlength)
        seqpg <- rep(0, length(seqx))
        for (p in 1:n.mixture)
          seqpg <- seqpg + do.call(pdistr, args = modifyList(list(q=seqx, lower.tail = TRUE,
                                                                  log.p = FALSE), as.list(getparcutter(parX, set=p)))) * pparX[p]
        
        xp <- sum((seqpg[2:nlength]-seqpg[1:(nlength-1)])*(seqx[1:(nlength-1)]+(seqx[2]-seqx[1])/2))/(sum((seqpg[2:nlength]-seqpg[1:(nlength-1)])))
        LDL_totp <- c(LDL_totp, "mean"=xp)
        prob_mcmc_LDL_tot <- c(prob_mcmc_LDL_tot, do.call(pdistr, args = modifyList(list(q=LDL, lower.tail = TRUE,
                                                                                         log.p = FALSE), as.list(getparcutter(parX, set=p)))) * pparX[p])
        
      }
      LDL_totp <- matrix(LDL_totp, nrow=1, byrow=FALSE)
      colnames(LDL_totp) <- as.character(LDL_x)
      rownames(LDL_totp) <- c("mean")
      result2 <- c(result2, list(LDL=LDL_totp))
      
      prob_mcmc_LDL_tot <- matrix(prob_mcmc_LDL_tot, nrow=1, byrow=FALSE)
      colnames(prob_mcmc_LDL_tot) <- as.character(LDL_x)
      rownames(prob_mcmc_LDL_tot) <- c("mean")
      result2 <- c(result2, list(prob.LDL=prob_mcmc_LDL_tot))
      
    } 
    
    
    
    UDL_x <- unique(observations[, "UDL"])
    
    if (all(!is.na(UDL_x))) {
      UDL_totp <- NULL
      prob_mcmc_UDL_tot <- NULL
      
      for (UDL in UDL_x) {
        # Alternative plus correcte
        nlength <- 100
        seqx <- seq(from=UDL, to=mx, length.out=nlength)
        
        seqpg <- rep(0, length(seqx))
        for (p in 1:n.mixture)
          seqpg <- seqpg + do.call(pdistr, args = modifyList(list(q=seqx, lower.tail = TRUE,
                                                                  log.p = FALSE), as.list(getparcutter(parX, set=p)))) * pparX[p]
        
        xp <- sum((seqpg[2:nlength]-seqpg[1:(nlength-1)])*(seqx[1:(nlength-1)]+(seqx[2]-seqx[1])/2))/(sum((seqpg[2:nlength]-seqpg[1:(nlength-1)])))
        UDL_totp <- c(UDL_totp, "mean"=xp)
        
        prob_mcmc_UDL_tot <- c(prob_mcmc_UDL_tot, do.call(pdistr, args = modifyList(list(q=UDL, lower.tail = FALSE,
                                                                                         log.p = FALSE), as.list(getparcutter(parX, set=p)))) * pparX[p])
        
      }
      UDL_totp <- matrix(UDL_totp, nrow=1, byrow=FALSE)
      colnames(UDL_totp) <- as.character(UDL_x)
      rownames(UDL_totp) <- c("mean")
      result2 <- c(result2, list(UDL=UDL_totp))
      
      prob_mcmc_UDL_tot <- matrix(prob_mcmc_UDL_tot, nrow=1, byrow=FALSE)
      colnames(prob_mcmc_UDL_tot) <- as.character(UDL_x)
      rownames(prob_mcmc_UDL_tot) <- c("mean")
      result2 <- c(result2, list(prob.UDL=prob_mcmc_UDL_tot))
    }
    
  }
  
  
  if (progress.bar & is.null(session)) cat("\n")
  
  result2 <- addS3Class(result2, class=c("cutter", "CompareAIC"))
  return(result2)
}

.getparcutter <- function(par, set=1) {
  # Retourne le set de paramètres set
  # Ou si c'est NULL, le set des p classé
  par_ec <- unlist(par)
  if (is.null(set)) {
    parX <- par_ec[grepl("^p[0-9]*", names(par_ec))]
    parX <- parX[order(as.numeric(gsub(".+([0-9]+)", "\\1", names(parX))))]
    parX <- c(parX, structure(1, names=paste0("p", as.character(length(parX)+1))))
    parX <- parX/sum(parX)
  } else {
    parX <- par_ec[grepl(as.character(set), names(par_ec)) & names(par_ec) != paste0("p", as.character(set))]
    names(parX) <- gsub("[0-9]+", "", names(parX))
  }
  if (is.list(par)) parX <- as.list(parX)
  return(parX)
}



