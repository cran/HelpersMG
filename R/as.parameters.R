#' @title Extract parameters from mcmcComposite object
#' @author Marc Girondot
#' @return A vector with parameters at maximum likelihood or index position
#' @param x A mcmcComposite obtained as a result of \code{MHalgoGen()} function
#' @param index At which iteration the parameters must be taken, see description
#' @param probs Quantiles to be returned, see description
#' @param chain The number of the chain in which to get parameters
#' @description Take a mcmcComposite object and create a vector object with parameter value at specified iteration.\cr
#' If \code{index="best"}, the function will return the parameters for the highest likelihood. It also indicates at which iteration the maximum lihelihood has been observed.\cr
#' If \code{index="last"}, the function will return the parameters for the last likelihood.\cr
#' If \code{index="median"}, the function will return the median value of the parameter.\cr
#' if \code{index="quantile"}, the function will return the \code{probs} defined by quantiles parameter.\cr
#' If \code{index="mode"}, the function will return the mode value of the parameter based on Asselin de Beauville (1978) method.\cr
#' \code{index} can also be a numeric value.\cr
#' This function uses the complete iterations available except the adaptation part, 
#' even if thin parameter is not equal to 1.
#' @references Asselin de Beauville J.-P. (1978). Estimation non paramétrique de la densité et du mode, 
#' exemple de la distribution Gamma. Revue de Statistique Appliquée, 26(3):47-70.
#' @family mcmcComposite functions
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' require(coda)
#' x <- rnorm(30, 10, 2)
#' dnormx <- function(data, x) {
#'  data <- unlist(data)
#'  return(-sum(dnorm(data, mean=x['mean'], sd=x['sd'], log=TRUE)))
#' }
#' parameters_mcmc <- data.frame(Density=c('dnorm', 'dlnorm'), 
#'                               Prior1=c(10, 0.5), 
#'                               Prior2=c(2, 0.5), 
#'                               SDProp=c(1, 1), 
#'                               Min=c(-3, 0), 
#'                               Max=c(100, 10), 
#'                               Init=c(10, 2), 
#'                               stringsAsFactors = FALSE, 
#'                               row.names=c('mean', 'sd'))
#' mcmc_run <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, data=x, 
#'                       likelihood=dnormx, n.chains=1, n.adapt=100, 
#'                       thin=1, trace=1)
#' plot(mcmc_run, xlim=c(0, 20))
#' plot(mcmc_run, xlim=c(0, 10), parameters="sd")
#' mcmcforcoda <- as.mcmc(mcmc_run)
#' #' heidel.diag(mcmcforcoda)
#' raftery.diag(mcmcforcoda)
#' autocorr.diag(mcmcforcoda)
#' acf(mcmcforcoda[[1]][,"mean"], lag.max=20, bty="n", las=1)
#' acf(mcmcforcoda[[1]][,"sd"], lag.max=20, bty="n", las=1)
#' batchSE(mcmcforcoda, batchSize=100)
#' 
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmcforcoda)$statistics[,"Time-series SE"]
#' summary(mcmc_run)
#' as.parameters(mcmc_run)
#' lastp <- as.parameters(mcmc_run, index="last")
#' parameters_mcmc[,"Init"] <- lastp
#' 
#' # The n.adapt set to 1 is used to not record the first set of parameters
#' # then it is not duplicated (as it is also the last one for 
#' # the object mcmc_run)
#' mcmc_run2 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, data=x, 
#'                        likelihood=dnormx, n.chains=1, n.adapt=1, 
#'                        thin=1, trace=1)
#' mcmc_run3 <- merge(mcmc_run, mcmc_run2)
#' 
#' ####### no adaptation, n.adapt must be 0
#' parameters_mcmc[,"Init"] <- c(mean(x), sd(x))
#' mcmc_run3 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, data=x, 
#'                        likelihood=dnormx, n.chains=1, n.adapt=0, 
#'                        thin=1, trace=1)
#' # With index being median, it returns the median value for each parameter
#' as.parameters(mcmc_run3, index="median")
#' as.parameters(mcmc_run3, index="mode")
#' as.parameters(mcmc_run3, index="best")
#' as.parameters(mcmc_run3, index="quantile", probs=0.025)
#' as.parameters(mcmc_run3, index="quantile", probs=0.975)
#' as.parameters(mcmc_run3, index="quantile", probs=c(0.025, 0.975))
#' }
#' @export


as.parameters <- function(x, index="best", chain=1, probs=0.025) {
  
  p <- x$resultMCMC.total[[chain]]
  if (is.null(p)) p <- x$resultMCMC[[chain]]
  
  if (index == "median") {
    pml <- apply(p, MARGIN=2, FUN = median)
    names(pml) <- colnames(p)
  } else {
    if (index=="quantile") {
      pml <- apply(p, MARGIN=2, FUN = quantile, probs=probs)
      names(pml) <- colnames(p)
  } else {
    
    if (index == "mode") {
      asselin <- function (x, bw = NULL, ...) {
        if (is.null(bw)) 
          bw <- 1
        nx <- length(x)
        kmax <- floor(ifelse(nx < 30, 10, 15) * log(nx))
        y <- sort(x)
        ok1 <- FALSE
        while (!ok1) {
          ny <- length(y)
          if (ny == 1) 
            return(y)
          qy <- stats::quantile(y, probs = c(0.1, 0.25, 0.5, 0.75, 
                                             0.9), names = FALSE, ...)
          delta <- min(qy[5] - qy[4], qy[2] - qy[1])
          a <- qy[1] - 3 * delta
          b <- qy[5] + 3 * delta
          yab <- y[y >= a & y <= b]
          k <- kmax
          ok2 <- FALSE
          while (!ok2) {
            b <- seq(from = min(yab), to = max(yab), length = k + 
                       1)
            n <- c(tabulate(findInterval(yab, b[-(k + 1)])), 
                   0)
            N <- sum(n)
            v <- as.numeric(n >= N/k)
            w <- which.max(v)
            v2 <- v[w:(k + 1)]
            w2 <- which.min(v2) + w - 1
            v3 <- v[w2:(k + 1)]
            nc <- sum(n[w:(w2 - 1)])
            if (any(v3 == 1) && nc == 1) {
              if (k > 3) {
                k <- k - 1
              }
              else if (k == 3) {
                if (n[3] > 1) {
                  w <- 3
                  w2 <- 4
                }
                ok2 <- TRUE
              }
              else {
                stop("k < 3", call. = FALSE)
              }
            }
            else if (any(v3 == 1) && nc > 1) {
              if (k > 3) {
                k <- k - 1
              }
              else if (k == 3) {
                if (n[3] > 1) {
                  p1 <- (1/n[1]) * prod(diff(yab[yab >= b[w] & 
                                                   yab <= b[w2]]))
                  p2 <- (1/n[3]) * prod(diff(yab[yab >= b[3] & 
                                                   yab <= b[4]]))
                  if (p1 > p2) {
                    w <- 3
                    w2 <- 4
                  }
                }
                ok2 <- TRUE
              }
              else {
                stop("k < 3", call. = FALSE)
              }
            }
            else if (!any(v3 == 1)) {
              ok2 <- TRUE
            }
          }
          nc <- sum(n[w:(w2 - 1)])
          d <- abs((qy[4] + qy[2] - 2 * qy[3])/(qy[4] - qy[2]))
          nc2 <- ny * (1 - d)
          y <- yab[yab >= b[w] & yab <= b[w2]]
          if (nc == ny) {
            ok1 <- TRUE
          }
          else if (nc <= ifelse(nx < 30, nx/3, bw * nc2)) {
            ok1 <- TRUE
          }
          else {
            ok1 <- FALSE
          }
        }
        stats::median(y)
      }
      
      pml <- apply(p, MARGIN=2, FUN = asselin)
      names(pml) <- colnames(p)
      
    } else {
      
      
      L <- x$resultLnL.total[[chain]]
      if (is.null(L)) L <- x$resultLnL[[chain]]
      
      pos <- NULL
      
      if (is.numeric(index)) {
        pos <- index
      } else {
        
        if (index=="best") {
          pos <- which.max(L)
          message(paste("The best likelihood has been observed at iteration", pos))
        }
        
        if (index=="last") {
          pos <- length(L)
        }
      }
      
      if (is.null(pos)) {
        stop("index is not recognized")
      }
      
      pml <- as.numeric(p[pos,])
      names(pml) <- names(p[pos,])
    }
  }
}
return(pml)

}
