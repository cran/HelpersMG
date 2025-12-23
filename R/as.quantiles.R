#' @title Extract quantile distribution from mcmcComposite object
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A data.frame with quantiles
#' @param x A mcmcComposite obtained as a result of \code{MHalgoGen()} function
#' @param chain The number of the chain in which to get parameters or "all"
#' @param fun The function to apply the parameters
#' @param probs The probability to get quantiles
#' @param xlim The values to apply in fun
#' @param nameparxlim The name of the parameter for xlim
#' @param namepar The name of parameters from mcmc object to be used in fun
#' @description Extract quantile distribution from mcmcComposite object
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
#' Prior1=c(10, 0.5), Prior2=c(2, 0.5), SDProp=c(1, 1), 
#' Min=c(-3, 0), Max=c(100, 10), Init=c(10, 2), stringsAsFactors = FALSE, 
#' row.names=c('mean', 'sd'))
#' mcmc_run <- MHalgoGen(n.iter=10000, parameters=parameters_mcmc, data=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=1)
#' k <- as.quantiles(x=mcmc_run, namepar="mean")
#' k <- as.quantiles(x=mcmc_run, namepar="mean", 
#'                  xlim=c(1:5), nameparxlim="sd", 
#'                  fun=function(...) return(sum(as.numeric(list(...)))))
#' }
#' @export


as.quantiles <- function(x, chain="all", fun=function(...) return(as.numeric(list(...))), 
                            probs=c(0.025, 0.975), xlim=NULL, 
                            nameparxlim=NULL, namepar=NULL) {
  
  # chain=1; fun=function(x) {x}; probs=c(0.025, 0.975); xlim=NULL; nameparxlim=NULL; namepar=NULL
  if (chain[1] == "all") chain <- seq_along(x$resultMCMC)
  
  # if (total) p <- x$resultMCMC.total[[chain[1]]] else p <- x$resultMCMC[[chain[1]]]
  
  p <- NULL
  
  if (length(chain) >= 1) {
    for (i in chain) {
      p <- rbind(p, x$resultMCMC[[i]])
    }
  }
  
  if (is.null(p)) stop("No data are available.")
  
	# p <- x$resultMCMC[[chain]]
	df <- matrix(NA, ncol=max(1, length(xlim)))
	df <- df[-1, ]
	
	for (l in 1:nrow(p)) {
	  if (length(xlim) == 0) {
	    df <- rbind(df, rep(NA, 1))
	    par <- p[l, ]
	    if (!is.null(namepar)) par <- par[namepar]
	    parx <- as.list(par)
	    h <- do.call(fun, args=parx)
	    df[l, 1] <- h
	  } else {
	    par <- p[l, ]
	    if (!is.null(namepar)) par <- par[namepar]
	    df <- rbind(df, rep(NA, length(xlim)))
	    for (i in 1:length(xlim)) {
	      m <- xlim[i]
	      names(m) <- nameparxlim
	      parx <- as.list(c(m, par))
	      h <- do.call(fun, args=parx)
	      df[l, i] <- h
	    }
	  }
	}

	q <- apply(df, MARGIN = 2, FUN = function(x) quantile(x, probs))
	colnames(q) <- as.character(xlim)
	return(q)

}
