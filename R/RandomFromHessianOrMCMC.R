#' RandomFromHessianOrMCMC returns random numbers based on Hessian matrix or MCMC
#' @title Random numbers based on Hessian matrix or MCMC
#' @author Marc Girondot
#' @return Returns a list with two data.frames named df_random and df_fn
#' @param Hessian An Hessian matrix
#' @param mcmc A result from MHalgogen()
#' @param chain MCMC chain to be used
#' @param fitted.parameters The fitted parameters
#' @param fixed.parameters The fixed parameters
#' @param replicates Number of replicates to generate
#' @param fn The function to apply to each replicate
#' @param probs Probability for quantiles
#' @param silent Should the function display some information
#' @param ... Parameters send to fn function
#' @description A data.frame with one column for each parameter
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' val <- rnorm(100, mean=20, sd=5)+(1:100)/10
#' # Return -ln L of values in val in Gaussian distribution with mean and sd in par
#' fitnorm <- function(par, data) {
#'   -sum(dnorm(data, par["mean"], abs(par["sd"]), log = TRUE))
#' }
#' # Initial values for search
#' p<-c(mean=20, sd=5)
#' # fit the model
#' result <- optim(par=p, fn=fitnorm, data=val, method="BFGS", hessian=TRUE)
#' # Using Hessian
#' df <- RandomFromHessianOrMCMC(Hessian=result$hessian, fitted.parameters=result$par)$df_random
#' hist(df[, 1], main="mean")
#' hist(df[, 2], main="sd")
#' plot(df[, 1], df[, 2], xlab="mean", ylab="sd", las=1, bty="n")
#' 
#' # Using MCMC
#' parameters_mcmc <- data.frame(Density=c('dnorm', 'dlnorm'), 
#' Prior1=c(10, 0.5), Prior2=c(2, 0.5), SDProp=c(0.35, 0.2), 
#' Min=c(-3, 0), Max=c(100, 10), Init=c(10, 2), stringsAsFactors = FALSE, 
#' row.names=c('mean', 'sd'))
#' # Use of trace and traceML parameters
#' # trace=1 : Only one likelihood is printed
#' mcmc_run <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, data=val, 
#' parameters_name = "par", 
#' likelihood=fitnorm, n.chains=1, n.adapt=100, thin=1, trace=1)
#' df <- RandomFromHessianOrMCMC(mcmc=mcmc_run, fitted.parameters=NULL)$df_random
#' hist(df[, 1], main="mean")
#' hist(df[, 2], main="sd")
#' plot(df[, 1], df[, 2], xlab="mean", ylab="sd", las=1, bty="n")
#' 
#' # Using a function
#' fitnorm <- function(par, data, x) { 
#'   y=par["a"]*(x)+par["b"]
#'   -sum(dnorm(data, y, abs(par["sd"]), log = TRUE))
#' }
#' p<-c(a=0.1, b=20, sd=5)
#' # fit the model
#' x <- 1:100
#' result <- optim(par=p, fn=fitnorm, data=val, x=x, method="BFGS", hessian=TRUE)
#' # Using Hessian
#' df <- RandomFromHessianOrMCMC(Hessian=result$hessian, fitted.parameters=result$par, 
#'                               fn=function(par) (par["a"]*(x)+par["b"]))
#' plot(1:100, val)
#' lines(1:100, df$quantile["50%", ])
#' lines(1:100, df$quantile["2.5%", ], lty=2)
#' lines(1:100, df$quantile["97.5%", ], lty=2)
#' }
#' @export

RandomFromHessianOrMCMC <- function(Hessian=NULL, mcmc=NULL, chain=1, 
                                    fitted.parameters=NULL, 
                                    fixed.parameters=NULL, 
                                    probs=c(0.025, 0.5, 0.975), 
                                    replicates=10000, fn=NULL, silent=FALSE, ...) {
  
  if (is.null(Hessian) & (is.null(mcmc))) stop("Both Hessian and mcmc cannot be NULL")
  if (!is.null(Hessian) & (!is.null(mcmc))) stop("Both Hessian and mcmc cannot be provided")
  
  if (!is.null(Hessian) & (is.null(fitted.parameters))) stop("Fitted.parameters cannot be NULL with Hessian")
  
  if (!is.null(Hessian)) {
    sigma <- try(solve(Hessian), silent = TRUE)
    if (all(class(sigma) != "try-error")) {
      if (!silent) cat("Estimation using variance-covariance matrix")
      s. <- svd(sigma)
      R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
      df_random <- matrix(rnorm(replicates * ncol(sigma)), nrow = replicates, byrow = TRUE) %*% R
      df_random <- sweep(df_random, 2, fitted.parameters[rownames(Hessian)], "+")
      colnames(df_random) <- rownames(Hessian)
    } else {
      # J'ai une erreur sur l'inversion de la matrice hessienne; 
      # Je prends les SE
      if (!silent) cat("Estimation using variances")
      se <- SEfromHessian(Hessian)
      
      df_random <- matrix(data = NA, ncol=length(se), nrow=replicates)
      colnames(df_random) <- names(se)
      
      for (i in names(se)) {
        df_random[, i] <- rnorm(replicates, fitted.parameters[i], se[i])
      }
    }
  }
  
  if (!is.null(mcmc)) {
    if (replicates<nrow(mcmc$resultMCMC[[chain]])) {
      df_random <- mcmc$resultMCMC[[chain]][sample(x=1:nrow(mcmc$resultMCMC[[chain]]), size=replicates), ]
    } else {
      df_random <- mcmc$resultMCMC[[chain]]
    }
  }
  
  if (!is.null(fixed.parameters)) {
    ajouter <- matrix(rep(fixed.parameters, replicates), 
                      nrow=replicates, byrow=TRUE,
                      dimnames = list(c(NULL), names(fixed.parameters)))
    df_random <- cbind(df_random, ajouter)
  }
  
  if (!is.null(fn)) {
    parx <- list(...)
    # print(str(parx))
    df_fn <- apply(df_random, MARGIN=1, 
                   FUN=function(par) do.call(fn, modifyList(list(par=par), parx)))
    if (is.null(dim(df_fn))) {
      df_fn <- matrix(df_fn, ncol=1)
    } else {
      df_fn <- t(df_fn)
    }
    
    q <- apply(df_fn, MARGIN=2, FUN=function(x) quantile(x, probs = probs))
  } else {
    df_fn <- NULL
    q <- NULL
  }
  
  
  
  
  return(list(random=df_random, fn=df_fn, quantile=q))
}
