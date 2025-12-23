#' RandomFromHessianOrMCMC returns random numbers based on Hessian matrix or MCMC
#' @title Random numbers based on Hessian matrix or MCMC
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return Returns a list with three data.frames named random, fn, and quantiles
#' @param Hessian An Hessian matrix
#' @param se A named vector with SE of parameters
#' @param mcmc A result from MHalgogen()
#' @param chain MCMC chain to be used or "all"
#' @param regularThin If TRUE, use regular thin for MCMC or use a number
#' @param MinMax A data.frame with at least two columns: Min and Max and rownames being the variable names
#' @param fitted.parameters The fitted parameters
#' @param fixed.parameters The fixed parameters
#' @param replicates Number of replicates to generate the randoms
#' @param method Can be NULL, "SE", "Hessian", "MCMC", or "PseudoHessianFromMCMC"
#' @param fn The function to apply to each replicate
#' @param probs Probability for quantiles
#' @param silent Should the function display some information
#' @param ParTofn Name of the parameter to send random values to fn
#' @param ... Parameters send to fn function
#' @description If it is very long, use silent parameter to check if something goes wrong.\cr
#' If replicates is NULL or is 0, or if method is NULL, parameters are just copied into data.frame.\cr
#' If method is NULL, replicate.CI is set to 0.\cr
#' If method is hessian, it will generate replicate.CI random numbers using 
#' Hessian matrix with covariance.\cr
#' If method is se, it will generate replicate.CI random numbers using 
#' SE values then without covariance.\cr
#' If method is mcmc, it will generate replicate.CI random numbers using 
#' random samples of the MCMC or the regularThin number if it is a number.\cr
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
#' df <- RandomFromHessianOrMCMC(Hessian=result$hessian, 
#'                               fitted.parameters=result$par, 
#'                               method="Hessian")$random
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
#' df <- RandomFromHessianOrMCMC(mcmc=mcmc_run, fitted.parameters=NULL, 
#'                               method="MCMC")$random
#' hist(df[, 1], main="mean")
#' hist(df[, 2], main="sd")
#' plot(df[, 1], df[, 2], xlab="mean", ylab="sd", las=1, bty="n")
#' 
#' # Return the two first elements of the MCMC
#' df <- RandomFromHessianOrMCMC(mcmc=mcmc_run, fitted.parameters=NULL, 
#'                method="MCMC", replicates = 2, regularThin = c(1, 2))$random
#' 
#' # Using a function fn
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
#'                               method="Hessian", 
#'                               fn=function(par) (par["a"]*(x)+par["b"]))
#' plot(1:100, val)
#' lines(1:100, df$quantiles["50%", ])
#' lines(1:100, df$quantiles["2.5%", ], lty=2)
#' lines(1:100, df$quantiles["97.5%", ], lty=2)
#' }
#' @export

RandomFromHessianOrMCMC <- function(se=NULL                          , 
                                    Hessian=NULL                     , 
                                    mcmc=NULL                        , 
                                    chain="all"                      , 
                                    regularThin=TRUE                 , 
                                    MinMax=NULL                      , 
                                    fitted.parameters=NULL           , 
                                    fixed.parameters=NULL            , 
                                    method=NULL                      , 
                                    probs=c(0.025, 0.5, 0.975)       , 
                                    replicates=10000                 , 
                                    fn=NULL                          , 
                                    silent=FALSE                      , 
                                    ParTofn="par"                    , 
                                    ...) {
  
  # se=NULL; Hessian=NULL; mcmc=NULL; chain=1
  # method=NULL
  # regularThin=TRUE
  # MinMax=NULL
  # fitted.parameters=NULL
  # fixed.parameters=NULL
  # probs=c(0.025, 0.5, 0.975)
  # replicates=10000; fn=NULL; silent=FALSE; ParTofn="par"
  
  if (is.null(replicates)) replicates <- 0
  if (is.null(method)) method <- "null"
  method <- tolower(method)
  
  method <- match.arg(method, choices=c("null", "se", "hessian", "mcmc", "pseudohessianfrommcmc"))
  
  df_random <- NULL
  
  if (method == "null" | replicates == 0)  {
    if (!is.null(fitted.parameters)) {
      # Je retourne un dataframe avec simplement des recopies
      df_random <- as.data.frame(matrix(data = rep(fitted.parameters, ifelse(replicates==0, 1, replicates)), 
                                        nrow = ifelse(replicates==0, 1, replicates), 
                                        byrow = TRUE))
      colnames(df_random) <- names(fitted.parameters)
    }
  }
  
  if (!is.null(Hessian) & (method == "hessian") & (is.null(fitted.parameters))) stop("Fitted.parameters or Hessian cannot be NULL with Hessian method.")
  if (!is.null(se) & (method == "se") & (is.null(fitted.parameters))) stop("Fitted.parameters or se cannot be NULL with se method.")
  
  if (method == "mcmc") {
    
    if (chain[1] == "all") chain <- seq_along(mcmc$resultMCMC)
    df_random <- NULL
    if (length(chain) >= 1) {
      for (i in chain) {
        df_random <- rbind(df_random, mcmc$resultMCMC[[i]])
      }
    } else {
      stop("No data have been selected.")
    }
    
    if (isTRUE(regularThin)) {
      if (replicates <= nrow(df_random)) {
        df_random <- df_random[seq(from=1, to=nrow(df_random), length.out=replicates), , drop=FALSE]
      } else {
        stop("When regularThin is TRUE, replicates must be lower or equal to number of MCMC iterations.")
      }
    } else {
      if (isFALSE(regularThin)) {
      if (replicates < nrow(df_random)) {
        df_random <- df_random[sample(x=1:nrow(df_random), size=replicates), , drop=FALSE]
      } else {
        df_random <- df_random[sample(x=1:nrow(df_random), size=replicates, replace = TRUE), , drop=FALSE]
      }
      } else {
        if (is.numeric(regularThin)) {
          regularThin <- rep(regularThin, replicates)[1:replicates]
          df_random <- df_random[regularThin, , drop=FALSE]
        } else {
          stop("I don't understand the format of regularThin.")
        }
      }
    }
  }
  
  if (method == "pseudohessianfrommcmc") {
    Hessian <- solve(cov(mcmc$resultMCMC[[chain]]))
    fitted.parameters <- as.parameters(mcmc)
    method <- "hessian"
  }
  
  if (method == "hessian") {
    
    sigma  <- try(solve(Hessian), silent=TRUE)
    # Add all: 22/4/2020; any !
    if (inherits(sigma, "try-error")) {
      if (!silent) warning("Error in Hessian matrix inversion")
      Hessianx <- try(as.matrix(getFromNamespace("nearPD", ns="Matrix")(Hessian)$mat), silent=TRUE)
      # 29/1/2021
      if (inherits(Hessianx, "try-error")) {
        if (!silent) warning("Error in estimation of the Nearest Positive Definite Matrix. Calculates the Moore-Penrose generalized inverse. Use result with caution.")
        sigma  <- try(ginv(Hessian), silent=TRUE)
      } else {
        if (!silent) warning("Calculates the Nearest Positive Definite Matrix. Use result with caution.")
        sigma  <- try(solve(Hessianx), silent=TRUE)
      }
    } 
    
    if (!inherits(sigma, "try-error")) {
      if (!silent) cat("Estimation using variance-covariance matrix\n")
      df_random <- matrix(data=NA, nrow=1, ncol=nrow(Hessian))
      df_random <- as.data.frame(df_random[-1, ])
      colnames(df_random) <- rownames(Hessian)
      
      reste.replicates <- replicates - nrow(df_random)
      repeat {
        if (!silent) cat(paste0("Still ", as.character(reste.replicates), " replicates\n"))
        s. <- svd(sigma)
        R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
        df_random_int <- matrix(rnorm(reste.replicates * ncol(sigma)), 
                                nrow = reste.replicates, byrow = TRUE) %*% R
        df_random_int <- sweep(df_random_int, 2, fitted.parameters[rownames(Hessian)], "+")
        colnames(df_random_int) <- rownames(Hessian)
        if (!is.null(MinMax)) {
          L <- rep(TRUE, nrow(df_random_int))
          for (i in colnames(df_random_int)) {
            if ((!is.na(MinMax[i, "Min"])) & (!is.na(MinMax[i, "Max"]))) {
              # Je veux TRUE si c'est bon
              L <- L & (MinMax[i, "Min"] <= df_random_int[, i]) & 
                (MinMax[i, "Max"] >= df_random_int[, i])
            }
          }
          df_random_int <- df_random_int[L, , drop=FALSE]
        }
        df_random <- rbind(df_random, df_random_int)
        reste.replicates <- replicates - nrow(df_random)
        if (reste.replicates == 0) break
      }
      
    } else {
      if (!silent) cat("Error using Hessian; Estimation using standard errors estimated from Hessian")
      se <- SEfromHessian(Hessian)
      method <- "se"
    }
  }
  
  if ((method == "se") & (is.null(se)) & (!is.null(Hessian))) {
    if (!silent) cat("Estimation using standard errors estimated from Hessian")
    se <- SEfromHessian(Hessian)
  }
  
  if (method == "se") {
    # J'ai une erreur sur l'inversion de la matrice hessienne; ou bien j'ai des se
    # Je prends les SE
    df_random <- matrix(data=NA, nrow=1, ncol=length(se))
    df_random <- as.data.frame(df_random[-1, ])
    colnames(df_random) <- names(se)
    reste.replicates <- replicates - nrow(df_random)
    repeat {
      df_random_int <- matrix(data = NA, ncol=length(se), nrow=reste.replicates)
      colnames(df_random_int) <- names(se)
      for (i in colnames(df_random_int)) {
        df_random_int[, i] <- rnorm(reste.replicates, fitted.parameters[i], se[i])
      }
      if (!is.null(MinMax)) {
        L <- rep(TRUE, nrow(df_random_int))
        for (i in colnames(df_random_int)) {
          # Je veux TRUE si c'est bon
          L <- L & (MinMax[i, "Min"] <= df_random_int[, i]) & 
            (MinMax[i, "Max"] >= df_random_int[, i])
        }
        df_random_int <- df_random_int[L, ]
      }
      df_random <- rbind(df_random, df_random_int)
      reste.replicates <- replicates - nrow(df_random)
      if (reste.replicates == 0) break
    } 
    
  }
  
  
  if (!is.null(fixed.parameters)) {
    ajouter <- matrix(rep(fixed.parameters, ifelse(replicates==0, 1, replicates)), 
                      nrow=ifelse(replicates==0, 1, replicates), byrow=TRUE,
                      dimnames = list(c(NULL), names(fixed.parameters)))
    df_random <- cbind(df_random, ajouter)
  }
  
  if (is.null(df_random)) {
    stop("Impossible to generate random numbers")
  }
  
  if (!is.null(fn)) {
    parx <- list(...)
    # print(str(parx))
    df_fn <- apply(df_random, MARGIN=1, 
                   FUN=function(par_dfr) {
                     lst <- list(par=unlist(par_dfr))
                     names(lst) <- ParTofn
                     do.call(fn, modifyList(lst, parx))
                   })
    # df_fn peut être un vecteur ou une liste
    if (is.list(df_fn)) {
      df_fn <- unlist(df_fn)
      df_fn <- matrix(df_fn, nrow=ifelse(replicates==0, 1, replicates), byrow = TRUE)
    } else {
      
      if (is.null(dim(df_fn))) {
        df_fn <- matrix(df_fn, ncol=1)
      } else {
        df_fn <- t(df_fn)
      }
    }
    
    q <- apply(df_fn, MARGIN=2, FUN=function(x) quantile(x, probs = probs, na.rm=TRUE))
    if (any(is.na(df_fn))) warning("Some NA are produced; take care because outputs are biased")
  } else {
    df_fn <- NULL
    q <- NULL
  }
  
  
  
  
  return(list(random=df_random, fn=df_fn, quantiles=q))
}
