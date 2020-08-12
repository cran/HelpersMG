#' MHalgoGen is a function to use mcmc with Metropolis-Hastings algorithm
#' @title Monte-Carlo Markov-chain with Metropolis-Hastings algorithm
#' @author Marc Girondot
#' @return A mcmcComposite object with all characteristics of the model and mcmc run
#' @param n.iter Number of iterations for each chain
#' @param parameters A data.frame with priors; see description and examples
#' @param likelihood The function that returns -ln likelihood using data and parameters
#' @param n.chains Number of chains
#' @param n.adapt Number of iteration to stabilize likelihood
#' @param thin Interval for thinning likelihoods
#' @param trace Or FALSE or period to show progress
#' @param traceML TRUE or FALSE to show ML
#' @param intermediate Or NULL of period to save intermediate result
#' @param filename Name of file in which intermediate results are saved
#' @param adaptive Should an adaptive process for SDProp be used
#' @param adaptive.lag  Lag to analyze the SDProp value in an adaptive context
#' @param adaptive.fun Function used to change the SDProp
#' @param previous The content of the file in which intermediate results are saved
#' @param parameters_name The name of the parameters in the likelihood function, default is "x"
#' @param ... Parameters to be transmitted to likelihood function
#' @description The parameters must be stored in a data.frame with named rows for each parameter with the following columns:\cr
#' \itemize{
#'   \item Density. The density function name, example \code{dnorm}, \code{dlnorm}, \code{dunif}, \code{dbeta}
#'   \item Prior1. The first parameter to send to the \code{Density} function
#'   \item Prior2. The second parameter to send to the \code{Density} function
#'   \item SDProp. The standard error from new proposition value of this parameter
#'   \item Min. The minimum value for this parameter
#'   \item Max. The maximum value for this parameter
#'   \item Init. The initial value for this parameter
#' }
#' This script has been deeply modified from a MCMC script provided by Olivier Martin (INRA, Paris-Grignon).\cr
#' The likelihood function must use a parameter named parameters_name for the nammed parameters.\cr
#' For adaptive mcmc, see:\cr
#' Rosenthal, J. S. 2011. Optimal Proposal Distributions and Adaptive MCMC. Pages 93-112 in S. Brooks, A. Gelman, 
#' G. Jones, and X.-L. Meng, editors. MCMC Handbook. Chapman and Hall/CRC.
#' @family mcmcComposite functions
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' require(coda)
#' val <- rnorm(30, 10, 2)
#' dnormx <- function(data, x) {
#'  data <- unlist(data)
#'  return(-sum(dnorm(data, mean=x['mean'], sd=x['sd'], log=TRUE)))
#' }
#' parameters_mcmc <- data.frame(Density=c('dnorm', 'dlnorm'), 
#' Prior1=c(10, 0.5), Prior2=c(2, 0.5), SDProp=c(0.35, 0.2), 
#' Min=c(-3, 0), Max=c(100, 10), Init=c(10, 2), stringsAsFactors = FALSE, 
#' row.names=c('mean', 'sd'))
#' # Use of trace and traceML parameters
#' # trace=1 : Only one likelihood is printed
#' mcmc_run <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, data=val, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=1)
#' # trace=10 : 10 likelihoods are printed
#' mcmc_run <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, data=val, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=10)
#' # trace=TRUE : all likelihoods are printed
#' mcmc_run <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, data=val, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=TRUE)
#' # trace=FALSE : No likelihood is printed
#' mcmc_run <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, data=val, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=FALSE)
#' # traceML=TRUE : values when likelihood is better are shown
#' mcmc_run <- MHalgoGen(n.iter=100, parameters=parameters_mcmc, data=val, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=TRUE, traceML=TRUE)
#' mcmc_run <- MHalgoGen(n.iter=100, parameters=parameters_mcmc, data=val, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=FALSE, traceML=TRUE)
#' 
#' plot(mcmc_run, xlim=c(0, 20))
#' plot(mcmc_run, xlim=c(0, 10), parameters="sd")
#' library(graphics)
#' library(fields)
#' # show a scatter plot of the result
#' x <- mcmc_run$resultMCMC[[1]][, 1]
#' y <- mcmc_run$resultMCMC[[1]][, 2]
#' marpre <- par(mar=c(4, 4, 2, 6)+0.4)
#' smoothScatter(x, y)
#' # show a scale
#' n <- matrix(0, ncol=128, nrow=128)
#' xrange <- range(x)
#' yrange <- range(y)
#' for (i in 1:length(x)) {
#'   posx <- 1+floor(127*(x[i]-xrange[1])/(xrange[2]-xrange[1]))
#'   posy <- 1+floor(127*(y[i]-yrange[1])/(yrange[2]-yrange[1]))
#'   n[posx, posy] <- n[posx, posy]+1
#' }
#' image.plot(legend.only=TRUE, zlim= c(0, max(n)), nlevel=128, 
#'  col=colorRampPalette(c("white", blues9))(128))
#' # Compare with a heatmap
#' x <- seq(from=8, to=12, by=0.2)
#' y <- seq(from=1, to=4, by=0.2)
#' df <- expand.grid(mean=x, sd=y)
#' df <- cbind(df, L=rep(0, length(nrow(df))))
#' for (i in 1:nrow(df)) df[i, "L"] <- -sum(dnorm(val, df[i, 1], df[i, 2], log = TRUE))
#' hm <- matrix(df[, "L"], nrow=length(x))
#' par(mar = marpre)
#' image.plot(x=x, y=y, z=hm, las=1)
#' # Diagnostic function from coda library
#' mcmcforcoda <- as.mcmc(mcmc_run)
#' #' heidel.diag(mcmcforcoda)
#' raftery.diag(mcmcforcoda)
#' autocorr.diag(mcmcforcoda)
#' acf(mcmcforcoda[[1]][,"mean"], lag.max=20, bty="n", las=1)
#' acf(mcmcforcoda[[1]][,"sd"], lag.max=20, bty="n", las=1)
#' batchSE(mcmcforcoda, batchSize=100)
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmcforcoda)$statistics[,"Time-series SE"]
#' summary(mcmc_run)
#' as.parameters(mcmc_run)
#' lastp <- as.parameters(mcmc_run, index="last")
#' parameters_mcmc[,"Init"] <- lastp
#' # The n.adapt set to 1 is used to not record the first set of parameters
#' # then it is not duplicated (as it is also the last one for 
#' # the object mcmc_run)
#' mcmc_run2 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, x=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=1, thin=1, trace=1)
#' mcmc_run3 <- merge(mcmc_run, mcmc_run2)
#' ####### no adaptation, n.adapt must be 0
#' parameters_mcmc[,"Init"] <- c(mean(x), sd(x))
#' mcmc_run3 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, x=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=0, thin=1, trace=1)
#' # Here is how to use adaptive mcmc
#' mcmc_run <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, data=val, adaptive = FALSE, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=1)
#' 1-rejectionRate(as.mcmc(mcmc_run))
#' mcmc_run <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, data=val, adaptive = TRUE,  
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=1)
#' 1-rejectionRate(as.mcmc(mcmc_run))
#' # To see the dynamics :
#' var <- "mean"
#' par(mar=c(4, 4, 1, 1)+0.4)
#' plot(1:nrow(mcmc_run$resultMCMC[[1]]), mcmc_run$resultMCMC[[1]][, var], type="l", 
#'        xlab="Iterations", ylab=var, bty="n", las=1)
#' }
#' @export

# Algo Metropolis-Hastings
# ------------------------

MHalgoGen<-function(likelihood=stop("A likelihood function must be supplied"), 
                    parameters_name="x",
                    parameters=stop("Priors  must be supplied"), ..., 
                    n.iter=10000, n.chains = 1, n.adapt = 100, thin=30, 
                    trace=FALSE, traceML=FALSE, 
                    adaptive = FALSE, adaptive.lag = 500, 
                    adaptive.fun = function(x) {ifelse(x>0.234, 1.3, 0.7)},
                    intermediate=NULL, filename="intermediate.Rdata",
                    previous=NULL)
  
{
  
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("coda package is necessary for this function")
  }
  
  previousML <- -Inf
  
  # likelihood=NULL; parameters_name="x"; parameters=NULL; n.iter=10000; n.chains = 1; n.adapt = 100; thin=30; trace=FALSE; traceML=FALSE; intermediate=NULL; filename="intermediate.Rdata"; previous=NULL; adaptive = FALSE; adaptive.lag = 500; adaptive.fun = function(x) {ifelse(x>0.234, 1.3, 0.7)}
  # datax <- list(temperatures=result$data, derivate=result$derivate, test=result$test, M0=result$M0, fixed.parameters=result$fixed.parameters, weight=result$weight, out="Likelihood",  progress=FALSE, warnings=FALSE, likelihood=getFromNamespace("info.nests", ns = "embryogrowth"))
  
  if (is.null(previous)) {
    nbvar <- dim(parameters)[1]
    deb_kk <- 1
    deb_i <- 2
    res <- as.list(NULL)
    deb_varp2 <- matrix(rep(NA, (nbvar+1)*(n.adapt+n.iter+2)), ncol=nbvar+1)
    deb_varp<-matrix(rep(NA, (nbvar+1)*(n.adapt+n.iter+2)), ncol=nbvar+1)
    colnames(deb_varp)<-c(rownames(parameters), "Ln L")
    colnames(deb_varp2)<-c(rownames(parameters), "Ln L")
    t <- as.character(trace)
    cpt_trace <- 0
    res<-as.list(NULL)
    resL<-as.list(NULL)
    # datax <- list(temperatures=result$data, derivate=result$derivate, test=result$test, M0=result$M0, fixed.parameters=result$fixed.parameters, weight=result$weight, out="Likelihood",  progress=FALSE, warnings=FALSE, likelihood=getFromNamespace("info.nests", ns = "embryogrowth"))
    # datax <- list(data)
    datax <- list(...)
  } else {
    n.iter <- previous$n.iter
    parameters <- previous$parameters
    # Initialisation
    nbvar<-dim(parameters)[1]
    datax <- previous$datax
    likelihood <- previous$likelihood
    parameters_name <- previous$parameters_name
    n.chains <- previous$n.chains
    n.adapt <- previous$n.adapt
    thin <- previous$thin
    trace <- previous$trace
    traceML <- previous$traceML
    intermediate <- previous$intermediate
    filename <- previous$filename
    deb_kk <- previous$chain
    deb_i <- previous$iter
    res <- previous$res
    deb_varp2 <- previous$varp2
    deb_varp <- previous$varp
    cpt <- previous$cpt
    cpt_trace <- previous$cpt_trace
    resL <- previous$resL
    sdg <- previous$sdg
    MaxL <- previous$MaxL
    t <- as.character(previous$trace)
    adaptive <- previous$adaptive
    adaptive.lag <- previous$adaptive.lag
    adaptive.fun <- previous$adaptive.fun
  }
  
  pt <- NULL
  if (t=="TRUE") {pt <- 1;tf <- TRUE}
  if (t=="FALSE") {pt <- 0;tf <- FALSE}
  if (is.null(pt)) {
    tf <- TRUE
    pt <- floor((n.adapt+n.iter)/trace)
  }
  
  
  for (kk in deb_kk:n.chains) {
    
    varp <- deb_varp
    varp2 <- deb_varp2
    
    if (is.null(previous)) {
      
      deb_varp2 <- matrix(rep(NA, (nbvar+1)*(n.adapt+n.iter+2)), ncol=nbvar+1)
      deb_varp <- matrix(rep(NA, (nbvar+1)*(n.adapt+n.iter+2)), ncol=nbvar+1)
      colnames(deb_varp)<-c(rownames(parameters), "Ln L")
      colnames(deb_varp2)<-c(rownames(parameters), "Ln L")
      
      
      varp[1, 1:nbvar] <- as.numeric(parameters[1:nbvar, 'Init'])
      
      param <- list(varp[1, 1:nbvar])
      names(param) <- parameters_name
      
      # ll <- c(datax, param)
      # names(ll) <- c("data", parameters_name)
      
      varp[1, "Ln L"] <- -do.call(likelihood, modifyList(datax, param))
      # likelihood(par=param$par, fixed.parameters = datax$fixed.parameters, res_1=datax$res_1, res_2=datax$res_2)
      cpt <- 1
      varp2[cpt, 1:nbvar] <- varp[1, 1:nbvar]
      varp2[cpt, "Ln L"] <- varp[1, "Ln L"]
      cpt <- 2
      
      
      if (trace) {
        cat(paste("Chain ", kk, ": [", 1, "] ",as.numeric(varp[1, nbvar+1]), "\n", sep=""))
      # } else {
      #   cat(paste("Chain ", kk, "\n", sep=""))
      }
      
      if (traceML) {
        if (as.numeric(varp[1, "Ln L"]) > previousML) {
          previousML <- varp[1, "Ln L"]
          cat(d(varp[1, 1:nbvar]))
          cat("\n")
        }
      }
      MaxL <- varp[1, ]
      
      # sdg=NULL
      
      
      # 18/1/2013
      # for(i in 1:nbvar) sdg<-c(sdg, as.numeric(parameters[i, 'SDProp']))
      # sdg <- c(sdg, as.numeric(parameters[1:nbvar, 'SDProp']))
      sdg <- as.numeric(parameters[1:nbvar, 'SDProp'])
      
      
      # previous <- NULL
    }
    
    
    if (is.data.frame(parameters)) {
      Prior <- as.matrix(parameters[,c("Prior1", "Prior2")])
      Limites <- as.matrix(parameters[,c("Min", "Max")])
    } else {
      Prior<-matrix(as.numeric(parameters[,c("Prior1", "Prior2")]), ncol=2)
      Limites<-matrix(as.numeric(parameters[,c("Min", "Max")]), ncol=2)
    }
    
    dfun <- as.character(parameters[,"Density"])
    
    # Itérations
    for (i in deb_i:(n.adapt+n.iter)) {
      
      # est-ce que je sauve où j'en suis
      if (!is.null(intermediate))
        if (i %% intermediate==0) {
          itr <- list(chain=kk, iter=i, varp2=varp2, res=res, 
                      n.iter=n.iter,
                      parameters=parameters,
                      datax=datax,
                      likelihood=likelihood,
                      parameters_name=parameters_name,
                      n.chains=n.chains,
                      n.adapt=n.adapt,
                      thin=thin,
                      trace=trace,
                      traceML=traceML,
                      intermediate=intermediate,
                      filename=filename, 
                      varp=varp, 
                      cpt=cpt,
                      cpt_trace=cpt_trace, 
                      resL=resL, 
                      sdg=sdg, MaxL=MaxL, 
                      adaptive=adaptive,
                      adaptive.lag=adaptive.lag,
                      adaptive.fun=adaptive.fun)
          save(itr, file=filename)
        }
      
      
      ###### Pas clair si previous
      deb_i <- 2
      Lprevious <- LpreviousT <- varp[i-1, "Ln L"]
      newvarp <- varp[i-1, 1:nbvar]
      
      for (j in 1:nbvar) {	
        # Nouveau paramètre
        propvarp <- newvarp
        propvarp[j] <- propvarp[j]+rnorm(1, mean=0, sd=sdg[j])
        
        if (propvarp[j]<=Limites[j,2] && propvarp[j]>=Limites[j,1]) {
          param <- list(propvarp)
          names(param) <- parameters_name
          Lprevious2 <- -do.call(likelihood, modifyList(datax, param))
          logratio <- get(dfun[j])(propvarp[j],Prior[j,1],Prior[j,2],log=TRUE) + Lprevious2 -
            (get(dfun[j])(newvarp[j],Prior[j,1],Prior[j,2],log=TRUE)+LpreviousT)
          alpha<-min(c(1,exp(logratio)))
          # 15/2/2015 Pour éviter des erreurs
          if (!is.finite(alpha)) alpha <- -1
          if (runif(1, min=0, max=1)<=alpha) {
            newvarp <- propvarp
            LpreviousT <- Lprevious2
          } 
        }
      }		
      varp[i, 1:nbvar] <- newvarp
      varp[i, "Ln L"] <- LpreviousT
      
      if (MaxL["Ln L"] < varp[i, "Ln L"]) {
        MaxL <- varp[i, ]
      }
      
      
      # 6/10/2012: Je stocke tout	
      #	if ((i>n.adapt) && ((i%%thin)==0)) {
      varp2[cpt, 1:nbvar] <- newvarp
      varp2[cpt, "Ln L"] <- varp[i, "Ln L"]
      cpt <- cpt + 1
      #	}
      
      if (tf) {
        cpt_trace <- cpt_trace+1
        if (cpt_trace>=pt) {
          cat(paste("Chain ", kk, ": [", i, "] ", as.numeric(varp[i, "Ln L"]), "\n", sep=""))
          cpt_trace <- 0
        }
      }
      
      if (traceML) {
        if (as.numeric(varp[i, "Ln L"]) > previousML) {
          previousML <- varp[i, "Ln L"]
          cat(d(varp[i, 1:nbvar]))
          cat("\n")
        }
      }
      
      
      
      if (adaptive) {
        if ((i %% adaptive.lag) == 0) {
          
          #		    sdp <- 
          for (j in 1:nbvar) {
            dta <- varp[(i-adaptive.lag):i, j]
            txaccept <- sum(diff(dta)!=0)/(length(dta)-1)
            sdg[j] <- adaptive.fun(txaccept)*sdg[j]
          }
          # print(sdg)
        }
      }
      
    }
    
    lp <- getFromNamespace("as.mcmc", ns="coda")(varp2[(n.adapt+1):(cpt-1), 1:nbvar, drop=FALSE])
    lp <- getFromNamespace("mcmc", ns="coda")(data=lp, start=n.adapt+1, end=n.iter+n.adapt, thin=thin)
    
    # Ca c'est uniquement pour les chaînes
    # Ce n'est pas ca qui fait perdre du temps
    
    res<-c(res, list(lp))
    resL <- c(resL, list(varp2[(n.adapt+1):(cpt-1), "Ln L"]))
  }
  
  names(res) <- 1:n.chains
  res <- getFromNamespace("as.mcmc.list", ns="coda")(res)
  
  if (trace) {
  cat("Best likelihood for: \n")
  for (j in 1:nbvar) {cat(paste(names(MaxL[j]), "=", MaxL[j], "\n"))}
  }
  
  # 16/9/2019: Si un seul paramètre, renvoie quand même une matrix
  if (is.null(dim(res[[1]]))) {
    for (i in 1:length(res)) {
      res[[i]] <- matrix(data = res[[i]], ncol = 1, dimnames = list(c(NULL), rownames(parameters)))
    }
  }
  
  out <- (list(resultMCMC=res, resultLnL=resL, 
               parametersMCMC=list(parameters=parameters, n.iter=n.iter, n.chains=n.chains, n.adapt=n.adapt, thin=thin, 
                                   SDProp.end=structure(sdg, .Names=rownames(parameters)), 
                                   control=datax)))
  class(out) <- "mcmcComposite"
  return(out)
  
}
