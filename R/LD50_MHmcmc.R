#' LD50_MHmcmc runs the Metropolis-Hastings algorithm for LD50 (Bayesian MCMC)
#' @title Metropolis-Hastings algorithm for LD50
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
#' @param n.iter Number of iterations for each step
#' @param parametersMCMC A set of parameters used as initial point for searching with information on priors
#' @param result An object obtained after a SearchR fit
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param trace True or False, shows progress
#' @param batchSize Number of observations to include in each batch fo SE estimation
#' @param adaptive Should an adaptive process for SDProp be used
#' @param adaptive.lag  Lag to analyze the SDProp value in an adaptive content
#' @param adaptive.fun Function used to change the SDProp
#' @param filename If intermediate is not NULL, save intermediate result in this file
#' @param intermediate Period for saving intermediate result, NULL for no save
#' @param previous Previous result to be continued. Can be the filename in which intermediate results are saved.
#' @description Run the Metropolis-Hastings algorithm for tsd.\cr
#' Deeply modified from a MCMC script by Olivier Martin (INRA, Paris-Grignon).\cr
#' The number of iterations is n.iter+n.adapt+1 because the initial likelihood is also displayed.\cr
#' I recommend that thin=1 because the method to estimate SE uses resampling.\cr
#' If initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
#' To get the SE from result_mcmc <- tsd_MHmcmc(result=try), use:\cr
#' result_mcmc$BatchSE or result_mcmc$TimeSeriesSE\cr
#' The batch standard error procedure is usually thought to be not as accurate as the time series methods.\cr
#' Based on Jones, Haran, Caffo and Neath (2005), the batch size should be equal to sqrt(n.iter).\cr
#' Jones, G.L., Haran, M., Caffo, B.S. and Neath, R. (2006) Fixed Width Output Analysis for Markov chain Monte Carlo , Journal of the American Statistical Association, 101:1537-1547.\cr
#' coda package is necessary for this function.\cr
#' The parameters intermediate and filename are used to save intermediate results every 'intermediate' iterations (for example 1000). Results are saved in a file of name filename.\cr
#' The parameter previous is used to indicate the list that has been save using the parameters intermediate and filename. It permits to continue a mcmc search.\cr
#' These options are used to prevent the consequences of computer crash or if the run is very very long and processes at time limited.\cr
#' @family LD50 functions
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' data <- data.frame(Doses=c(80, 120, 150, 150, 180, 200),
#' Alive=c(10, 12, 8, 6, 2, 1),
#' Dead=c(0, 1, 5, 6, 9, 15))
#' LD50_logistic <- LD50(data, equation="logistic")
#' pMCMC <- LD50_MHmcmc_p(LD50_logistic, accept=TRUE)
#' # Take care, it can be very long
#' result_mcmc_LD50 <- LD50_MHmcmc(result=LD50_logistic, 
#' 		parametersMCMC=pMCMC, n.iter=10000, n.chains = 1,  
#' 		n.adapt = 0, thin=1, trace=1000, adaptive=TRUE, )
#' # summary() permits to get rapidly the standard errors for parameters
#' summary(result_mcmc_LD50)
#' plot(x=result_mcmc_LD50, parameters="S", scale.prior=TRUE, las=1)
#' plot(result_mcmc_LD50, parameters="S", scale.prior=TRUE, las=1, xlim=c(-20, 20))
#' plot(result_mcmc_LD50, parameters="P", scale.prior=TRUE, las=1)
#' 1-rejectionRate(as.mcmc(result_mcmc_LD50))
#' raftery.diag(as.mcmc(result_mcmc_LD50))
#' heidel.diag(as.mcmc(result_mcmc_LD50))
#' 
#' #### Example with Uniforms priors
#' 
#' pMCMC <- structure(list(Density = c("dunif", "dunif"), 
#' Prior1 = c(77.6216005852911, -31.0438095277258), 
#' Prior2 = c(310.486402341165, 31.0438095277258), 
#' SDProp = c(2, 0.5), 
#' Min = c(77.6216005852911, -31.0438095277258), 
#' Max = c(310.486402341165, 31.0438095277258), 
#' Init = c(155.243201170582, -15.5219047638629)), 
#' row.names = c("P", "S"), class = "data.frame")
#' result_mcmc_LD50 <- LD50_MHmcmc(result=LD50_logistic, 
#' 		parametersMCMC=pMCMC, n.iter=10000, n.chains = 1,  
#' 		n.adapt = 0, thin=1, trace=1000, adaptive=TRUE, )
#' # summary() permits to get rapidly the standard errors for parameters
#' summary(result_mcmc_LD50)
#' plot(x=result_mcmc_LD50, parameters="S", scale.prior=TRUE, las=1)
#' plot(result_mcmc_LD50, parameters="S", scale.prior=TRUE, las=1, xlim=c(-40, 40))
#' plot(result_mcmc_LD50, parameters="P", scale.prior=TRUE, las=1)
#' 1-rejectionRate(as.mcmc(result_mcmc_LD50))
#' raftery.diag(as.mcmc(result_mcmc_LD50))
#' heidel.diag(as.mcmc(result_mcmc_LD50))
#' 
#' }
#' @export

LD50_MHmcmc <- function(result=stop("A result of LD50() fit must be provided"), n.iter=10000, 
                        parametersMCMC=NULL, n.chains = 1, n.adapt = 0, 
                        thin=1, trace=FALSE, batchSize=sqrt(n.iter), 
                        adaptive=FALSE, adaptive.lag=500, 
                        adaptive.fun=function(x) {ifelse(x>0.234, 1.3, 0.7)},
                        intermediate=NULL, filename="intermediate.Rdata", previous=NULL) {
  
  # result=eo_logistic; parametersMCMC=NULL; 
  # n.iter=10000; n.chains = 1;  n.adapt = 0; thin=1; trace=TRUE; batchSize=sqrt(n.iter);intermediate=NULL; filename="intermediate.Rdata"; previous=NULL; adaptive=FALSE; adaptive.lag=500; adaptive.fun=function(x) {ifelse(x>0.234, 1.3, 0.7)}
  
  if (is.character(previous)) {
    itr <- NULL
    load(previous)
    previous <- itr
    rm(itr)
    print("Continue previous mcmc run")
  } else {
    print(parametersMCMC)
  }
  
  
  # 29/1/2014; Ajout de result$weight
  # 30/1/2015 Ajout de fixedparameters
  out <- MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, 
                   parameters_name = "par", 
                   n.chains = n.chains, n.adapt = n.adapt, thin=thin, trace=trace, 
                   alive=result$alive, N=result$N, doses=result$doses, 
                   equation=result$equation, fixed.parameters=result$fixed.parameters, 
                   likelihood= getFromNamespace(".LD50_fit", ns="HelpersMG"), 
                   adaptive=adaptive, adaptive.lag=adaptive.lag, adaptive.fun=adaptive.fun,
                   intermediate=intermediate, filename=filename, previous=previous)
  
  if (batchSize>=n.iter/2) {
    print("batchSize cannot be larger than half the number of iterations.")
    rese <- rep(NA, dim(parametersMCMC)[1])
    names(rese) <- rownames(parametersMCMC)
    out <- c(out, SE=list(rese))
  } else {
    out <- c(out, BatchSE=list(coda::batchSE(out$resultMCMC, batchSize=batchSize)))
  }
  
  # pas clair si c'est utile
  out <- addS3Class(out, "mcmcComposite")
  
  fin <- try(summary(out), silent=TRUE)
  
  if (inherits(fin, "try-error")) {
    lp <- rep(NA, nrow(out$parametersMCMC$parameters))
    names(lp) <- rownames(out$parametersMCMC$parameters)
    out <- c(out, TimeSeriesSE=list(lp))
    out <- c(out, SD=list(lp))
  } else {
    out <- c(out, TimeSeriesSE=list(fin$statistics[,4]))
    out <- c(out, SD=list(fin$statistics[,"SD"]))
  }
  
  out <- addS3Class(out, "mcmcComposite")
  
  return(out)
                        }
