#' from_min_max returns standard deviation and/or mean from minimum and maximum
#' @title Distribution from minimum and maximum
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return from_min_max returns a list with output, ML, Bayesian results
#' @param observed.Minimum The observed minimum
#' @param observed.Maximum The observed maximum
#' @param observed.Mean The observed mean (can be omitted)
#' @param observed.Median The observed median (can be omitted)
#' @param observed.Quantiles Observed quantiles with names being the values of quantiles with Q being first letter (can be omitted; see description)
#' @param priors The priors (see MHAlgogen()) or character "dnorm" or "dunif"
#' @param n Number of observations
#' @param fitted.parameters The initial value to fit
#' @param fixed.parameters The fixed parameters
#' @param D The distribution to fit as character. ex "norm" or "lnorm"
#' @param n.iter Number of iterations for each chain
#' @param n.chains Number of chains
#' @param n.adapt Number of iteration to stabilize likelihood
#' @param thin Interval for thinning likelihoods
#' @param adaptive Should an adaptive process for SDProp be used
#' @param trace Or FALSE or period to show progress
#' @param replicates Number of replicates to model D
#' @param silent If TRUE, do not print information
#' @description Bayesian estimate of distribution when minimum, maximum, median or mean are known.\cr
#' If D="norm**" or "lnorm**", it will use the approximation of Gumbel based on D.\cr
#' If D="norm*" or "lnorm*", it will generate D distribution using replicates number of random numbers, 
#' and estimate Gumbel parameters from the simulated D distribution. \cr
#' Otherwise it will estimate parameters of Gumbel distribution based on maximum likelihood.\cr
#' For D="pois", "beta" or "chisq" only second (D="pois*", "beta*" or "chisq*") and third solutions are available.\cr
#' The observed.Quantiles parameter must be a named value, for example observed.Quantiles = c(Q0.025=17, Q0.975=26)\cr
#' @examples
#' \dontrun{
#' minobs <- 5
#' maxobs <- 25
#' # These two values are not mandatory
#' meanobs <- 15
#' medianobs <- 16
#' n <- 10
#' # To estimate only the sd of the distribution; mean is fixed
#' # Note that there is an obligation to have mean even if it
#' # is in fixed.parameters
#' # By defaut a normal distribution is fitted
#' 
#' out_sd <- from_min_max(n=n                                , 
#'                        observed.Minimum=minobs                    , 
#'                        observed.Maximum=maxobs                    , 
#'                        fitted.parameters=c(sd=5)                  , 
#'                        fixed.parameters=c(mean=meanobs)           ,
#'                        n.iter=10000                               ,
#'                        trace=TRUE                                 )
#'                             
#' plot(out_sd, what="MarkovChain", parameters="sd")
#' plot(out_sd, what="posterior", parameters="sd")
#' as.parameters(out_sd, index="quantile")
#' 
#' # To estimate both the sd and mean of the distribution
#' 
#' out_sd_mean_norm <- from_min_max(n=n, 
#'                             observed.Minimum=minobs                    , 
#'                             observed.Maximum=maxobs                    , 
#'                             fitted.parameters=c(mean=15, sd=5)         , 
#'                             fixed.parameters=NULL                      ,
#'                             n.iter=10000                               ,
#'                             D="norm**"                                 ,
#'                             trace=FALSE                                )
#'                             
#' plot(out_sd_mean_norm, what="MarkovChain", parameters="sd", ylim=c(0, 10))
#' plot(out_sd_mean_norm, what="MarkovChain", parameters="mean", ylim=c(10, 20))
#' plot(out_sd_mean_norm, what="posterior", parameters="mean",  
#'      breaks=seq(from=0, to=100, by=1), xlim=c(10, 20))
#' as.parameters(out_sd_mean_norm, index="quantile")
#' 
#' # Let see what's happened for a lognormal distribution
#' 
#' out_sd_mean_lnorm <- from_min_max(n=n, 
#'                             observed.Minimum=minobs                    , 
#'                             observed.Maximum=maxobs                    , 
#'                             fitted.parameters=c(meanlog=15, sdlog=5)   , 
#'                             fixed.parameters=NULL                      ,
#'                             n.iter=10000                               ,
#'                             D="lnorm**"                                ,
#'                             trace=FALSE                                )
#'                             
#' plot(out_sd_mean_lnorm, what="MarkovChain", parameters="sdlog", ylim=c(0, 10))
#' plot(out_sd_mean_lnorm, what="MarkovChain", parameters="meanlog", ylim=c(10, 20))
#' plot(out_sd_mean_lnorm, what="posterior", parameters="meanlog", xlim=c(0, 20), 
#'      breaks=seq(from=0, to=100, by=1))
#' as.parameters(out_sd_mean_lnorm, index="quantile")
#' 
#' # To be compared with the rule of thumb:
#' print(paste0("mean = ", as.character((maxobs + minobs) / 2))) # Mean Not so bad
#' print(paste0("sd = ", as.character((maxobs - minobs) / 4))) # SD Clearly biased
#' 
#' # Covariation of sd and mean is nearly NULL
#' cor(x=as.parameters(out_sd_mean_norm, index="all")[, "mean"], 
#'     y=as.parameters(out_sd_mean_norm, index="all")[, "sd"])^2
#' plot(x=as.parameters(out_sd_mean_norm, index="all")[, "mean"], 
#'      y=as.parameters(out_sd_mean_norm, index="all")[, "sd"], 
#'      xlab="mean", ylab="sd")
#'      
#' # Example when minimum, maximum and mean are known
#' 
#' out_sd_mean2 <- from_min_max(n=n                                   , 
#'                              observed.Minimum=minobs               , 
#'                              observed.Maximum=maxobs               , 
#'                              observed.Mean=meanobs                 ,
#'                              fitted.parameters=c(mean=15, sd=5)    , 
#'                              fixed.parameters=NULL                 ,
#'                              n.iter=10000                          ,
#'                              trace=FALSE                            )
#'                             
#' # Example when minimum, maximum, mean and median are known
#' 
#' out_sd_mean3 <- from_min_max(n=n                                   , 
#'                              observed.Minimum=minobs               , 
#'                              observed.Maximum=maxobs               , 
#'                              observed.Mean=meanobs                 ,
#'                              observed.Median=medianobs             ,
#'                              fitted.parameters=c(mean=15, sd=5)    , 
#'                              fixed.parameters=NULL                 ,
#'                              n.iter=10000                          ,
#'                              trace=FALSE                            )
#'                              
#' plot(out_sd_mean2, what="MarkovChain", parameters="sd")
#' plot(out_sd_mean2, what="MarkovChain", parameters="mean")
#' plot(out_sd_mean2, what="posterior", parameters="mean", xlim=c(0, 100), 
#'      breaks=seq(from=0, to=100, by=5))
#' as.parameters(out_sd_mean2, index="quantile")
#' 
#' # Example of GEV density function  
#' # Parametrisation from https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution
#' dGEV <- getFromNamespace(".dGEV", ns="HelpersMG")
#' x <- seq(from=-4, to=4, by=0.1)
#' plot(x, y=dGEV(x=x, 
#'                location=0, scale=1, shape=-1/2, log=FALSE, sum=FALSE), 
#'      type="l", col="green", xlab="x", ylab="Density")
#' lines(x, y=dGEV(x=x, 
#'                 location=0, scale=1, shape=0, log=FALSE, sum=FALSE), col="red")
#' lines(x, y=dGEV(x=x, 
#'                 location=0, scale=1, shape=1/2, log=FALSE, sum=FALSE), col="blue")
#' legend("topleft", legend=c("shape=-1/2", "shape=0", "shape=1/2"), 
#'        lty=1, col=c("green", "red", "blue"))
#'        
#' # Note the different parametrisation about shape
#' dGEV <- getFromNamespace("dgevd", ns="EnvStats")
#' x <- seq(from=-4, to=4, by=0.1)
#' plot(x, y=dGEV(x=x, 
#'                location=0, scale=1, shape=-1/2), 
#'      type="l", col="green", xlab="x", ylab="Density")
#' lines(x, y=dGEV(x=x, 
#'                 location=0, scale=1, shape=0), col="red")
#' lines(x, y=dGEV(x=x, 
#'                 location=0, scale=1, shape=1/2), col="blue")
#' legend("topleft", legend=c("shape=-1/2", "shape=0", "shape=1/2"), 
#'        lty=1, col=c("green", "red", "blue"))
#' 
#' # Compute dn using the approximation from Wan et al. (2014)
#'    get_dn <- function(n) {
#'      if (n < 2) {
#'        stop("Sample size n must be at least 2.")
#'      }
#'      qnorm((n - 0.375) / (n + 0.25)) * 2
#'    }
#'    
#'    # Estimate standard deviation from min and max
#'    estimate_sd_from_range <- function(min_val, max_val, n) {
#'      dn <- get_dn(n)
#'      range <- max_val - min_val
#'      sd_estimate <- range / dn
#'      return(sd_estimate)
#'    }
#'    
#'    # Example usage:
#'    n <- 10
#'    min_val <- 5
#'    max_val <- 25
#'    
#'    dn_value <- get_dn(n)
#'    sd_estimate <- estimate_sd_from_range(min_val, max_val, n)
#'    
#'    cat("dn =", dn_value, "\n")
#'    cat("Estimated SD =", sd_estimate, "\n")
#'    
#'    # To generate data from publication
#'    
#'    library(parallel)
#'    library(embryogrowth)
#'    
#'    # Values for the prior of SD
#'    outSD <- subset(DatabaseTSD, subset = (((!is.na(IP.SD)) | 
#'                    (!is.na(IP.SE))) & (!is.na(Hatched))), 
#'                    select=c("Hatched", "IP.SE", "IP.SD", "IP.mean"))
#'    outSD$IP.SD <- ifelse(is.na(outSD$IP.SD), outSD$IP.SE*sqrt(outSD$Hatched), outSD$IP.SD)
#'    
#'    # Model estimation
#'    
#'    Example <- subset(DatabaseTSD, subset = (!is.na(IP.min)) & 
#'                        ((is.na(IP.SE)) & (is.na(IP.SD)) & (!is.na(Hatched)) & 
#'                           (is.na(IP.mean))), select=c("Species", "Incubation.temperature.set", 
#'                                                       "Hatched", "IP.min", "IP.max", 
#'                                                       "Reference"))
#'    
#'    out <- universalmclapply(X=1:nrow(Example), FUN=function(i) {
#'      n <- Example[i, "Hatched"]
#'      
#'      priors <- structure(list(Density = c("dunif", "dlnorm"), 
#'                               Prior1 = c(30, log(mean(outSD$IP.SD))), 
#'                               Prior2 = c(120, log(sd(outSD$IP.SD))), 
#'                               SDProp = c(1, 1), 
#'                               Min = c(30, 0.1), 
#'                               Max = c(120, 6), 
#'                               Init = c((Example[i, "IP.min"]+ Example[i, "IP.max"])/2, log(2))), 
#'                          row.names = c("mean", "sd"), 
#'                          class = c("PriorsmcmcComposite", "data.frame"))
#'      
#'      
#'      out_sd_mean_mcmc <- from_min_max(n=n, observed.Minimum=Example[i, "IP.min"], 
#'                                       observed.Maximum=Example[i, "IP.max"], 
#'                                       fitted.parameters=c(mean=(Example[i, "IP.min"]+ 
#'                                                              Example[i, "IP.max"])/2, 
#'                                                           sd=log(2)), 
#'                                       priors = priors, 
#'                                       D="norm**", 
#'                                       n.iter = 10000, n.adapt=15000, thin=30, 
#'                                       trace=100, adaptive = TRUE)
#'      
#'      # plot(out_sd_mean_mcmc, what = "MarkovChain", parameters = "sd")
#'      
#'      assign(paste0("out_sd_mean_mcmc_", as.character(i)), out_sd_mean_mcmc)
#'      save(list = paste0("out_sd_mean_mcmc_", as.character(i)), 
#'           file = file.path("dataOut", paste0("out_sd_mean_mcmc_", as.character(i), ".Rdata")))
#'      # rm(list=paste0("out_sd_mean_mcmc_", as.character(i)))
#'    }, progressbar = TRUE)
#'    
#'    # Generate table with all results
#'    
#'    Example <- subset(DatabaseTSD, subset = (!is.na(IP.min)) & 
#'    ((is.na(IP.SE)) & (is.na(IP.SD)) & (!is.na(Hatched)) & 
#'      (is.na(IP.mean))), select=c("Species", "Incubation.temperature.set", 
#'                                  "Hatched", "IP.min", "IP.max", 
#'                                  "Reference"))
#'    
#'    Example <- cbind(Example, dn=NA)
#'    Example <- cbind(Example, "SD(Hozo 2005)"=NA)
#'    Example <- cbind(Example, "SD(Wan 2014)"=NA)
#'    Example <- cbind(Example, "mean(Wan 2014)"=NA)
#'    Example <- cbind(Example, "median(SD)"=NA)
#'    Example <- cbind(Example, "2.5%(SD)"=NA)
#'    Example <- cbind(Example, "97.5%(SD)"=NA)
#'    Example <- cbind(Example, "25%(SD)"=NA)
#'    Example <- cbind(Example, "75%(SD)"=NA)
#'    Example <- cbind(Example, "median(mean)"=NA)
#'    Example <- cbind(Example, "2.5%(mean)"=NA)
#'    Example <- cbind(Example, "97.5%(mean)"=NA)
#'    Example <- cbind(Example, "25%(mean)"=NA)
#'    Example <- cbind(Example, "75%(mean)"=NA)
#'    Example <- cbind(Example, "z(mean)"=NA)
#'    Example <- cbind(Example, "z(SD)"=NA)
#'    
#'    
#'    library(coda)
#'    
#'    for (i in 1:nrow(Example)) {
#'      
#'      n <- Example[i, "Hatched"]
#'      if (n<= 15) {
#'        Example[i, "SD(Hozo 2005)"] <- (1/sqrt(12))*sqrt(((Example[i, "IP.max"] - 
#'         Example[i, "IP.min"]))^2+((Example[i, "IP.max"] - Example[i, "IP.min"]))^2/4)
#'      } else {
#'        if (n<=70) {
#'          Example[i, "SD(Hozo 2005)"] <- (Example[i, "IP.max"] - Example[i, "IP.min"])/4
#'        } else {
#'          Example[i, "SD(Hozo 2005)"] <- (Example[i, "IP.max"] - Example[i, "IP.min"])/6
#'        }
#'      }
#'      
#'      Example[i, "dn"] <- get_dn(n)
#'      Example[i, "SD(Wan 2014)"] <- estimate_sd_from_range(Example[i, "IP.min"], 
#'             Example[i, "IP.max"], n)
#'      Example[i, "mean(Wan 2014)"] <- (Example[i, "IP.min"]+ Example[i, "IP.max"])/2
#'      load(file = file.path("dataOut", paste0("out_sd_mean_mcmc_", as.character(i), ".Rdata")))
#'      out_sd_mean_mcmc <- get(paste0("out_sd_mean_mcmc_", as.character(i)))
#'      k <- as.parameters(out_sd_mean_mcmc, index="quantile", probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
#'      outgk <- geweke.diag(as.mcmc(out_sd_mean_mcmc))
#'      rm(out_sd_mean_mcmc)
#'      rm(list=paste0("out_sd_mean_mcmc_", as.character(i)))
#'      Example[i, "median(SD)"] <- k["50%", "sd"]
#'      Example[i, "2.5%(SD)"] <- k["2.5%", "sd"]
#'      Example[i, "97.5%(SD)"] <- k["97.5%", "sd"]
#'      Example[i, "25%(SD)"] <- k["25%", "sd"]
#'      Example[i, "75%(SD)"] <- k["75%", "sd"]
#'      Example[i, "median(mean)"] <- k["50%", "mean"]
#'      Example[i, "2.5%(mean)"] <- k["2.5%", "mean"]
#'      Example[i, "97.5%(mean)"] <- k["97.5%", "mean"]
#'      Example[i, "25%(mean)"] <- k["25%", "mean"]
#'      Example[i, "75%(mean)"] <- k["75%", "mean"]
#'      Example[i, "z(mean)"] <- outgk$`1`$z["mean"]
#'      Example[i, "z(SD)"] <- outgk$`1`$z["sd"]
#'    }
#'    
#'    rownames(Example) <- as.character(1:nrow(Example))
#'   
#'   # Figure 1
#'   layout(mat = matrix(1:4, nrow=2))
#'   par(mar=c(4, 4, 0, 0))
#'   plot(out_sd_mean_mcmc_11, what = "MarkovChain", parameters = "mean", ylim=c(70, 76))
#'   text(x=ScalePreviousPlot(x=0.05, y=0.95)$x, y=ScalePreviousPlot(x=0.85, y=0.95)$y, 
#'      labels = "A: mean", cex=1.5, pos=4)
#'   plot(out_sd_mean_mcmc_8, what = "MarkovChain", parameters = "mean", ylim=c(50, 52))
#'   text(x=ScalePreviousPlot(x=0.05, y=0.95)$x, y=ScalePreviousPlot(x=0.85, y=0.95)$y, 
#'       labels = "C: mean", cex=1.5, pos=4)
#'   plot(out_sd_mean_mcmc_11, what = "MarkovChain", parameters = "sd")
#'   text(x=ScalePreviousPlot(x=0.05, y=0.95)$x, y=ScalePreviousPlot(x=0.85, y=0.95)$y, 
#'       labels = "B: sd", cex=1.5, pos=4)
#'   plot(out_sd_mean_mcmc_8, what = "MarkovChain", parameters = "sd")
#'   text(x=ScalePreviousPlot(x=0.05, y=0.95)$x, y=ScalePreviousPlot(x=0.85, y=0.95)$y, 
#'       labels = "D: sd", cex=1.5, pos=4)
#'   
#'   # Figure 2
#'   dtafigure2 <- matrix(NA, nrow=nrow(as.parameters(out_sd_mean_mcmc, index = "all")), 
#'                        ncol=nrow(Example))
#'   
#'   for (i in 1:nrow(Example)) {
#'     # i <- 1
#'     load(file=file.path("dataOut", paste0("out_sd_mean_mcmc_", as.character(i), ".Rdata")))
#'     out_sd_mean_mcmc <- get(paste0("out_sd_mean_mcmc_", as.character(i)))
#'     PPM <- rnorm(nrow(as.parameters(out_sd_mean_mcmc, index = "all")), 
#'        mean = as.parameters(out_sd_mean_mcmc, index = "all")[, "mean"], 
#'        sd=as.parameters(out_sd_mean_mcmc, index = "all")[, "sd"])
#'     dtafigure2[, i] <- PPM
#'   }
#'   
#'   layout(mat = 1)
#'   par(mar=c(3, 4, 0, 0))
#'   boxplot(dtafigure2, outline=FALSE, las=1, bty="n", xaxt="n", frame=FALSE, ylim=c(40, 90), 
#'     col=sapply(as.character(Example$Species), 
#'       FUN=function(i) switch(i, "Caretta caretta"="white", 
#'                                  "Chelonia mydas"="green", 
#'                                  "Dermochelys coriacea"="lightblue")), 
#'     ylab="Incubation duration in days")
#'   axis(1, at=1:30, labels = rep(NA, 30))
#'   Cc <- sum(as.character(Example$Species) == "Caretta caretta")
#'   CcCm <- sum(as.character(Example$Species) == "Caretta caretta" | 
#'               as.character(Example$Species) == "Chelonia mydas")
#'   segments(x0= Cc + 0.5, 
#'            x1= Cc + 0.5, 
#'            y0=40, y1=90, lty = 2)
#'   
#'   segments(x0= CcCm + 0.5, 
#'            x1= CcCm + 0.5, 
#'            y0=40, y1=90, lty = 2)
#'   
#'   par(xpd=TRUE)
#'   text(x=Cc/2, y=33, labels = expression(italic("Caretta caretta")), cex=0.9)
#'   
#'   text(x=Cc+(CcCm-Cc)/2, y=32, labels = expression(italic("Chelonia\n  mydas")), cex=0.9)
#'   text(x=CcCm+(30-CcCm)/2, y=32, labels = expression(italic("Dermochelys\n  coriacea")), cex=0.9)
#'   
#'   # With Lognormal
#'    # To generate data from publication
#'    
#'    library(parallel)
#'    library(embryogrowth)
#'    
#'    # Values for the prior of SD
#'    outSD <- subset(DatabaseTSD, subset = (((!is.na(IP.SD)) | 
#'                    (!is.na(IP.SE))) & (!is.na(Hatched))), 
#'                    select=c("Hatched", "IP.SE", "IP.SD", "IP.mean"))
#'    outSD$IP.SD <- ifelse(is.na(outSD$IP.SD), outSD$IP.SE*sqrt(outSD$Hatched), outSD$IP.SD)
#'    
#'    # Model estimation
#'    
#'    Example <- subset(DatabaseTSD, subset = (!is.na(IP.min)) & 
#'                        ((is.na(IP.SE)) & (is.na(IP.SD)) & (!is.na(Hatched)) & 
#'                           (is.na(IP.mean))), select=c("Species", "Incubation.temperature.set", 
#'                                                       "Hatched", "IP.min", "IP.max", 
#'                                                       "Reference"))
#'    
#'    out <- universalmclapply(X=1:nrow(Example), FUN=function(i) {
#'      n <- Example[i, "Hatched"]
#'      
#'      priors <- structure(list(Density = c("dunif", "dlnorm"), 
#'                               Prior1 = c(30, log(mean(outSD$IP.SD))), 
#'                               Prior2 = c(120, log(sd(outSD$IP.SD))), 
#'                               SDProp = c(1, 1), 
#'                               Min = c(30, 0.1), 
#'                               Max = c(120, 6), 
#'                               Init = c((Example[i, "IP.min"]+ Example[i, "IP.max"])/2, log(2))), 
#'                          row.names = c("meanlog", "sdlog"), 
#'                          class = c("PriorsmcmcComposite", "data.frame"))
#'      
#'      
#'      out_sd_mean_mcmc_LN <- from_min_max(n=n, observed.Minimum=Example[i, "IP.min"], 
#'                                       observed.Maximum=Example[i, "IP.max"], 
#'                                       fitted.parameters=c(mean=(Example[i, "IP.min"]+ 
#'                                                              Example[i, "IP.max"])/2, 
#'                                                           sd=log(2)), 
#'                                       priors = priors, 
#'                                       D="lnorm**", 
#'                                       n.iter = 10000, n.adapt=15000, thin=30, 
#'                                       trace=100, adaptive = TRUE)
#'      
#'      # plot(out_sd_mean_mcmc, what = "MarkovChain", parameters = "sd")
#'      
#'      assign(paste0("out_sd_mean_mcmc_LN_", as.character(i)), out_sd_mean_mcmc_LN)
#'      save(list = paste0("out_sd_mean_mcmc_LN_", as.character(i)), 
#'           file = file.path("dataOut", paste0("out_sd_mean_mcmc_LN_", as.character(i), ".Rdata")))
#'      # rm(list=paste0("out_sd_mean_mcmc_LN_", as.character(i)))
#'    }, progressbar = TRUE)
#'    
#' }
#' @export

from_min_max <- function(n = stop("n must be known.")                                             , 
                         fitted.parameters=stop("At least one parameter must be supplied.")       , 
                         observed.Minimum                                                         , 
                         observed.Maximum                                                         , 
                         observed.Median=NULL                                                     ,
                         observed.Mean=NULL                                                       , 
                         observed.Quantiles=NULL                                                  , 
                         priors="dnorm"                                                           , 
                         fixed.parameters=NULL                                                    ,
                         D="norm**"                                                               ,
                         n.iter=10000                                                             ,
                         n.chains = 1                                                             ,
                         n.adapt = 100                                                            ,
                         thin=30                                                                  ,
                         adaptive = FALSE                                                         , 
                         trace = 100                                                              , 
                         replicates = 10000                                                       , 
                         silent = FALSE                                                           ) {
  
  
  # x sont les paramètres
  L_min_max <- function(x=NULL                               , 
                        fixed.parameters=NULL                , 
                        n=NULL                               , 
                        D="norm"                             ,
                        meanobs=NULL                         , 
                        medianobs=NULL                       ,
                        minobs=NULL                          , 
                        maxobs=NULL                          , 
                        quantilesobs=NULL                    ,
                        rep=10000                            ) {
    
    # J'envoie les paramètres de D (* ou **)
    # ainsi que les valeurs observées meanobs, medianobs, minobs, maxobs, quantilesobs
    # Et je retourne la vraisemblance
    
    # x = c(mean=100, sd=10); fixed.parameters=NULL; n=10; D="norm**"; meanobs=100; medianobs=100; minobs=90; maxobs=110; quantilesobs=c(Q0.95=106); rep=10000      
    
    
    # par <- c(mean=50, sd=3)
    # n <- 4; minobs <- 43; maxobs <- 57
    # print(fixed.parameters)
    x <- c(x, fixed.parameters)
    # print(d(x))
    
    # if (is.null(meanobs)) {
    #   meanobs <- (minobs + maxobs) / 2
    # }
    
    D_ec <- D
    if (D_ec == "norm**") D <- "norm"
    if (D_ec == "norm*") D <- "norm"
    if (D_ec == "lnorm**") D <- "lnorm"
    if (D_ec == "lnorm*") D <- "lnorm"
    
    Lmin <- Lmax <- NULL
    kkk <- NULL
    
    # Si norm** et lnorm**, je calcule Lmin et Lmax sur l'approximation
    
    if (D_ec == "norm**") {
      # x <- c(mean=50, sd=2)
      # hist(rnorm(1000, mean=x["mean"], sd=x["sd"]))
      # n <- 10
      # minobs <- 42
      # maxobs <- 56
      # D <- "norm"; rep <- 10000
      
      mean <- - x["mean"]
      sd <- x["sd"]
      scale <- sd/sqrt(2*log(n))
      location <- mean + sd*sqrt(2*log(n)) - ((log(log(n))+log(4*pi))*sd/(2*sqrt(2*log(n))))
      Lmin <- - getFromNamespace(".dGEV", ns="HelpersMG")(x = - minobs, location= location, scale= scale, shape=0, log=TRUE, sum=TRUE)
      
      mean <- x["mean"]
      # sd <- x["sd"] # J'ai le même SD
      # scale <- sd/sqrt(2*log(n)) # J'ai le même scale
      location <- mean + sd*sqrt(2*log(n)) - ((log(log(n))+log(4*pi))*sd/(2*sqrt(2*log(n))))
      Lmax <- - getFromNamespace(".dGEV", ns="HelpersMG")(x = maxobs, location= location, scale= scale, shape=0, log=TRUE, sum=TRUE)
    }
    
    if (D_ec == "lnorm**") {
      # x <- c(meanlog=log(50), sdlog=log(2))
      # hist(log(rlnorm(1000, meanlog=x["meanlog"], sdlog=x["sdlog"])))
      # hist(rnorm(1000, mean=x["meanlog"], sd=x["sdlog"]))
      # n <- 10
      # minobs <- 47
      # maxobs <- 54
      # D_ec <- "lnorm**"; D <- "lnorm"; rep <- 100000
      
      mean <- x["meanlog"]
      sd <- x["sdlog"]
      location <- mean - sd * sqrt(2*log(n))
      scale <- sd/sqrt(2*log(n))
      
      Lmin <- - getFromNamespace(".dGEV", ns="HelpersMG")(x= minobs, location=location, scale=scale, shape=0, log=TRUE, sum=TRUE)
      
      mean <- x["meanlog"]
      sd <- x["sdlog"]
      location <- mean + sd * sqrt(2*log(n))
      scale <- sd/sqrt(2*log(n))
      
      Lmax <- - getFromNamespace(".dGEV", ns="HelpersMG")(x=maxobs, location=location, scale=scale, shape=0, log=TRUE, sum=TRUE)
    }
    
    if ((D_ec == "norm*") | (D_ec == "lnorm*") | (D_ec == "pois*") | (D_ec == "beta*") | (D_ec == "chisq*")) {
      
      rD <- get(paste0("r", D))
      dsit <- sapply(1:rep, FUN=function(i) {k <- do.call(rD, args = as.list(c(n=n, x))); return(c(max(-k), max(k)))}) # mean(k), median(k), 
      kk <- matrix(dsit, ncol=2, byrow = TRUE)
      
      kkk <- apply(kk, MARGIN = 2, FUN=function(x) return(c(sd(x), median(x)))) # mean(x), 
      # colnames(kkk) <- c("min", "max")
      # rownames(kkk) <- c("sd", "median")
      
      # Je dois aussi sortir les quantiles si nécessaire
      
      
      # v <- kkk[2, 1]^2
      s <- kkk[1, 1] # sd de min
      # m <- kkk[1, 1]
      md <- kkk[2, 1] # median de min
      # print(c(m, md, v))
      # kpi <- 0.607927101854027 # (6/pi^2)
      sqkpi <- 0.779696801233676 # sqrt((6/pi^2))
      kll2 <- -0.366512920581664 # log(log(2))
      
      # sigma^2*(pi^2/6) = v
      # sigma^2 = v*(6/pi^2)
      scale <- s * sqkpi # sqrt(v*kpi)
      
      # md = mu - sigma * log(log(2))
      location <- md + scale * kll2
      
      shape <- 0
      # dGEV <- getFromNamespace(".dGEV", ns="HelpersMG")
      # dGEV(x=kk[, 1], par=c(location=location, scale=scale, shape=0), sum=TRUE, log=TRUE)
      
      # print(c(location=location, scale=scale, shape=0))
      # outGEV_min <- optim(par=c(location=location, scale=scale, shape=0), 
      #                     x=kk[, 1], 
      #                     sum=TRUE, log=TRUE, silent=TRUE, 
      #                     fn = getFromNamespace(".dGEV", ns="HelpersMG"), 
      #                     # lower = c(location=-Inf, scale=1E-6, shape=-Inf), 
      #                     # upper = c(location=Inf, scale=Inf, shape=Inf), 
      #                     method = "Nelder-Mead")$par
      # outGEV_min[2] <- abs(outGEV_min[2])
      Lmin <- - getFromNamespace(".dGEV", ns="HelpersMG")(x=-minobs, location=location, scale=scale, shape=0, log=TRUE, sum=TRUE)
      # print(outGEV_min)
      # print(Lmin)
      # v <- kkk[2, 4]^2
      s <- kkk[1, 2] # sd de max
      # m <- kkk[1, 4]
      md <- kkk[2, 2] # median de max
      
      # sigma^2*(pi^2/6) = v
      # sigma^2 = v*(6/pi^2)
      scale <- s * sqkpi # sqrt(v*kpi)
      
      # md = mu - sigma * log(log(2))
      location <- md + scale * kll2
      
      # outGEV_max <- optim(par=c(location=location, scale=scale, shape=0), 
      #                     x=kk[, 3], 
      #                     sum=TRUE, log=TRUE, silent=TRUE, 
      #                     fn = getFromNamespace(".dGEV", ns="HelpersMG"), 
      #                     # lower = c(location=-Inf, scale=1E-6, shape=-Inf), 
      #                     # upper = c(location=Inf, scale=Inf, shape=Inf), 
      #                     method = "Nelder-Mead")$par
      # outGEV_max[2] <- abs(outGEV_max[2])
      Lmax <- - getFromNamespace(".dGEV", ns="HelpersMG")(x=maxobs, location=location, scale=scale, shape=0, log=TRUE, sum=TRUE)
      
    }
    
    # Ca signifie que je n'ai rien trouvé
    if (is.null(Lmin)) {
      
      rD <- get(paste0("r", D))
      dsit <- sapply(1:rep, FUN=function(i) {k <- do.call(rD, args = as.list(c(n=n, x))); return(c(max(-k), max(k)))}) # mean(k), median(k), 
      kk <- matrix(dsit, ncol=2, byrow = TRUE)
      
      kkk <- apply(kk, MARGIN = 2, FUN=function(x) return(c(sd(x), median(x)))) # mean(x), 
      # colnames(kkk) <- c("min", "max")
      # rownames(kkk) <- c("sd", "median")
      
      # v <- kkk[2, 1]^2
      s <- kkk[1, 1] # sd de min
      # m <- kkk[1, 1]
      md <- kkk[2, 1] # median de min
      # print(c(m, md, v))
      # kpi <- 0.607927101854027 # (6/pi^2)
      sqkpi <- 0.779696801233676 # sqrt((6/pi^2))
      kll2 <- -0.366512920581664 # log(log(2))
      
      # sigma^2*(pi^2/6) = v
      # sigma^2 = v*(6/pi^2)
      scale <- s * sqkpi # sqrt(v*kpi)
      
      # md = mu - sigma * log(log(2))
      location <- md + scale * kll2
      
      shape <- 0
      # dGEV <- getFromNamespace(".dGEV", ns="HelpersMG")
      # dGEV(x=kk[, 1], par=c(location=location, scale=scale, shape=0), sum=TRUE, log=TRUE)
      
      # print(c(location=location, scale=scale, shape=0))
      outGEV_min <- optim(par=c(location=location, scale=scale, shape=0),
                          x=kk[, 1],
                          sum=TRUE, log=TRUE, silent=TRUE,
                          fn = getFromNamespace(".dGEV", ns="HelpersMG"),
                          # lower = c(location=-Inf, scale=1E-6, shape=-Inf),
                          # upper = c(location=Inf, scale=Inf, shape=Inf),
                          method = "Nelder-Mead")$par
      outGEV_min[2] <- abs(outGEV_min[2])
      Lmin <- - getFromNamespace(".dGEV", ns="HelpersMG")(x=-minobs, par=outGEV_min, log=TRUE, sum=TRUE)
      # print(outGEV_min)
      # print(Lmin)
      # v <- kkk[2, 1]^2
      s <- kkk[1, 2] # sd de max
      # m <- kkk[1, 4]
      md <- kkk[2, 2] # median de max
      # print(c(m, md, v))
      # kpi <- 0.607927101854027 # (6/pi^2)
      # sqkpi <- 0.779696801233676 # sqrt((6/pi^2))
      # kll2 <- -0.366512920581664 # log(log(2))
      
      # sigma^2*(pi^2/6) = v
      # sigma^2 = v*(6/pi^2)
      scale <- s * sqkpi # sqrt(v*kpi)
      
      # md = mu - sigma * log(log(2))
      location <- md + scale * kll2
      
      outGEV_max <- optim(par=c(location=location, scale=scale, shape=0),
                          x=kk[, 3],
                          sum=TRUE, log=TRUE, silent=TRUE,
                          fn = getFromNamespace(".dGEV", ns="HelpersMG"),
                          # lower = c(location=-Inf, scale=1E-6, shape=-Inf),
                          # upper = c(location=Inf, scale=Inf, shape=Inf),
                          method = "Nelder-Mead")$par
      outGEV_max[2] <- abs(outGEV_max[2])
      Lmax <- - getFromNamespace(".dGEV", ns="HelpersMG")(x=maxobs, par=outGEV_max, log=TRUE, sum=TRUE)
      
    }
    
    Lmean <- NULL
    Lmedian <- NULL
    Lquantiles <- NULL
    
    if (D_ec == "norm**") {
      # Je suis avec norm** 
      
      mean <- x["mean"]
      sd <- x["sd"]
      
      if (!is.null(meanobs)) {
        Lmean <- dnorm(x=meanobs, mean = mean, sd = sd/sqrt(n), log = TRUE)
      } else {
        Lmean <- 0
      }
      
      # Si median et norm
      if (!is.null(medianobs)) {
        
        likelihood_median_exact <- function(M, n, mean, sd) {
          if (n %% 2 == 0) {
            # --- Approximate likelihood (normal approximation) ---
            sd_median <- sd * sqrt(pi / (2 * n))
            logL <- dnorm(M, mean, sd_median, log = TRUE)
          } else {
            # --- Exact likelihood for odd n ---
            k <- (n - 1) / 2
            logf <- dnorm(M, mean, sd, log=TRUE)
            F <- pnorm(M, mean, sd)
            logcoef <- log(factorial(n)) - 2*log(factorial(k))
            logL <- logcoef + k*log(F) + k*log(1 - F) + logf
          }
          return(logL)
        }
        
        Lmedian <- likelihood_median_exact(M=medianobs, n, mean, sd)
      } else {
        Lmedian <- 0
      }
      
      if (!is.null(quantilesobs)) {
        
        likelihood_quantile <- function(Q, n, p = 0.5, mean = 0, sd = 1) {
          # choose k using the "floor((n+1)p)" convention
          k <- floor((n + 1) * p)
          k <- ifelse(k>1, 1, k)
          k <- ifelse(k>n, n, k)
          
          # Normal pdf and cdf at M
          logfM <- dnorm(Q, mean = mean, sd = sd, log = TRUE)
          FM <- pnorm(Q, mean = mean, sd = sd)
          
          # Exact PDF of the k-th order statistic
          logcoef <- lchoose(n - 1, k - 1)
          logL_exact <- logcoef +  (k - 1)*log(FM) + (n - k)*log((1 - FM)) + logfM
          
          # Asymptotic approx (large n): variance = p(1-p) / (n * f(q_p)^2)
          if (FALSE) {
            zp <- qnorm(p)
            q_p <- mean + sd * zp
            f_qp <- dnorm(q_p, mean = mean, sd = sd)
            var_qp <- (p * (1 - p)) / (n * (f_qp^2))
            sd_qp <- sqrt(var_qp)
            L_approx <- dnorm(Q, mean = q_p, sd = sd_qp)
          }
          return(logL_exact)
        }
        
        Lquantiles <- 0
        vp <- as.numeric(substr(names(quantilesobs), 2, nchar(names(quantilesobs))))
        vQ <- unname(quantilesobs)
        Lquantiles <- Lquantiles + sum(likelihood_quantile(Q=vQ, n, p = vp, mean = mean, sd = sd))
        
        
      } else {
        Lquantiles <- 0
      }
      
    }
    
    if (D_ec == "lnorm**") {
      
      meanlog <- x["meanlog"]
      sdlog <- x["sdlog"]
      
      if (!is.null(meanobs)) {
        Lmean <- dnorm(x=meanobs, mean = exp(meanlog+sdlog^2/2), sd = sqrt((exp(sdlog^2)-1)*exp(2*meanlog+sdlog^2)) / sqrt(n), log = TRUE)
      } else {
        Lmean <- 0
      }
      
      if (!is.null(medianobs) | !is.null(quantilesobs)) {
        lognormal_quantile_likelihood <- function(Q, n, p=0.5, meanlog, sdlog) {
          # Determine the order statistic index for the p-th quantile
          r <- ceiling(p * n)   # or floor(p * n) depending on convention
          
          # PDF and CDF of lognormal
          f_Q <- dlnorm(Q, meanlog = meanlog, sdlog = sdlog)
          F_Q <- plnorm(Q, meanlog = meanlog, sdlog = sdlog)
          # Numerical safety: avoid log(0)
          F_Q <- pmin(pmax(F_Q, 1e-12), 1 - 1e-12)
          # Order statistic likelihood
          logL <- lchoose(n - 1, r - 1) +
            log(f_Q) +
            (r - 1) * log(F_Q) +
            (n - r) * log(1 - F_Q)
          
          return(logL)
        }
      }
      
      if (!is.null(medianobs)) {
        Lmedian <- lognormal_quantile_likelihood(Q=medianobs, n, p=0.5, meanlog, sdlog)
      } else {
        Lmedian <- 0
      }
      
      if (!is.null(quantilesobs)) {
        Lquantiles <- 0
        vp <- as.numeric(substr(names(quantilesobs), 2, nchar(names(quantilesobs))))
        vQ <- unname(quantilesobs)
        Lquantiles <- Lquantiles + lognormal_quantile_likelihood(Q=vQ, n, p = vp, meanlog, sdlog)

      } else {
        Lquantiles <- 0
      }
      
    }
    
    if (is.null(Lmean) | is.null(Lmedian)) {
      if (((!is.null(meanobs)) | (!is.null(medianobs))) & (is.null(kkk))) {
        rD <- get(paste0("r", D))
        dsit <- sapply(1:rep, FUN=function(i) {k <- do.call(rD, args = as.list(c(n=n, x))); return(c(max(-k), mean(k), median(k), max(k)))})
        kk <- matrix(dsit, ncol=4, byrow = TRUE)
        kkk <- apply(kk, MARGIN = 2, FUN=function(x) return(c(mean(x), sd(x), median(x))))
      }
      
      if (!is.null(meanobs)) {
        Lmean <- dnorm(x=meanobs, mean = kkk[1, 2], sd = kkk[2, 2]/ sqrt(n), log = TRUE)
      } else {
        Lmean <- 0
      }
      if (!is.null(medianobs)) {
        Lmedian <- dnorm(x=medianobs, mean = kkk[1, 3], sd = kkk[2, 3]/ sqrt(n), log = TRUE)
      } else {
        Lmedian <- 0
      }
    }
    
    L <- - ( Lmin + Lmax + Lmean + Lmedian + Lquantiles)
    
    if (is.infinite(L)) {
      stop("Problem during likelihood estimation")
      print(d(x))
      print(n)
      print(D_ec)
      print(minobs)
      print(maxobs)
    }
    # print(L)
    return(L)
  }
  
  # https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution
  # n <- 4; observed.Minimum <- 43; observed.Maximum <- 57
  # fitted.parameters <- c(mean=50, sd=4); fixed.parameters <- NULL
  
  if (!silent) cat(paste0("The data are supposed to be generated from a ", gsub("\\*", "", D), " distribution.\n"))
  
  # if (D == "lnorm**") {
  #   warning("The lnorm** model is not available; it has been changed to lnorm*.")
  #   D <- "lnorm*"
  # }
  
  if (is.null(priors) | is.character(priors)) {
    if ((D == "norm") | (D == "norm*") | (D == "norm**")) {
      if (priors == "dunif") {
        if (length(fitted.parameters) == 2) {
          priors <- setPriors(par=fitted.parameters, 
                              se=c(sd=1, mean=1), 
                              density = "dunif", 
                              rules = rbind(data.frame(Name="sd", Min=1E-6, Max=10*unname(fitted.parameters["sd"])), 
                                            data.frame(Name="mean", Min=-10*unname(fitted.parameters["mean"]), 
                                                       Max=10*unname(fitted.parameters["mean"]))))
          priors$SDProp <- unname(fitted.parameters)/10
        } else {
          priors <- setPriors(par=fitted.parameters, 
                              se=c(sd=1), 
                              density = "dunif", 
                              rules = data.frame(Name="sd", Min=1E-6, Max=10*unname(fitted.parameters)))
          priors$SDProp <- unname(fitted.parameters)/10
        }
      } else {
        if (length(fitted.parameters) == 2) {
          priors <- setPriors(par=fitted.parameters, 
                              se=c(sd=1, mean=1), 
                              density = "dnorm", 
                              rules = rbind(data.frame(Name="sd", Min=1E-6, Max=10*unname(fitted.parameters["sd"])), 
                                            data.frame(Name="mean", Min=-10*unname(fitted.parameters["mean"]), 
                                                       Max=10*unname(fitted.parameters["mean"]))))
          priors$SDProp <- unname(fitted.parameters)/10
        } else {
          priors <- setPriors(par=fitted.parameters, 
                              se=c(sd=1), 
                              density = "dnorm", 
                              rules = data.frame(Name="sd", Min=1E-6, Max=10*unname(fitted.parameters)))
          priors$SDProp <- unname(fitted.parameters)/10
        }
      }
    } else {
      if ((D == "lnorm") | (D == "lnorm*") | (D == "lnorm**")) {
        
        if (priors == "dunif") {
          if (length(fitted.parameters) == 2) {
            priors <- setPriors(par=fitted.parameters, 
                                se=c(sd=1, mean=1), 
                                density = "dunif", 
                                rules = rbind(data.frame(Name="sdlog", Min=1E-6, Max=10*unname(fitted.parameters["sdlog"])), 
                                              data.frame(Name="meanlog", Min=-10*unname(fitted.parameters["meanlog"]), 
                                                         Max=10*unname(fitted.parameters["meanlog"]))))
            priors$SDProp <- unname(fitted.parameters)/10
          } else {
            priors <- setPriors(par=fitted.parameters, 
                                se=c(sd=1), 
                                density = "dunif", 
                                rules = data.frame(Name="sdlog", Min=1E-6, Max=10*unname(fitted.parameters)))
            priors$SDProp <- unname(fitted.parameters)/10
          }
        } else {
          if (length(fitted.parameters) == 2) {
            priors <- setPriors(par=fitted.parameters, 
                                se=c(sd=1, mean=1), 
                                density = "dnorm", 
                                rules = rbind(data.frame(Name="sdlog", Min=1E-6, Max=10*unname(fitted.parameters["sdlog"])), 
                                              data.frame(Name="meanlog", Min=-10*unname(fitted.parameters["meanlog"]), 
                                                         Max=10*unname(fitted.parameters["meanlog"]))))
            priors$SDProp <- unname(fitted.parameters)/10
          } else {
            priors <- setPriors(par=fitted.parameters, 
                                se=c(sd=1), 
                                density = "dnorm", 
                                rules = data.frame(Name="sd", Min=1E-6, Max=10*unname(fitted.parameters)))
            priors$SDProp <- unname(fitted.parameters)/10
          }
        }
        
        
      } else {
        
        stop("Automatic priors is available only for lnorm and norm distributions.")
        
      }
      
    }
  }
  
  # try <- L_min_max(x=fitted.parameters                  , 
  #                  n=n                                  ,
  #                  fixed.parameters=fixed.parameters    , 
  #                  minobs=observed.Minimum              , 
  #                  maxobs=observed.Maximum              , 
  #                  medianobs=observed.Median            , 
  #                  meanobs=observed.Mean                , 
  #                  quantiles=observed.Quantiles         ,
  #                  rep=replicates                       , 
  #                  D=D                                  )
  
  out_mcmc <- MHalgoGen(likelihood = L_min_max                                    , 
                        parameters = priors                                       , 
                        fixed.parameters=fixed.parameters                         , 
                        meanobs=observed.Mean                                     , 
                        n=n                                                       , 
                        medianobs=observed.Median                                 ,
                        minobs=observed.Minimum                                   , 
                        maxobs=observed.Maximum                                   ,
                        quantiles=observed.Quantiles                              ,
                        D=D                                                       ,
                        n.iter=n.iter                                             , 
                        n.chains = n.chains                                       , 
                        n.adapt = n.adapt                                         , 
                        thin=thin                                                 , 
                        trace=trace                                               , 
                        rep=replicates                                            , 
                        adaptive=adaptive                                         )
  
  print(as.parameters(out_mcmc, index="quantile"))
  
  return(invisible(out_mcmc))
  
}

.dGEV <- function (x, par = NULL, location = 0, scale = 1, shape = 0, log=FALSE, sum=FALSE, silent=TRUE){
  if (!is.null(par)) {
    if (!is.na(par["location"])) location <- par["location"]
    if (!is.na(par["scale"])) scale <- par["scale"]
    if (!is.na(par["shape"])) shape <- par["shape"]
  }
  
  scale <- abs(scale)
  
  if (shape == 0) {
    # tx <- exp(-(x-location)/scale)
    # y <- (1/scale) * tx^(shape+1) * exp(-tx)
    # 
    logy <- - log(scale) - (x-location)/scale - exp(-(x-location)/scale)
    y <- exp(logy)
    if (any(is.infinite(logy))) {
      y <- ifelse(is.infinite(logy), 1E-99, y)
      logy <- log(y)
    }
  } else {
    tx <- (1+shape*(x-location)/scale)^(-1/shape)
    y <- (1/scale) * tx^(shape+1) * exp(-tx)
    logy <- log(y)
  }
  
  if (shape < 0) {
    y <- ifelse(x > (location + scale / abs(shape)), 1E-99, y)
    logy <- log(y)
  }
  if (shape > 0) {
    y <- ifelse(x < (location - scale/shape), 1E-99, y)
    logy <- log(y)
  }
  
  if (log) {
    y <- logy
    
    if (sum) {
      y <- unname(-sum(y))
      y <- ifelse(is.infinite(y), 1E99, y)
    }
  }
  y
  
}

# .dGEV_global <- function (x, par = NULL, location = 0, scale = 1, shape = 0, log=FALSE, sum=FALSE, silent=TRUE) 
# {
#   if (!is.null(par)) {
#     if (!is.na(par["location"])) location <- par["location"]
#     if (!is.na(par["scale"])) scale <- par["scale"]
#     if (!is.na(par["shape"])) shape <- par["shape"]
#   }
#   
#   scale <- abs(scale)
#   
#   if (!silent) {
#     print(c(location, scale, shape))
#   }
#   
#   # scale <- abs(scale)
#   
#   names.x <- names(x)
#   arg.mat <- suppressWarnings(cbind(x = as.vector(x), location = as.vector(location), 
#                                     scale = as.vector(scale), shape = as.vector(shape)))
#   na.index <- apply(arg.mat, MARGIN = 1, FUN = function(x) any(is.na(x)))
#   if (all(na.index)) {
#     y <- rep(NA, nrow(arg.mat))
#   } else {
#     y <- numeric(nrow(arg.mat))
#     y[na.index] <- NA
#     y.no.na <- y[!na.index]
#     for (i in c("x", "location", "scale", "shape")) assign(i, 
#                                                            arg.mat[!na.index, i])
#     if (any(scale < .Machine$double.eps)) 
#       stop("All values of 'scale' must be positive.")
#     z <- (x - location)/scale
#     shape.eq.0 <- shape == 0
#     if (any(shape.eq.0)) 
#       y.no.na[shape.eq.0] <- (exp(-exp(-z[shape.eq.0])) * 
#                                 exp(-z[shape.eq.0]))/scale[shape.eq.0]
#     shape.gt.0 <- shape > 0
#     if (any(shape.gt.0)) {
#       z.out <- z >= 1/shape
#       y.no.na[shape.gt.0 & z.out] <- 0
#       index <- shape.gt.0 & !z.out
#       if (any(index)) 
#         y.no.na[index] <- exp(-(1 - shape[index] * z[index])^(1/shape[index])) * 
#         (1/scale[index]) * ((1 - shape[index] * z[index])^((1/shape[index]) - 
#                                                              1))
#     }
#     shape.lt.0 <- shape < 0
#     if (any(shape.lt.0)) {
#       z.out <- z <= 1/shape
#       y.no.na[shape.lt.0 & z.out] <- 0
#       index <- shape.lt.0 & !z.out
#       if (any(index)) 
#         y.no.na[index] <- exp(-(1 - shape[index] * z[index])^(1/shape[index])) * 
#         (1/scale[index]) * ((1 - shape[index] * z[index])^((1/shape[index]) - 
#                                                              1))
#     }
#     y[!na.index] <- y.no.na
#   }
#   if (!is.null(names.x)) 
#     names(y) <- rep(names.x, length = length(x))
#   else names(y) <- NULL
#   
#   if (log) y <- log(y)
#   if (sum & log) {
#     y <- unname(-sum(y))
#     y <- ifelse(is.infinite(y), 1E9, y)
#   }
#   y
# }




