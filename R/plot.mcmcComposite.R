#' plot.mcmcComposite plots the result of a MCMC search
#' @title Plot the result of a mcmcComposite object
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return None
#' @param x A mcmcComposite object
#' @param chains The chains to use
#' @param parameters Name of parameters or "all"
#' @param transform Function to be used to transform the variable
#' @param legend If FALSE, the legend is not shown; see description
#' @param show.prior Should the prior be shown?
#' @param show.posterior.density Should the posterior density be shown?
#' @param scale.prior If TRUE, the prior is scaled at the same size as posterior
#' @param bty Design of box for Markov Chain plot
#' @param col.prior Color for prior curve 
#' @param lty.prior Type of line for prior curve 
#' @param lwd.prior Width of line for prior curve
#' @param col.posterior Color for posterior histogram 
#' @param lty.posterior Type of line for posterior histogram
#' @param lwd.posterior Width of line for posterior histogram
#' @param ylab y-label for posterior
#' @param ylab.prior y-label for prior
#' @param show.yaxis.prior Should the y-axis for prior be shown
#' @param las las parameter (orientation of y-axis graduation)
#' @param ... Graphical parameters to be sent to hist() or plot()
#' @param what can be Posterior, MarkovChain or LnL
#' @family mcmcComposite functions
#' @description Plot the results within a mcmcComposite object.\cr
#' If scale.prior is TRUE, another scale is shown at right.\cr
#' legend can take these values: \cr
#' FALSE, TRUE, topleft, topright, bottomleft, bottomright, c(x=, y=)
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
#'                               Prior1=c(10, 0.5), Prior2=c(2, 0.5), 
#'                               SDProp=c(1, 1), 
#'                               Min=c(-3, 0), Max=c(100, 10), 
#'                               Init=c(10, 2), stringsAsFactors = FALSE, 
#'                               row.names=c('mean', 'sd'))
#'                               
#' mcmc_run <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, data=x, 
#'                       adaptive = TRUE,
#'                       likelihood=dnormx, n.chains=4, 
#'                       n.adapt=100, thin=1, trace=1)
#' plot(mcmc_run, xlim=c(0, 20), parameters="mean")
#' plot(mcmc_run, xlim=c(0, 10), parameters="sd")
#' plot(x=mcmc_run, what="MarkovChain", ylim=c(0, 15), parameters="mean", 
#'      col.posterior = colorRampPalette(c("blue", "grey"), alpha = 0.001)
#'                      (mcmc_run$parametersMCMC$n.chains))
#' plot(mcmc_run, what="MarkovChain", ylim=c(0, 10), parameters="sd", 
#'      col.posterior = colorRampPalette(c("blue", "grey"), alpha = 0.001)
#'                      (mcmc_run$parametersMCMC$n.chains))
#' plot(mcmc_run, what="LnL", 
#'      col.posterior = colorRampPalette(c("blue", "grey"), alpha = 0.001)
#'                      (mcmc_run$parametersMCMC$n.chains))
#' mcmcforcoda <- as.mcmc(mcmc_run)
#' heidel.diag(mcmcforcoda)
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
#' mcmc_run2 <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, data=x, 
#' adaptive = TRUE,
#' likelihood=dnormx, n.chains=1, n.adapt=1, thin=1, trace=1)
#' mcmc_run3 <- merge(mcmc_run, mcmc_run2)
#' ####### no adaptation, n.adapt must be 0
#' parameters_mcmc[,"Init"] <- c(mean(x), sd(x))
#' mcmc_run3 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, data=x, 
#' adaptive = TRUE,
#' likelihood=dnormx, n.chains=1, n.adapt=0, thin=1, trace=1)
#' 
#' #########################
#' ## Example with transform
#' #########################
#' 
#' x.1<-rnorm(6000, 2.4, 0.6)
#' x.2<-rlnorm(10000, 1.3,0.1)
#' 
#' X<-c(x.1, x.2)
#' hist(X,100,freq=FALSE, ylim=c(0,1.5))
#' Lnormlnorm <- function(par, val) {
#'   p <- invlogit(par["p"])
#'   return(-sum(log(p*dnorm(val, par["m1"], abs(par["s1"]), log = FALSE)+
#'                     (1-p)*dlnorm(val, par["m2"], abs(par["s2"]), log = FALSE))))
#' }
#' # Mean 1
#' m1=2.3; s1=0.5
#' # Mean 2
#' m2=1.3; s2=0.1
#' # proportion of category 1 - logit transform
#' p=0
#' 
#' par<-c(m1=m1, s1=s1, m2=m2, s2=s2, p=p)
#' 
#' result2<-optim(par, Lnormlnorm, method="BFGS", val=X, 
#'               hessian=FALSE, control=list(trace=1))
#'               
#' lines(seq(from=0, to=5, length=100), 
#' dnorm(seq(from=0, to=5, length=100), 
#'       result2$par["m1"], abs(result2$par["s1"])), col="red")
#' 
#' lines(seq(from=0, to=5, length=100), 
#'       dlnorm(seq(from=0, to=5, length=100), 
#'              result2$par["m2"], abs(result2$par["s2"])), col="green")
#' 
#' p <- invlogit(result2$par["p"])
#' 
#' paste("Proportion of Gaussian data",  p)
#' 
#' lines(seq(from=0, to=5, length=100), 
#'       p*dnorm(seq(from=0, to=5, length=100), 
#'               result2$par["m1"], result2$par["s1"])+
#'         (1-p)*dlnorm(seq(from=0, to=5, length=100), 
#'                      result2$par["m2"], result2$par["s2"]), col="blue")               
#'
#' parameters_mcmc <- data.frame(Density=c('dunif', 'dunif', 'dunif', 'dunif', 'dunif'), 
#'                                         Prior1=c(0, 0.001, 0, 0.001, -3), 
#'                                         Prior2=c(10, 10, 10, 10, 3), 
#'                                         SDProp=c(1, 1, 1, 1, 1), 
#'                                         Min=c(0, 0.001, 0, 0.001, -3), 
#'                                         Max=c(10, 10, 10, 10, 3), 
#'                                         Init=result2$par, stringsAsFactors = FALSE, 
#'                                         row.names=c('m1', 's1', 'm2', 's2', 'p'))
#'                                         
#' mcmc_run <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, val=X, 
#'                       parameters_name = "par",
#'                       adaptive = TRUE,
#'                       likelihood=Lnormlnorm, n.chains=1, 
#'                       n.adapt=100, thin=1, trace=100)
#' plot(mcmc_run, parameters="m1", breaks=seq(from=0, to =10, by=0.1), 
#'      legend=c(x=6, y=3.10))
#' plot(mcmc_run, parameters="p", transform=invlogit, xlim=c(0,1), 
#'      breaks=seq(from=0, to=1, by=0.01), legend=c(x=0.6, y=10))
#' plot(mcmc_run, parameters="p", xlim=c(-3,3), 
#'      breaks=seq(from=-3, to =3, by=0.05), legend=c(x=1, y= 3.10))
#' plot(mcmc_run, parameters="p", what="MarkovChain")
#' plot(mcmc_run, parameters="p", transform=invlogit, what="MarkovChain", 
#'      ylim=c(0, 1))
#' plot(mcmc_run, what="lnl", col.posterior = "black")
#'      
#' parameters_mcmc <- data.frame(Density=c('dunif', 'dunif', 'dunif', 'dunif', 'dnorm'), 
#'                                         Prior1=c(0, 0.001, 0, 0.001, 0.5), 
#'                                         Prior2=c(10, 10, 10, 10, 1), 
#'                                         SDProp=c(1, 1, 1, 1, 1), 
#'                                         Min=c(0, 0.001, 0, 0.001, -3), 
#'                                         Max=c(10, 10, 10, 10, 3), 
#'                                         Init=result2$par, stringsAsFactors = FALSE, 
#'                                         row.names=c('m1', 's1', 'm2', 's2', 'p'))
#'                                         
#' mcmc_run_pnorm <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, val=X, 
#'                       parameters_name = "par",
#'                       adaptive = TRUE,
#'                       likelihood=Lnormlnorm, n.chains=1, 
#'                       n.adapt=100, thin=1, trace=100)
#' plot(mcmc_run_pnorm, parameters="m1", breaks=seq(from=0, to =10, by=0.1), 
#'      legend=c(x=6, y=0.10))
#' plot(mcmc_run_pnorm, parameters="p", transform=invlogit, xlim=c(0,1), 
#'      breaks=seq(from=0, to=1, by=0.01), legend=c(x=0.6, y=0.10))
#' plot(x=mcmc_run_pnorm, parameters="p", xlim=c(-3,3), 
#'      breaks=seq(from=-3, to =3, by=0.05), legend=c(x=1, y= 0.10))
#'      
#'      
#' # Note that it is more logic to use beta distribution for p as a  
#' # proportion. However p value must be checked to be used in optim
#' # The use of logit transform can be a problem because it can stuck 
#' # the p value to 1 or 0 during fit.
#' 
#' Lnormlnorm <- function(par, val) {
#'   p <- par["p"]
#'   return(-sum(log(p*dnorm(val, par["m1"], abs(par["s1"]), log = FALSE)+
#'                     (1-p)*dlnorm(val, par["m2"], abs(par["s2"]), log = FALSE))))
#' }
#' 
#' # Example of beta distribution
#' 
#' # Mean is alpha/(alpha+beta)
#' # Variance is (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
#' alpha = 5
#' beta = 9
#' plot(x = seq(0.0001, 1, by = .0001), 
#'     y = dbeta(seq(0.0001, 1, by = .0001), alpha, beta),
#'     type = "l", ylab="Density", xlab="p", bty="n")
#' points(x=alpha/(alpha+beta), y=0, pch=4)
#' segments(x0=alpha/(alpha+beta)-sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1))), 
#'         x1=alpha/(alpha+beta)+sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1))),
#'         y0=0, y1=0)
#'
#' # Use of optim with L-BFGS-B to limit p between 0 and 1 and s > 0
#' 
#' # Mean 1
#' m1=2.3; s1=0.5
#' # Mean 2
#' m2=1.3; s2=0.1
#' # proportion of category 1 - logit transform
#' p=0.5
#' 
#' par <- c(m1=m1, s1=s1, m2=m2, s2=s2, p=p)
#' 
#' result2 <- optim(par, Lnormlnorm, method="L-BFGS-B", val=X, 
#'               lower = c(-Inf, 0, -Inf, 0, 0), 
#'               upper = c(Inf, Inf, Inf, Inf, 1),
#'               hessian=FALSE, control=list(trace=1))
#' 
#' parameters_mcmc <- data.frame(Density=c('dunif', 'dunif', 'dunif', 'dunif', 'dbeta'), 
#'                                         Prior1=c(0, 0.001, 0, 0.001, 5), 
#'                                         Prior2=c(10, 10, 10, 10, 9), 
#'                                         SDProp=c(1, 1, 1, 1, 1), 
#'                                         Min=c(0, 0.001, 0, 0.001, 0), 
#'                                         Max=c(10, 10, 10, 10, 1), 
#'                                         Init=c('m1' = 2.4, 
#'                                                's1' = 0.6, 
#'                                                'm2' = 1.3, 
#'                                                's2' = 0.1, 
#'                                                'p' = 0.5), stringsAsFactors = FALSE, 
#'                                         row.names=c('m1', 's1', 'm2', 's2', 'p'))
#'                                         
#' mcmc_run_pbeta <- MHalgoGen(n.iter=50000, parameters=parameters_mcmc, val=X, 
#'                       parameters_name = "par",
#'                       adaptive = TRUE,
#'                       likelihood=Lnormlnorm, n.chains=1, 
#'                       n.adapt=100, thin=1, trace=100)
#' plot(mcmc_run_pbeta, parameters="m1", breaks=seq(from=0, to =10, by=0.1), 
#'      legend=c(x=6, y=2.10))
#' plot(mcmc_run_pbeta, parameters="p", xlim=c(0,1), 
#'      breaks=seq(from=0, to=1, by=0.01), legend=c(x=0.6, y=6))
#' 
#' 
#' }
#' @method plot mcmcComposite
#' @export

plot.mcmcComposite <- function(x, 
                               ... , 
                               chains="all"                                                      , 
                               parameters=1                                                      , 
                               transform=NULL                                                    , 
                               scale.prior=TRUE                                                  , 
                               legend="topright"                                                 , 
                               ylab="Posterior density"                                          ,
                               las = 1                                                           ,
                               bty = "n"                                                         ,
                               show.prior=TRUE                                                   , 
                               show.posterior.density=TRUE                                       ,
                               col.prior = "red"                                                 , 
                               lty.prior = 1                                                     , 
                               lwd.prior = 1                                                     , 
                               what = "Posterior"                                                ,
                               col.posterior = colorRampPalette(c("blue", "grey"), alpha=0.001)  , 
                               lty.posterior = 1                                                 , 
                               lwd.posterior = 1                                                 ,
                               show.yaxis.prior=TRUE                                              ,
                               ylab.prior="Prior density"                                        ) {
  
  # x<-NULL;transform=NULL;chains="all";parameters=1;scale.prior=TRUE;legend=TRUE;ylab="Posterior density";las=1;show.prior=TRUE;col.prior="red";lty.prior=1;lwd.prior=1;col.posterior="white";lty.posterior=1;lwd.posterior=1;ylab.prior="Prior density"; what="posterior"
  # col.posterior <- colorRampPalette(c("blue", "grey"))(x$parametersMCMC$n.chains)
  # print(col.posterior)
  
  if (all(chains == "all")) chains <- 1:x$parametersMCMC$n.chains
  
  what <- match.arg(tolower(what), choices=c("posterior", "markovchain", "lnl"))
  chains <- chains[chains <= x$parametersMCMC$n.chains]
  
  if (is.null(transform)) {
    transformx <- function(x) {x}
  } else {
    transformx <- transform
  }
  
  pre.mar <- par("mar")
  
  # chains=1; parameters=1; scale.prior=FALSE; legend=TRUE
  
  Parameters_L <- x[["parametersMCMC"]]
  # Parameters avec une majuscule, c'est le tableau pMCMC
  # avec une minuscule, c'est le paramètres à afficher !
  Parameters <- Parameters_L[["parameters"]]
  n.iter <- Parameters_L[["n.iter"]]
  n.chains <- Parameters_L[["n.chains"]]
  n.adapt <- Parameters_L[["n.adapt"]]
  thin <- Parameters_L[["thin"]]
  
  possible <- rownames(Parameters)
  NbTS <- length(possible)
  
  if (parameters[[1]]=="all") {
    parameters <- possible
  }
  
  tpt <- list(...) # tpt <- list()
  
  if (what == "lnl") {
    
    if (is.function(col.posterior)) col.posterior <- col.posterior(length(chains))
    
    L <- modifyList(modifyList(list(xlab="Iterations", 
                                    ylab="Ln L", 
                                    las=las, main="", 
                                    col=col.posterior[1],
                                    lty=lty.posterior[1], 
                                    lwd=lwd.posterior[1]), modifyList(list(x=1:length(x$resultLnL[[1]]), y=x$resultLnL[[1]], 
                                                                           bty=bty, type="l"), tpt)), list(type="n"))
    do.call(plot, L)
    
    col.posterior <- rep(col.posterior, length(chains))[1:length(chains)]
    lty.posterior <- rep(lty.posterior, length(chains))[1:length(chains)]
    lwd.posterior <- rep(lwd.posterior, length(chains))[1:length(chains)]
    
    for (chain in chains) {
      vals <- x[["resultLnL"]][[chain]]
      lines(x=1:length(vals), y=vals, col=col.posterior[chain], lty=lty.posterior[chain], lwd=lwd.posterior[chain])
    }
    
    
  }
  
  nitercorrige <- floor(n.iter/thin)
  if (is.numeric(parameters)) parameters <- possible[parameters]
  
  
  if (what == "posterior") {
    if (show.prior) {
      if (scale.prior) {
        pre.mar.2 <- pre.mar
        if (pre.mar.2[4] < 5.1)  pre.mar.2[4] <- 5.1
        if (pre.mar.2[2] < 5.1)  pre.mar.2[2] <- 5.1
        par(mar=pre.mar.2)
      }
    }
    
    if (is.function(col.posterior)) col.posterior <- col.posterior(1)
    
    for (variable in parameters) {
      
      vals <- NULL
      for (chain in chains)
        if (NbTS == 1) vals <- c(vals, x[["resultMCMC"]][[chain]]) else vals <- c(vals, x[["resultMCMC"]][[chain]][ ,variable])
        # c'est quoi ça ??? 6/10/2012
        # vals <- vals[(length(vals)-nitercorrige):length(vals)]
        
        rg <- range(vals)
        vals <- transformx(vals)
        
        
        if (Parameters[variable, "Density"]=="dunif") {
          # 22/8/2014
          # je suis en uniforme. Il faut que je limite les breaks
          # 9/3/2015 Si l'un est négatif ou nul, je dois faire un shift d'origine
          shift <- as.numeric(Parameters[variable, "Prior1"])
          mx <- as.numeric(Parameters[variable, "Prior2"])-shift
          mn <- as.numeric(Parameters[variable, "Prior1"])-shift
          mxlog <- 10^(floor(log10(mx))-1)
          # C'est quoi ???
          xl <- as.numeric(c(mn, mx)+shift)+c(-mxlog, mxlog)
          br1 <- mn
          br2 <- mx
          interval <- (mx-mn)/20
          decalage <- br1 %% interval
          
          br <- seq(from=decalage+shift, to=mx+shift, by=interval)
        } else {
          shift <- as.numeric(Parameters[variable, "Min"])
          mx <- as.numeric(Parameters[variable, "Max"])-shift
          mn <- as.numeric(Parameters[variable, "Min"])-shift
          mxlog <- 10^(floor(log10(mx))-1)
          # C'est quoi ???
          xl <- as.numeric(c(mn, mx)+shift)+c(-mxlog, mxlog)
          br1 <- mn
          br2 <- mx
          interval <- (mx-mn)/20
          decalage <- br1 %% interval
          
          br <- seq(from=decalage+shift, to=mx+shift, by=interval)
          
        }
        
        # tpt <- list(las=1, xlim=c(0,30), breaks=c(0, 1.00095, 2.0009, 3.00085, 4.0008, 5.00075, 6.0007, 7.00065, 8.0006, 9.00055, 10.0005, 11.00045, 12.0004, 13.00035, 14.0003, 15.00025, 16.0002, 17.00015, 18.0001, 19.00005, 20))
        
        L <- modifyList(list(ylab="", xlab=variable, 
                             las=las, main="", freq=FALSE, 
                             xlim=xl, breaks=br, 
                             col=col.posterior[1], 
                             lty=lty.posterior[1], 
                             lwd=lwd.posterior[1]), modifyList(list(x=vals), tpt)) 
        
        do.call(hist, L) 
        mtext(ylab, side=2, line=3)
        
        
        
        if (show.prior) {
          
          scl <- ScalePreviousPlot()
          sequence <- seq(from=scl$xlim[1], to=scl$xlim[2], length=200)
          
          p1 <- as.numeric(Parameters[variable, "Prior1"])
          p2 <- as.numeric(Parameters[variable, "Prior2"])
          y <- get(as.character(Parameters[variable, "Density"]))(sequence, p1, p2)
          yl <- c(0, max(y[is.finite(y)]))
          
          sequence <- transformx(sequence)
          
          
          if (Parameters[variable, "Density"] != "dunif") {
            
            
            if (scale.prior) {
              lines(x=sequence, y=(y/max(y[is.finite(y)]))*scl$ylim["end"], 
                    col=col.prior, lty=lty.prior, lwd=lwd.prior)
              if (show.yaxis.prior) {
                axis(side = 4, las=las)
                mtext(text=ylab.prior, side=4, line = 4)
              }
            } else {
              lines(x=sequence, y=y, 
                    col=col.prior, lty=lty.prior, lwd=lwd.prior)
            }
          } else {
            if (scale.prior) {
              # par(new=TRUE)
              # plot(c(scl$xlim["begin"], scl$xlim["end"]), c(0, max(y)), type="n", axes=FALSE, bty="n",
              #      xlab="", ylab="", main="")
              
              segments(x0=scl$xlim[1], x1=transformx(p1), y0=0, y1=0, col=col.prior, lty=lty.prior, lwd=lwd.prior)
              segments(x0=transformx(p1), x1=transformx(p1), y0=0, y1=scl$ylim[2], col=col.prior, lty=lty.prior, lwd=lwd.prior)
              segments(x0=transformx(p1), x1=transformx(p2), y0=scl$ylim[2], y1=scl$ylim[2], col=col.prior, lty=lty.prior, lwd=lwd.prior)
              segments(x0=transformx(p2), x1=transformx(p2), y0=scl$ylim[2], y1=0, col=col.prior, lty=lty.prior, lwd=lwd.prior)
              segments(x0=transformx(p2), x1=scl$xlim[2], y0=0, y1=0, col=col.prior, lty=lty.prior, lwd=lwd.prior)
              if (show.yaxis.prior) {
                axis(side = 4, las=las)
                mtext(text=ylab.prior, side=4, line = 4)
              }
              
            } else {
              segments(x0=scl$xlim[1], x1=transformx(p1), y0=0, y1=0, col=col.prior, lty=lty.prior, lwd=lwd.prior)
              segments(x0=transformx(p1), x1=transformx(p1), y0=0, y1=yl[2], col=col.prior, lty=lty.prior, lwd=lwd.prior)
              segments(x0=transformx(p1), x1=transformx(p2), y0=yl[2], y1=yl[2], col=col.prior, lty=lty.prior, lwd=lwd.prior)
              segments(x0=transformx(p2), x1=transformx(p2), y0=yl[2], y1=0, col=col.prior, lty=lty.prior, lwd=lwd.prior)
              segments(x0=transformx(p2), x1=scl$xlim[2], y0=0, y1=0, col=col.prior, lty=lty.prior, lwd=lwd.prior)
            }
          }
        }
    }
    
    
    if (legend[1] != FALSE) {
      if (length(legend) == 2) {
        if(show.prior) {
          legend(x=legend["x"], y=legend["y"], legend=c("Prior", "Posterior"), lty=c(lty.prior, lty.posterior), 
                 col=c(col.prior, col.posterior), bty="n", pch=c(NA, 0))
        } else {
          legend(x=legend["x"], y=legend["y"], legend=c("Posterior"), lty=c(lty.posterior), 
                 col=c(col.posterior), bty="n", pch=c(0))
        }
      } else {
        if (isTRUE(legend)) legend <- "topright"
        if(show.prior) {
          legend(x=legend, legend=c("Prior", "Posterior"), lty=c(lty.prior, lty.posterior), 
                 col=c(col.prior, col.posterior), bty="n", pch=c(NA, 0))
        } else {
          legend(x=legend, legend=c("Posterior"), lty=c(lty.posterior), 
                 col=c(col.posterior), bty="n", pch=c(0))
        }
      }
    }
  }
  
  if (what == "markovchain") {
    
    if (is.function(col.posterior)) col.posterior <- col.posterior(length(chains))
    
    for (variable in parameters) {
      
      # 22/8/2014
      # je suis en uniforme. Il faut que je limite les breaks
      # 9/3/2015 Si l'un est négatif ou nul, je dois faire un shift d'origine
      shift <- as.numeric(Parameters[variable, "Min"])
      mx <- as.numeric(Parameters[variable, "Max"])-shift
      mn <- as.numeric(Parameters[variable, "Min"])-shift
      mxlog <- 10^(floor(log10(mx))-1)
      # C'est quoi ???
      xl <- as.numeric(c(mn, mx)+shift)+c(-mxlog, mxlog)
      xl <- transformx(xl)
      if (NbTS == 1) vals <- x$resultMCMC[[1]] else vals <- x$resultMCMC[[1]][ ,variable]
      maxx <- length(vals)
      if (show.prior) maxx <- maxx + 0.1*maxx
      vals <- transformx(vals)
      
      L <- modifyList(modifyList(list(xlab="Iterations", 
                                      ylab=variable, 
                                      las=las, main="", 
                                      ylim=xl, 
                                      col=col.posterior[1],
                                      lty=lty.posterior[1], 
                                      lwd=lwd.posterior[1]), modifyList(list(ylab=ylab, x=1:length(vals), y=rep(0, length(vals)), 
                                                                             xlim=c(1, maxx), bty=bty), tpt)), list(type="n"))
      L <- modifyList(L, tpt)
      
      do.call(plot, L)
      
      col.posterior <- rep(col.posterior, length(chains))[1:length(chains)]
      lty.posterior <- rep(lty.posterior, length(chains))[1:length(chains)]
      lwd.posterior <- rep(lwd.posterior, length(chains))[1:length(chains)]
      for (chain in chains) {
        if (NbTS == 1) vals <- x[["resultMCMC"]][[chain]] else vals <- x[["resultMCMC"]][[chain]][ ,variable]
        vals <- transformx(vals)
        lines(x=1:length(vals), y=vals, col=col.posterior[chain], lty=lty.posterior[chain], lwd=lwd.posterior[chain])
      }
      
      if (show.posterior.density) {
        segments(x0=length(vals), x1=length(vals), y0=ScalePreviousPlot()$ylim["begin"], 
                 y1=ScalePreviousPlot()$ylim["end"], col="black")
        for (chain in chains) {
          if (NbTS == 1) valsd <- x[["resultMCMC"]][[chain]] else valsd <- x[["resultMCMC"]][[chain]][ ,variable]
          valsd <- transformx(valsd)
        y <- density(valsd)
        lines(x=length(vals) + (y$y/max(y$y) * length(vals) *0.1) , y=y$x, 
              col=col.posterior[chain], lty=lty.posterior[chain], lwd=lwd.posterior[chain])
        }
      }
      
      if (show.prior) {
        scl <- ScalePreviousPlot()
        
        # sequence <- seq(from=scl$ylim[1], to=scl$ylim[2], length=200)
        
        sequence <- seq(from=as.numeric(Parameters[variable, "Min"]), to=as.numeric(Parameters[variable, "Max"]), length=200)
        
        p1 <- as.numeric(Parameters[variable, "Prior1"])
        p2 <- as.numeric(Parameters[variable, "Prior2"])
        y <- get(as.character(Parameters[variable, "Density"]))(sequence, p1, p2)
        
        # sequence <- transformx(sequence)
        lines(x=length(vals) + (y/max(y) * length(vals) *0.1) , y=transformx(sequence), 
              col=col.prior, lty=lty.prior, lwd=lwd.prior)
      }
      
    }
    
  }
  
  par(mar=pre.mar)
}
