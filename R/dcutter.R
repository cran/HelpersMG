#' dcutter returns the density of the cutter function
#' @title Distribution of the fitted distribution without cut. 
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return The density of the cutter function according to observations.
#' @param par Values for parameters of distribution 
#' @param observations The observations; see description.
#' @param distribution Can be gamma, normal, weibull, lognormal, or generalized.gamma.
#' @param n.mixture Number of distributions
#' @param debug If TRUE, show some information. If 2, show more information.
#' @param limits.lower Value for lower detection limit 
#' @param limits.upper Value for upper detection limit
#' @param log If TRUE, return the log likelihood
#' @description If observations must be a data.frame with 4 columns:\cr
#' observations: A column for the measurements;\cr
#' LDL: A column for the lower detection limit;\cr
#' UDL: A column for the upper detection limit;\cr
#' Cut: A column for the truncated of censored nature of the data.\cr
#' @family Distributions
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' par <- c('shape1' = 0.42265849507444225, 
#'          'scale1' = 14.139457094879594, 
#'          'shape2' = 1.667131542489706, 
#'          'scale2' = 0.10763344388223803, 
#'          'p1' = 0.12283307526788023)
#' obs <- data.frame(Observations=c(0.755, 1.013, 2.098, 6.265, 4.708, 0.078, 2.169, 0.403, 1.251, 
#'                                  0.008, 1.419, 1.078, 2.744, 81.534, 1.426, 13.486, 7.813, 0.165, 
#'                                  0.118, 0.864, 0.369, 7.159, 2.605, 1.579, 1.646, 0.484, 4.492, 
#'                                  0.139, 0.28, 0.154, 0.106, 0.104, 4.185, 0.735, 0.149, 0.183, 
#'                                  0.062, 8.246, 0.165, 0.121, 0.109, 0.092, 0.162, 0.108, 0.139, 
#'                                  0.141, 0.124, 0.124, 0.151, 0.141, 0.364, 0.295, 0.09, 0.135, 
#'                                  0.154, 0.218, 0.167, -Inf, 0.203, 0.228, 0.107, 0.162, 0.194, 
#'                                  0.322, 0.351, 0.17, 0.236, 0.176, 0.107, 0.12, 0.095, 0.27, 0.194, 
#'                                  0.125, 0.123, 0.085, 0.164, 0.106, 0.079, 0.162), 
#'                  LDL=0.001, UDL=NA, Cut="censored")
#' dcutter(par=par, observations=obs, distribution="gamma", 
#'         n.mixture=NULL, debug=FALSE, limits.lower=NULL, 
#'         limits.upper=NULL,log=FALSE)
#' dcutter(par=par, observations=obs, distribution="gamma", 
#'         n.mixture=NULL, debug=FALSE, limits.lower=NULL, 
#'         limits.upper=NULL, log=TRUE)
#'         
#' }
#' @export

dcutter <- function(par, observations=NULL, distribution="gamma", 
                    n.mixture=NULL, debug=FALSE, limits.lower=NULL, 
                    limits.upper=NULL, log = TRUE) {
  
  # par=NULL; observations=NULL; distribution="gamma"; n.mixture=NULL; debug=FALSE
  # limits.lower=NULL; limits.upper=NULL; log = TRUE
  
  # par=result$par
  # observations=result$observations
  # distribution=result$distribution
  # n.mixture=result$n.mixture
  # debug=2
  # limits.lower=NULL; limits.upper=NULL; log = TRUE
  
  
  getparcutter <- getFromNamespace(".getparcutter", ns="HelpersMG")
  debug <- as.numeric(debug)
  
  if (debug == 2) print(par)
  
  if (!is.null(limits.lower)) {
    if (any(limits.lower > par)) {
      if (debug >= 1) {
        np <- names(par[limits.lower > par])
        for (i in np) {
          cat(paste0("Parameter ", i," is lower than lower limit during search", "\n"))
          cat(paste0("The parameter was ", as.character(par[i]), " and "))
          cat(paste0("the lower limit was ", as.character(limits.lower[i]), 
                     " (difference was ", specify_decimal(limits.lower[i]-par[i]), ")\n"))
        }
      }
      return(-1E6)
    }
  }
  
  
  if (!is.null(limits.upper)) {
    if (any(limits.upper < par)) {
      if (debug >= 1) {
        np <- names(par[limits.upper < par])
        for (i in np) {
          cat(paste0("Parameter ", i," is higher than upper limit during search", "\n"))
          cat(paste0("The parameter was ", as.character(par[i]), " and "))
          cat(paste0("the upper limit was ", as.character(limits.upper[i]), 
                     " (difference was ", specify_decimal(-limits.upper[i]+par[i]), ")\n"))
        }
      }
      return(-1E6)
    }
  }
  
  # par <- c('shape1' = 1, 'scale1' = 3, 'shape2' = 4, 'scale2' = 6, 'p1' = 0.1)
  # par <- c('shape' = 1, 'scale' = 3)
  # observations <- data.frame(Observations=1:10, Cut=rep("censored", 10), LDL=rep(0.01, 10), UDL=rep(NA, 10))
  # distribution="gamma"
  # debug=FALSE
  
  # dcutter(par=result2$par, observations=result2$observations, distribution=result2$distribution, n.mixture=result2$n.mixture, debug=FALSE)
  
  
  
  # gsub(".+([0-9]).+", "\\1", names(par))
  
  if (is.null(n.mixture)) {
    np <- unique(as.numeric(gsub("[a-zA-Z]+([0-9]*)$", "\\1", names(par))))
    if (all(is.na(np))) np <- 0 else np <- max(np, na.rm=TRUE)
    
    if (np == 0) {
      names(par) <- paste0(names(par), "1")
      np <- 1
    }
  } else {
    np <- n.mixture
  }
  
  LDL <- observations[, "LDL"]
  UDL <- observations[, "UDL"]
  Cut <- observations[, "Cut"]
  obs <- observations[, "Observations"]
  
  obs_finite <- is.finite(obs)
  LDLX <- LDL[obs_finite]
  UDLX <- UDL[obs_finite]
  CutX <- Cut[obs_finite]
  obsX <- obs[obs_finite]
  
  parX <- as.list(par)
  pparX <- getparcutter(par, set=NULL)
  
  ddistr <- switch(EXPR = distribution, gamma=dgamma, lognormal=dlnorm, 
                   normal=dnorm, weibull=dweibull, 
                   generalized.gamma=dggamma)
  
  pdistr <- switch(EXPR = distribution, gamma=pgamma, lognormal=plnorm, 
                   normal=pnorm, weibull=pweibull, 
                   generalized.gamma=pggamma)
  
  
  # L <- Inf
  
  if (np == 1) {
    
    parX_mixture <- getparcutter(parX, set=1)
    
    L <- sum(do.call(ddistr, args = modifyList(list(x=obsX, log = TRUE), parX_mixture)))
    L <- L-sum(log(1 - ifelse(!is.na(LDLX), do.call(pdistr, args = modifyList(list(q=LDLX, lower.tail = TRUE, log.p = FALSE), parX_mixture)), 0) 
                   -ifelse(!is.na(UDLX), do.call(pdistr, args = modifyList(list(q=UDLX, lower.tail = FALSE, log.p = FALSE), parX_mixture)), 0)))
    
    if (debug >= 1) {
      cat(paste0("The log likelihood of quantifiable data is ", specify_decimal(L), "\n"))
    }
    
    
    if (any(observations$Cut == "censored")) {
      cs <- subset(observations, subset = (Cut == "censored"))
      fby1 <- factor(cs$LDL, exclude = "")
      fby2 <- factor(cs$UDL, exclude = "")
      
      d1 <- aggregate(x=cs, by=list(fby1, fby2), FUN=function(x) sum(is.finite(x)))
      dplus <- aggregate(x=cs, by=list(fby1, fby2), FUN=function(x) sum(x == +Inf))
      dmoins <- aggregate(x=cs, by=list(fby1, fby2), FUN=function(x) sum(x == -Inf))
      
      NA_character <- structure(1L, .Label = NA_character_, class = "factor")
      
      for (i in 1:nrow(d1)) {
        if ((!identical(d1[i, "Group.1"], NA_character)) | 
            (!identical(d1[i, "Group.2"], NA_character))) {
          # J'en ai au moins un qui a une valeur
          if ((!identical(d1[i, "Group.1"], NA_character)) & 
              (identical(d1[i, "Group.2"], NA_character))) {
            # J'ai une valeur sur LDL; donc je suis en left censored
            right <- d1[i, "Observations"]
            
            pright <- do.call(pdistr, modifyList(list(q=as.numeric(as.character(d1[i, "Group.1"])), lower.tail = FALSE, log.p = FALSE), parX_mixture))
            
            left <- dmoins[i, "Observations"]
            pleft <- 1-pright
            if (pleft <= 0) pleft <- 1E-6
            if (pleft >= 1) pleft <- 1-1E-6
            Lint <- dbinom(x=left, size=left+right, prob=pleft, log=TRUE)
            # if (!is.finite(Lint)) {
            #   message(paste0("L from binomial is infinite: ", d(par)))
            #   Lint <- 0
            # }
          } else {
            if ((identical(d1[i, "Group.1"], NA_character)) & 
                (!identical(d1[i, "Group.2"], NA_character))) {
              # J'ai une valeur sur UDL; donc je suis en right censored
              left <- d1[i, "Observations"]
              right <- dplus[i, "Observations"]
              pleft <- do.call(pdistr, modifyList(list(q=as.numeric(as.character(d1[i, "Group.2"])), lower.tail = TRUE, log.p = FALSE), parX_mixture))
              if (pleft <= 0) pleft <- 1E-6
              if (pleft >= 1) pleft <- 1-1E-6
              Lint <- dbinom(x=left, size=left+right, prob=pleft, log=TRUE)
              # if (!is.finite(Lint)) {
              #   message(paste0("L from binomial is infinite: ", d(par)))
              #   Lint <- 0
              # }
            } else {
              # if ((!identical(d1[i, "Group.1"], NA_character)) & 
              #     (!identical(d1[i, "Group.2"], NA_character))) {
              # J'ai une valeur sur LDL et UDL; donc je suis en left & right censored
              left <- dmoins[i, "Observations"]
              center <- d1[i, "Observations"]
              right <- dplus[i, "Observations"]
              
              pright <- do.call(pdistr, modifyList(list(q=as.numeric(as.character(d1[i, "Group.2"])), lower.tail = FALSE, log.p = FALSE), parX_mixture))
              pleft <- do.call(pdistr, modifyList(list(q=as.numeric(as.character(d1[i, "Group.1"])), lower.tail = TRUE, log.p = FALSE), parX_mixture))
              if (pright <= 0) pright <- 1E-6
              if (pright >= 1) pright <- 1-1E-6
              if (pleft <= 0) pleft <- 1E-6
              if (pleft >= 1) pleft <- 1-1E-6
              pcenter <- 1-pleft-pright
              if (pcenter <= 0) pcenter <- 1E-6
              if (pcenter >= 1) pcenter <- 1-1E-6
              Lint <- dmultinom(x=c(left, center, right), size = left+center+right, prob=c(pleft, pcenter, pright), log = TRUE)
            }
          }
          
          
          if (debug >= 1) {
            cat(paste0("The log likelihood of unquantifiable data LDL=", d1[i, "Group.1"]," and UDL=", d1[i, "Group.2"]," is ", specify_decimal(Lint), "\n"))
          }
          
          L <- L + Lint
          
        }
        
      }
      # fin du censored pour np = 1
    }
    # fin du np == 1
  } else {
    # Je fais chaque individu
    
    parX_mixture <- list()
    for (p in 1:np)
      parX_mixture <- c(parX_mixture, list(getparcutter(parX, set=p)))
    
    L <- sapply(1:length(obsX), FUN = function(n) {
      L_set <- sapply(1:np, FUN = function(p) {
        pr <- do.call(ddistr, args = modifyList(list(x=obsX[n], log = FALSE), parX_mixture[[p]]))
        pLDL <- 0
        if (!is.na(LDLX[n])) pLDL <- do.call(pdistr, args = modifyList(list(q=LDLX[n], lower.tail = TRUE, log.p = FALSE), parX_mixture[[p]]))
        pUDL <- 0
        if (!is.na(UDLX[n])) pUDL <- do.call(pdistr, args = modifyList(list(q=UDLX[n], lower.tail = FALSE, log.p = FALSE), parX_mixture[[p]]))
        
        prcond <- pr / (1 - pLDL - pUDL)
        return(pparX[p] * prcond)
      })
      return(log(sum(L_set)))
    })
    # Dans L j'ai log L
    L <- sum(L)
    
    if (debug >= 1) {
      cat(paste0("The log likelihood of quantifiable data is ", specify_decimal(L), "\n"))
    }
    
    # Là je dois traiter le censored
    
    if (any(observations$Cut == "censored")) {
      cs <- subset(observations, subset = (Cut == "censored"))
      fby1 <- factor(cs$LDL, exclude = "")
      fby2 <- factor(cs$UDL, exclude = "")
      
      d1 <- aggregate(x=cs, by=list(fby1, fby2), FUN=function(x) sum(is.finite(x)))
      dplus <- aggregate(x=cs, by=list(fby1, fby2), FUN=function(x) sum(x == +Inf))
      dmoins <- aggregate(x=cs, by=list(fby1, fby2), FUN=function(x) sum(x == -Inf))
      
      NA_character <- structure(1L, .Label = NA_character_, class = "factor")
      
      for (i in 1:nrow(d1)) {
        if ((!identical(d1[i, "Group.1"], NA_character)) | 
            (!identical(d1[i, "Group.2"], NA_character))) {
          # J'en ai au moins un qui a une valeur
          if ((!identical(d1[i, "Group.1"], NA_character)) & 
              (identical(d1[i, "Group.2"], NA_character))) {
            # J'ai une valeur sur LDL; donc je suis en left censored
            right <- d1[i, "Observations"]
            pright <- NULL
            for (p in 1:np) 
              pright <- c(pright, do.call(pdistr, modifyList(list(q=as.numeric(as.character(d1[i, "Group.1"])), lower.tail = FALSE, log.p = FALSE), parX_mixture[[p]] )))
            pright <- sum(pright * pparX)
            
            left <- dmoins[i, "Observations"]
            pleft <- 1-pright
            if (pleft <= 0) pleft <- 1E-6
            if (pleft >= 1) pleft <- 1-1E-6
            Lint <- dbinom(x=left, size=left+right, prob=pleft, log=TRUE)
            
          } else {
            if ((identical(d1[i, "Group.1"], NA_character)) & 
                (!identical(d1[i, "Group.2"], NA_character))) {
              # J'ai une valeur sur UDL; donc je suis en right censored
              left <- d1[i, "Observations"]
              right <- dplus[i, "Observations"]
              
              pleft <- NULL
              for (p in 1:np) 
                pleft <- c(pleft, do.call(pdistr, modifyList(list(q=as.numeric(as.character(d1[i, "Group.2"])), lower.tail = TRUE, log.p = FALSE), parX_mixture[[p]] )))
              pleft <- sum(pleft * pparX)
              
              
              # pright <- 1-pleft
              if (pleft <= 0) pleft <- 1E-6
              if (pleft >= 1) pleft <- 1-1E-6
              Lint <- dbinom(x=left, size=left+right, prob=pleft, log=TRUE)
              
            } else {
              left <- dmoins[i, "Observations"]
              center <- d1[i, "Observations"]
              right <- dplus[i, "Observations"]
              
              pright <- NULL
              for (p in 1:np) 
                pright <- c(pright, do.call(pdistr, modifyList(list(q=as.numeric(as.character(d1[i, "Group.2"])), lower.tail = FALSE, log.p = FALSE), parX_mixture[[p]] )))
              pright <- sum(pright * pparX)
              if (pright <= 0) pright <- 1E-6
              if (pright >= 1) pright <- 1-1E-6
              
              pleft <- NULL
              for (p in 1:np) 
                pleft <- c(pleft, do.call(pdistr, modifyList(list(q=as.numeric(as.character(d1[i, "Group.1"])), lower.tail = TRUE, log.p = FALSE), parX_mixture[[p]] )))
              pleft <- sum(pleft * pparX)
              
              if (pleft <= 0) pleft <- 1E-6
              if (pleft >= 1) pleft <- 1-1E-6
              pcenter <- 1-pleft-pright
              if (pcenter <= 0) pcenter <- 1E-6
              if (pcenter >= 1) pcenter <- 1-1E-6
              Lint <- dmultinom(x=c(left, center, right), size = left+center+right, prob=c(pleft, pcenter, pright), log = TRUE)
              
            }
          }
          
          
          
          if (debug >= 1) {
            cat(paste0("The log likelihood of unquantifiable data LDL=", d1[i, "Group.1"]," and UDL=", d1[i, "Group.2"]," is ", specify_decimal(Lint), "\n"))
          }
          
          L <- L + Lint
        }
      }
      
    }
  }
  
  
  if (as.numeric(debug) == 2) print(paste0("log Likelihood is ", specify_decimal(L)))
  if (is.na(L)) L <- -1E6
  if (is.infinite(L)) L <- -1E6
  # L est log L
  if (!log) L <- exp(L)
  
  return(L)
}




