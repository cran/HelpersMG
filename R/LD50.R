#' LD50 estimates the parameters that best describe LD50
#' @title Estimate the parameters that best describe LD50
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list with the LD50, Transitional Range of Doses and their SE
#' @param alive A vector with alive individuals at the end of experiment
#' @param dead A vector with dead individuals at the end of experiment
#' @param N A vector with total numbers of tested individuals
#' @param doses The doses
#' @param df A dataframe with at least two columns named alive, dead or N and doses columns
#' @param l The limit to define TRD (see Girondot, 1999)
#' @param parameters.initial Initial values for P, S or K search as a vector, ex. c(P=29, S=-0.3)
#' @param fixed.parameters Parameters that will not be changed during fit
#' @param SE Standard errors for parameters
#' @param equation Could be "logistic", "logit", "probit", Hill", "Richards", "Hulin", "flexit" or "Double-Richards"
#' @param replicates Number of replicates to estimate confidence intervals
#' @param range.CI The range of confidence interval for estimation, default=0.95
#' @param limit.low.TRD.minimum Minimum lower limit for TRD
#' @param limit.high.TRD.maximum Maximum higher limit for TRD
#' @param doses.plot Sequences of doses that will be used for plotting. If NULL, does not estimate them
#' @param print Do the results must be printed at screen? TRUE (default) or FALSE
#' @description Estimate the parameters that best describe LD50\cr
#' Logistic and logit models are the same but with different parametrization:\cr
#' logistic = 1/(1+exp((1/S)*(P-d)))\cr
#' logit = 1/(1+exp(P+d*S))\cr
#' See these publications for the description of equations:\cr
#' Girondot, M. 1999. Statistical description of temperature-dependent sex determination using maximum likelihood. Evolutionary Ecology Research, 1, 479-486.\cr
#' Godfrey, M.H., Delmas, V., Girondot, M., 2003. Assessment of patterns of temperature-dependent sex determination using maximum likelihood model selection. Ecoscience 10, 265-272.\cr
#' Hulin, V., Delmas, V., Girondot, M., Godfrey, M.H., Guillon, J.-M., 2009. Temperature-dependent sex determination and global change: are some species at greater risk? Oecologia 160, 493-506.\cr
#' The flexit equation is not still published :
#' \deqn{if dose < P then (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}{if dose < P then (1 + (2^K1 - 1) *  exp(4 * S1 * (P - x)))^(-1/K1)}
#' \deqn{if dose > P then 1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}{if dose > P then 1-((1 + (2^K2 - 1) * exp(4 * S2 * (x - P)))^(-1/K2)}
#' with:\cr
#'      \deqn{S1 = S/((4/K1)*(2^(-K1))^(1/K1+1)*(2^K1-1))}{S1 = S/((4/K1)*(2^(-K1))^(1/K1+1)*(2^K1-1))}
#'      \deqn{S2 = S/((4/K2)*(2^(-K2))^(1/K2+1)*(2^K2-1))}{S2 = S/((4/K2)*(2^(-K2))^(1/K2+1)*(2^K2-1))}
#' @family LD50 functions
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' data <- data.frame(Doses=c(80, 120, 150, 150, 180, 200),
#' Alive=c(10, 12, 8, 6, 2, 1),
#' Dead=c(0, 1, 5, 6, 9, 15))
#' LD50_logistic <- LD50(data, equation="logistic")
#' predict(LD50_logistic, doses=c(140, 170))
#' plot(LD50_logistic, xlim=c(0, 300), at=seq(from=0, to=300, by=50))
#' LD50_probit <- LD50(data, equation="probit")
#' predict(LD50_probit, doses=c(140, 170))
#' plot(LD50_probit)
#' LD50_logit <- LD50(data, equation="logit")
#' predict(LD50_logit, doses=c(140, 170))
#' plot(LD50_logit)
#' LD50_hill <- LD50(data, equation="hill")
#' predict(LD50_hill, doses=c(140, 170))
#' plot(LD50_hill)
#' LD50_Richards <- LD50(data, equation="Richards")
#' predict(LD50_Richards, doses=c(140, 170))
#' plot(LD50_Richards)
#' LD50_Hulin <- LD50(data, equation="Hulin")
#' predict(LD50_Hulin, doses=c(140, 170))
#' plot(LD50_Hulin)
#' LD50_DoubleRichards <- LD50(data, equation="Double-Richards")
#' predict(LD50_DoubleRichards, doses=c(140, 170))
#' plot(LD50_DoubleRichards)
#' LD50_flexit <- LD50(data, equation="flexit")
#' predict(LD50_flexit, doses=c(140, 170))
#' plot(LD50_flexit)
#' }
#' @export


LD50 <- function(df=NULL, alive=NULL, dead=NULL, N=NULL, 
                 doses=NULL, 
                 l=0.05, parameters.initial=NULL, 
                 fixed.parameters=NULL, SE=NULL,
                 equation="logistic", replicates=1000, range.CI=0.95, 
                 limit.low.TRD.minimum=5, limit.high.TRD.maximum=1000, print=TRUE, 
                 doses.plot=seq(from=0, to=1000, by=0.1)) {
  # df=NULL; alive=NULL; dead=NULL; N=NULL; doses=NULL; l=0.05; parameters.initial=NULL; fixed.parameters=NULL; SE=NULL;equation="logistic"; replicates=1000; range.CI=0.95; limit.low.TRD.minimum=5; limit.high.TRD.maximum=1000; print=TRUE; doses.plot=seq(from=0, to=1000, by=0.1)
  # df <- data.frame(Doses=c(80, 120, 150, 150, 180, 200), Alive=c(10, 12, 8, 6, 2, 1), Dead=c(0, 1, 5, 6, 9, 15))
  
  range.CI.qnorm <- qnorm(1-((1-range.CI)/2))
  
  if (!is.null(df)) {
    if ((!inherits(df, "data.frame")) & (!inherits(df, "matrix"))) {
      stop("df parameter must be a data.frame or a matrix")
    }
    
    namesdf <- tolower(colnames(df))
    alive <- NULL
    dead <- NULL
    N <- NULL
    
    if (any(namesdf=="alive")) alive <- df[, which(namesdf=="alive")]
    if (any(namesdf=="dead")) dead <- df[, which(namesdf=="dead")]
    if (any(namesdf=="n")) N <- df[, which(namesdf=="n")]
    if (any(namesdf=="doses")) doses <- df[, which(namesdf=="doses")]
  }
  
  if (is.null(doses) ) {
    stop("Doses must be provided")
  }
  
  
  cpt1 <- ifelse(is.null(alive), 0, 1)+ifelse(is.null(dead), 0, 1)+ifelse(is.null(N), 0, 1)
  if (cpt1<2) {
    stop("At least two informations from alive, dead or N must be provided")
  }
  
  if (is.null(alive)) alive <- N-dead
  if (is.null(dead)) dead <- N-alive
  if (is.null(N)) N <- alive+dead
  
  if (length(doses)!=length(N)) {
    stop("A dose must be provided for each experiment")
  }
  
  equation <- tolower(equation)
  
  equation <- match.arg(equation, choices = c("gsd", "hulin", "logistic", "logit", "probit", "hill", "flexit", "richards", "double-richards"))
  
  if (is.null(parameters.initial)) parameters.initial <- c(P=NA, S=NA, K=0, K1=1, K2=0)
  
  if (equation=="double-richards") parameters.initial["K2"] <- 1
  if (equation=="hulin") parameters.initial["K1"] <- 1/max(doses)
  
  
  par <- parameters.initial
  if (is.na(par["P"])) {
    if ((equation!="probit") & (equation!="logit")) {
      par["P"] <- doses[which.min(abs((alive/(alive+dead))-0.5))]
    } else {
      par["P"] <- 0
    }
  }
  
  if (is.na(par["S"])) {
    if ((equation!="probit")) {
    pente <- lm(alive/(alive+dead) ~ doses)
    par["S"] <- pente$coefficients["doses"]
    if (equation!="flexit") par["S"] <- 1/(4*par["S"])
    } else {
      par["S"] <- 0
    }
  }
  
  LD50_fit <- getFromNamespace(".LD50_fit", ns="HelpersMG")
  modelLD50 <- getFromNamespace(".modelLD50", ns="HelpersMG")
  
  if (equation!="richards") par <- par[which(names(par)!="K")]
  if (equation!="hulin" & equation!="double-richards" & equation!="flexit") par <- par[which(names(par)!="K1")]
  if (equation!="hulin" & equation!="double-richards" & equation!="flexit") par <- par[which(names(par)!="K2")]
  
  # Je commence toujours toujours par un glm
  if (equation=="probit") {
    inter <- glm(cbind(dead, alive) ~ doses, family=binomial(link="probit"))
    par["P"] <- unname(inter$coefficients["(Intercept)"])
    par["S"] <- unname(inter$coefficients["doses"])
  } else {
    inter <- glm(cbind(dead, alive) ~ doses, family=binomial(link="logit"))
    if (equation=="logit") {
      par["P"] <- unname(inter$coefficients["(Intercept)"])
      par["S"] <- unname(inter$coefficients["doses"])
    } else {
      a <- inter$coefficients["(Intercept)"]
      b <- inter$coefficients["doses"]
      par["P"] <- unname(-a/b)
      if (equation=="flexit") {
        par["S"] <- 1/(4*unname(-1/b))
      } else {
        par["S"] <- unname(-1/b)
      }
    }
  }
  
  repeat {
    result  <- optim(par, LD50_fit, fixed.parameters=fixed.parameters, alive=alive, 
                     N=N, doses=doses, equation=equation, method="BFGS", hessian=TRUE, 
                     control = list(maxit=1000))
    
    if (result$convergence==0) break
    par <- result$par
    if (print) print("Convergence is not achieved. Optimization continues !")
  }
  
  par <- c(result$par, fixed.parameters)
  
  if (!is.null(SE)) {
    res <- SE
  } else {
    mathessian <- result$hessian
    
    inversemathessian <- try(solve(mathessian), silent=TRUE)
    
    if (inherits(inversemathessian, "try-error")) {
      res <- rep(NA, length(par))
      names(res) <- names(par)
    } else {
      res <- abs(diag(inversemathessian))
    }
  }
  
  result$SE <- res
  result$AIC <- 2*result$value+2*length(result$par)
  
  
  result$range.CI <- range.CI
  result$l <- l
  
  l.l.TRD.c <- NULL
  l.h.TRD.c <- NULL
  TRD.c <- NULL
  #   plist <- NULL
  
  rep <- replicates-1
  df_par <- data.frame(P=unname(c(par["P"], rnorm(rep, par["P"], res["P"]))), 
                       S=unname(c(par["S"], rnorm(rep, par["S"], res["S"]))),
                       K=unname(ifelse(is.na(par["K"]), rep(NA, replicates), c(par["K"], rnorm(rep, par["K"], ifelse(is.na(res["K"]), 0, res["K"]))))),
                       K1=unname(ifelse(is.na(par["K1"]), rep(NA, replicates), c(par["K1"], rnorm(rep, par["K1"], ifelse(is.na(res["K1"]), 0, res["K1"]))))),
                       K2=unname(ifelse(is.na(par["K2"]), rep(NA, replicates), c(par["K2"], rnorm(rep, par["K2"], ifelse(is.na(res["K2"]), 0, res["K2"]))))))
  
  # sens==TRUE si descend
  sens <- (modelLD50(par, limit.low.TRD.minimum, equation)>modelLD50(par, limit.high.TRD.maximum, equation))
  
  out_LD50 <- apply(df_par, 1, function(par) {
    # marche en doses mais par en IP - corrige 9/9/2014
    limit.low.TRD <- c(limit.low.TRD.minimum, limit.high.TRD.maximum)
    limit.high.TRD <- c(limit.low.TRD.minimum, limit.high.TRD.maximum)
    P.TRD <- c(limit.low.TRD.minimum, limit.high.TRD.maximum)
    if (sens == (modelLD50(par, limit.low.TRD.minimum, equation)>modelLD50(par, limit.high.TRD.maximum, equation))) {
      
      
      
      repeat {
        doses.1 <- seq(from=limit.low.TRD[1], to=limit.low.TRD[2], length=3)
        p1 <- modelLD50(par, doses.1, equation)
        doses.2 <- seq(from=limit.high.TRD[1], to=limit.high.TRD[2], length=3)
        p2 <- modelLD50(par, doses.2, equation)
        doses.3 <- seq(from=P.TRD[1], to=P.TRD[2], length=3)
        p3 <- modelLD50(par, doses.3, equation)
        
        if (sens) {
          if (p1[2]<(1-l)) {
            limit.low.TRD <- doses.1[1:2]
          } else {
            limit.low.TRD <- doses.1[2:3]
          }
          if (p2[2]<l) {
            limit.high.TRD <- doses.2[1:2]
          } else {
            limit.high.TRD <- doses.2[2:3]
          }
          if (p3[2]<0.5) {
            P.TRD <- doses.3[1:2]
          } else {
            P.TRD <- doses.3[2:3]
          }
          
        } else {
          if (p1[2]<l) {
            limit.low.TRD <- doses.1[1:2]
          } else {
            limit.low.TRD <- doses.1[2:3]
          }
          if (p2[2]<(1-l)) {
            limit.high.TRD <- doses.2[1:2]
          } else {
            limit.high.TRD <- doses.2[2:3]
          }
          if (p3[2]>0.5) {
            P.TRD <- doses.3[1:2]
          } else {
            P.TRD <- doses.3[2:3]
          }
          
        }
        
        if (max(c(P.TRD[2]-P.TRD[1], limit.high.TRD[2]-limit.high.TRD[1], 
                  limit.low.TRD[2]-limit.low.TRD[1]))<1E-4) break
        
      }
      return(data.frame(l.l.TRD.c=mean(limit.low.TRD), l.h.TRD.c=mean(limit.high.TRD), 
                        TRD.c=mean(limit.high.TRD-limit.low.TRD), P.TRD=mean(P.TRD)))
      
    } else {
      return(data.frame(l.l.TRD.c=NA, l.h.TRD.c=NA, 
                        TRD.c=NA, P.TRD=NA))
      
    }
  })
  
  out_LD502 <- t(sapply(out_LD50, c))
  
  if (!is.null(doses.plot)) {
    out_LD50_plot <- apply(df_par, 1, function(par) {
      p <- modelLD50(par, doses.plot, equation)
      return(list(sr=p))
    }
    )
    
    
    out_LD502_plot <- sapply(out_LD50_plot, c, simplify = TRUE)
    # dans un df chaque colonne est une temperature
    out_LD503_plot <- t(sapply(out_LD502_plot, c, simplify = TRUE))
    #      dim(out_tsd3_plot)
    outquant <- apply(out_LD503_plot, 2, quantile, probs = c(0.025, 0.975))
    outquant <- rbind(outquant, mean=out_LD503_plot[1,], doses=doses.plot)
    result$CI <- outquant
    
  }
  
  
  result$LD50 <- as.numeric(out_LD502[1, "P.TRD"])
  result$SE_LD50 <- sd(as.numeric(out_LD502[, "P.TRD"]), na.rm = TRUE)
  result$TRD <- as.numeric(out_LD502[1, "TRD.c"])
  result$TRD_limits <- as.numeric(c(out_LD502[1, "l.l.TRD.c"], out_LD502[1, "l.h.TRD.c"]))
  result$SE_TRD <- sd(as.numeric(out_LD502[, "TRD.c"]), na.rm = TRUE)
  result$SE_TRD_limits <- c(sd(as.numeric(out_LD502[, "l.l.TRD.c"]), na.rm = TRUE), sd(as.numeric(out_LD502[, "l.h.TRD.c"]), na.rm = TRUE))
  limit.low.TRD <- result$TRD_limits[1]-range.CI.qnorm*result$SE_TRD_limits[1]
  limit.high.TRD <- result$TRD_limits[2]+range.CI.qnorm*result$SE_TRD_limits[2]
  
  
  saturated <- 0
  for (i in 1:length(N)) saturated <- saturated- dbinom(alive[i], N[i], alive[i]/N[i], log=TRUE)
  LRT <- 2*(result$value-saturated)
  result$GOF <- 1-pchisq(LRT, length(N)-length(par))
  
  result$alive <- alive
  result$dead <- dead
  result$N <- N
  result$doses <- doses
  result$equation <- equation
  result$l <- l
  result$fixed.parameters <- fixed.parameters
  
  if (print) {
    print(paste("The goodness of fit test is", sprintf("%.5f",result$GOF)))
    print(paste("The LD50 is", sprintf("%.3f",result$LD50), "SE", sprintf("%.3f",result$SE_LD50)))
    print(paste("The transitional range of doses is", sprintf("%.3f",result$TRD), "SE", sprintf("%.3f",result$SE_TRD)))
    print(paste("The lower limit of transitional range of doses is", sprintf("%.3f",result$TRD_limits[1]), "SE", sprintf("%.3f",result$SE_TRD_limits[1])))
    print(paste("The higher limit of transitional range of doses is", sprintf("%.3f",result$TRD_limits[2]), "SE", sprintf("%.3f",result$SE_TRD_limits[2])))
  }
  
  result <- addS3Class(result, "LD50")
  
  return(invisible(result))
}

