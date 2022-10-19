#' dSnbinom returns the density for the sum of random variable with negative binomial distributions
#' @title Density for the sum of random variable with negative binomial distributions. 
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @references Furman, E., 2007. On the convolution of the negative binomial random variables. Statistics & Probability Letters 77, 169-172.
#' @references Vellaisamy, P. & Upadhye, N.S. 2009. On the sums of compound negative binomial and gamma random variables. Journal of Applied Probability, 46, 272-283.
#' @return dSnbinom gives the density
#' @param x vector of (non-negative integer) quantiles.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param mu alternative parametrization via mean.
#' @param log	logical; if TRUE, probabilities p are given as log(p).
#' @param tol Tolerance for recurrence for Furman (2007) method
#' @param infinite Maximum level of recursion if tol is not reached for Furman (2007) method
#' @param method Can be Furman (default), Vellaisamy&Upadhye, or approximate.RandomObservations
#' @param n.random Number of random numbers used to estimate parameters of distribution for approximate methods
#' @param verbose Give information on the method
#' @description Density for the sum of random variable with negative binomial distributions.\cr
#' If all prob values are the same, infinite is automatically set to 0.\cr
#' Estimate using Vellaisamy&Upadhye method uses parallel computing depending on option mc.cores that can be 
#' defined using options(mc.cores = c) with c being the number of cores to be used. By default it will 
#' use all the available cores. Forking will be used in Unix system and no forking on Windows systems.
#' @family Distribution of sum of random variable with negative binomial distributions
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' alpha <- c(1, 2, 5, 1, 2)
#' p <- c(0.1, 0.12, 0.13, 0.14, 0.14)
#' dSnbinom(20, size=alpha, prob=p)
#' dSnbinom(20, size=alpha, prob=p, log=TRUE)
#' dSnbinom(20, size=2, mu=c(0.01, 0.02, 0.03))
#' dSnbinom(20, size=2, mu=c(0.01, 0.02, 0.03), log=TRUE)
#' # Test with a single distribution
#' dSnbinom(20, size=1, mu=20)
#' # when only one distribution is available, it is the same as dnbinom()
#' dnbinom(20, size=1, mu=20)
#' # If a parameter is supplied as only one value, it is supposed to be constant
#' dSnbinom(20, size=1, mu=c(14, 15, 10))
#' # The function is vectorized:
#' plot(0:200, dSnbinom(0:200, size=alpha, prob=p, method="furman"), bty="n", type="h", 
#'      xlab="x", ylab="Density")
#'      
#' # Comparison with simulated distribution using rep replicates
#' alpha <- c(2.1, 2.05, 2)
#' mu <- c(10, 30, 20)
#' rep <- 100000
#' distEmpirique <- rSnbinom(rep, size=alpha, mu=mu)
#' tabledistEmpirique <- rep(0, 301)
#' names(tabledistEmpirique) <- as.character(0:300)
#' tabledistEmpirique[names(table(distEmpirique))] <- table(distEmpirique)/rep
#' 
#' plot(0:300, dSnbinom(0:300, size=alpha, mu=mu, method="furman"), type="h", bty="n", 
#'    xlab="x", ylab="Density", ylim=c(0,0.02))
#' plot_add(0:300, tabledistEmpirique, type="l", col="red")
#' legend(x=200, y=0.02, legend=c("Empirical", "Theoretical"), 
#'    text.col=c("red", "black"), bty="n")
#' 
#' 
#' # Example from Vellaisamy, P. & Upadhye, N.S. (2009) - Table 1
#' # Note that value for k = 7 is not estimated because it is too long
#' k <- 2:6
#' x <- c(3, 5, 8, 10, 15)
#' table1_Vellaisamy <- matrix(NA, ncol=length(x), nrow=length(k))
#' rownames(table1_Vellaisamy) <- paste0("n = ", as.character(k))
#' colnames(table1_Vellaisamy) <- paste0("x = ", as.character(x))
#' table1_approximateObservations <- table1_Vellaisamy
#' table1_Furman <- table1_Vellaisamy
#' 
#' st_Furman <- rep(NA, length(k))
#' st_approximateObservations <- rep(NA, length(k))
#' st_Vellaisamy <- rep(NA, length(k))
#' 
#' for (n in k) {
#'     print(n)
#'     alpha <- 1:n
#'     p <- (1:n)/10
#'     st_Vellaisamy[which(n == k)] <- 
#'         system.time({
#'         table1_Vellaisamy[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="Vellaisamy&Upadhye", log=FALSE, verbose=TRUE)
#'         })[1]
#'     st_approximateObservations[which(n == k)] <- 
#'         system.time({
#'         table1_approximateObservations[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="approximate.RandomObservations", log=FALSE)
#'         })[1]
#'     st_Furman[which(n == k)] <- 
#'         system.time({
#'             table1_Furman[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="Furman", log=FALSE)
#'         })[1]
#' }
#' 
#' cbind(table1_Vellaisamy, st_Vellaisamy)
#' cbind(table1_Furman, st_Furman)
#' cbind(table1_approximateObservations, st_approximateObservations)
#' 
#' # Test of different methods
#' n <- 7
#' x <- 8
#' alpha <- 1:n
#' p <- (1:n)/10
#' dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, verbose=TRUE)
#' dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE)
#' dSnbinom(x=x, prob=p, size=alpha, method="approximate.RandomObservations", 
#'          log=FALSE)
#' 
#' # Test of different methods
#' n <- 2
#' x <- 15
#' alpha <- 1:n
#' p <- (1:n)/10
#' dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, verbose=TRUE)
#' dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, verbose=TRUE)
#' dSnbinom(x=x, prob=p, size=alpha, method="approximate.RandomObservations", 
#'          log=FALSE, verbose=TRUE)
#' 
#' 
#' # Test of different methods
#' alpha <- c(2.05, 2)
#' mu <- c(10, 30)
#' test <- rSnbinom(n=100000, size=alpha, mu=mu)
#' plot(0:200, table(test)[as.character(0:200)]/sum(table(test), na.rm=TRUE), 
#'      bty="n", type="h", xlab="x", ylab="Density")
#' lines(x=0:200, dSnbinom(0:200, size=alpha, mu=mu, log=FALSE, method="Furman"), col="blue")
#' lines(x=0:200, y=dSnbinom(0:200, size=alpha, mu=mu, log=FALSE, 
#'                           method="vellaisamy&upadhye"), col="red")
#' lines(x=0:200, y=dSnbinom(0:200, size=alpha, mu=mu, log=FALSE, 
#'                           method="approximate.randomobservations"), col="green")
#' 
#' # Another example
#' set.seed(2)
#' mutest <- c(56, 6.75, 1)
#' ktest <- c(50, 50, 50)
#' nr <- 100000
#' test <- rSnbinom(nr, size=ktest, mu=mutest)
#' pr_vellaisamy <- dSnbinom(0:150, size=ktest, mu=mutest, method = "vellaisamy&upadhye")
#' pr_furman <- dSnbinom(0:150, size=ktest, mu=mutest, method = "furman")
#' pr_approximateObservations <- dSnbinom(0:150, size=ktest, mu=mutest, 
#'                                        method = "approximate.randomobservations")
#' 
#' plot(table(test), xlab="N", ylab="Density", las=1, bty="n", ylim=c(0, 4000), xlim=c(0, 150))
#' lines(0:150, pr_vellaisamy*nr, col="red")
#' lines(0:150, pr_furman*nr, col="blue")
#' lines(0:150, pr_approximateObservations*nr, col="green")
#' 
#' dSnbinom(42, size=ktest, mu=mutest, method = "vellaisamy&upadhye")
#' dSnbinom(42, size=ktest, mu=mutest, method = "Furman")
#' dSnbinom(42, size=ktest, mu=mutest, method = "approximate.randomobservations")
#' 
#' # example to fit a distribution
#' data <- rnbinom(1000, size=1, mu=10)
#' hist(data)
#' ag <- rep(1:100, 10)
#' r <- aggregate(data, by=list(ag), FUN=sum)
#' hist(r[,2])
#' 
#' parx <- c(size=1, mu=10)
#' 
#' dSnbinomx <- function(x, par) {
#'   -sum(dSnbinom(x=x[,2], mu=rep(par["mu"], 10), size=par["size"], log=TRUE))
#' }
#' 
#' fit_mu_size <- optim(par = parx, fn=dSnbinomx, x=r, method="BFGS", control=c(trace=TRUE))
#' fit_mu_size$par
#' }
#' @export

dSnbinom <- function(x = stop("You must provide a x value")       , 
                     size = NULL                                  , 
                     prob = NULL                                  , 
                     mu = NULL                                    , 
                     log = FALSE                                  ,  
                     tol = 1E-6                                   ,
                     infinite = 1E6                               , 
                     method="Furman"                              ,
                     n.random = 1E6                               , 
                     verbose = FALSE                              ) {
  
  method <- tolower(method)
  method <- match.arg(arg=method, choices = c("furman", "vellaisamy&upadhye", 
                                              "approximate.randomobservations"))
  
  # prob=NULL; mu=NULL; log = FALSE; infinite=10
  if (is.null(mu) + is.null(size) + is.null(prob) != 1) stop("Two values among mu, size and prob must be provided")
  
  m <- max(c(length(size), length(prob), length(mu)))
  if (!is.null(mu)) mu <- rep(mu, m)[1:m]
  if (!is.null(size)) size <- rep(size, m)[1:m]
  if (!is.null(prob)) prob <- rep(prob, m)[1:m]
  
  if (is.null(prob)) {
    prob <- size/(size+mu)
    prob <- ifelse(prob>1-(1e-9), 1-1e-6, prob)
  }
  if (is.null(mu)) mu <- size/prob - size
  if (is.null(size))  size  <- (prob * mu) / (1 - prob)
  
  if (all(prob == prob[1])) {
    
    #    if (length(prob)<length(size)) prob <- rep(prob, length(size))[1:length(size)]
    #    if (length(size)<length(prob)) size <- rep(size, length(prob))[1:length(prob)]
    if (verbose) message("Exact method because all prob parameters are the same")
    return(dnbinom(x, size=sum(size), prob=prob[1], log=log))
  }
  if (method == "furman") {
    if (verbose) message("Furman (2007) method")
    alpha <- size
    p <- prob
    
    q <- 1-p
    p1 <- max(p)
    q1 <- 1-p1
    
    # R <- prod(((q*p1)/(q1*p))^(-alpha))
    R <- sum(log(((q*p1)/(q1*p))^(-alpha)))
    
    delta <- c(1, rep(NA, infinite))
    xi <- rep(NA, infinite)
    
    i <- 1
    repeat {
      xi[i] <- sum((alpha*(1-((q1*p)/(q*p1)))^i)/i)
      Ks <- 1:i
      delta[i+1] <- (1/i)*sum(Ks*xi[Ks]*delta[i-Ks+1])
      if ((abs(delta[i+1] - delta[i]) < tol) | (i == infinite)) {
        delta <- delta[!is.na(delta)]
        break
      }
      i <- i + 1
    }
    
    Pr <- sapply(x, function(S) {
      PrS <- sum(delta*dnbinom(S, size=sum(alpha)+seq_along(delta)-1, prob=p1))
      if (log) {
        PrS <- R + log(PrS)
      } else {
        PrS <- exp(R) * PrS
      }
      return(PrS)
    }
    )
    return(Pr)
  } 
  
  if (method == "approximate.randomobservations") {
    if (verbose) message("Approximate method with probabilities of observations.")
    test <- rSnbinom(n=n.random, size = size, mu=mu)
    return(sapply(X = x, FUN=function(y) sum(test==y)/n.random))
  }
  
  # if (method == "approximate.randomdistribution") {
  #   if (verbose) message("Approximate method with distribution modeled from observations.")
  #   test <- rSnbinom(n=n.random, size = size, mu=mu)
  #   mu <- mean(test)
  #   size <- mu^2 / (var(test) - mu)
  #   return(dnbinom(x=x, size = size, mu=mu, log = log))
  # }
  
  if (verbose) message("Vellaisamy & Upadhye (2009) method")
  if ((m > 6) & verbose) message("The Vellaisamy method with more than 6 summed distributions can be very slow and produces out of memory error.\n")
  TS <- NULL
  
  for (xec in x) {
    
    xx <- 0:xec
    nx <- xec + 1
    orep <- nx ^ m
    rep.fac <- 1L
    df <- matrix(data = NA, ncol=1, nrow = orep)
    rownames(df) <- as.character(1:orep)
    
    x2 <- xx
    orep <- orep/nx
    df[, 1]  <- x2[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, 
                                                        nx)), orep)]
    rep.fac <- rep.fac * nx
    if (m != 1) {
      for (i in 2:m) {
        x2 <- xx
        orep <- orep/nx
        
        x2 <- x2[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, 
                                                      nx)), orep)]
        df <- cbind(df, matrix(data = NA, ncol=1, nrow = nrow(df)))
        df[, i] <- x2[as.numeric(rownames(df))]
        
        if (dim(df)[1] == 1) {
          dd <- (sum(df[1, 1:i]) <= xec)
        } else {
          dd <- (rowSums(df[, 1:i]) <= xec)
        }
        df <- df[dd, , drop = FALSE]
        rep.fac <- rep.fac * nx
      }
    }
    if (dim(df)[1] == 1) {
      df <- df[sum(df[1, 1:i]) == xec, , drop = FALSE]
    } else {
      df <- df[rowSums(df) == xec, , drop = FALSE]
    }
    
    if (verbose) {
      message(paste0(as.character(nrow(df)), " combinations of ", as.character(m), " values produced a sum of ", as.character(xec), ". A total of ", as.character(nrow(df)*nx), " iterations will be necessary.\n"))
    }
    
    L <- universalmclapply(1:nrow(df), FUN = function(rw) {
      n <- df[rw, , drop=TRUE]
      sum(sapply(1:m, FUN = function(j) dnbinom(x=n[j], size = size[j], mu=mu[j], log=TRUE)))
    }, mc.cores = getOption("mc.cores", parallel::detectCores()), 
    clusterExport=list(varlist=c("df", "size", "mu"), envir=environment()))
    TS <- c(TS, sum(exp(unlist(L))))
  }
  TS <- unname(TS)
  if (log) TS <- log(TS)
  return(TS)
  
}

