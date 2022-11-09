#' Distribution of the Sum of Negative Binomial.
#' @title Distribution of the sum of random variable with negative binomial distributions. 
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @references Furman, E., 2007. On the convolution of the negative binomial random variables. Statistics & Probability Letters 77, 169-172.
#' @references Vellaisamy, P. & Upadhye, N.S. 2009. On the sums of compound negative binomial and gamma random variables. Journal of Applied Probability, 46, 272-283.
#' @references Girondot M, Barry J. Submitted. On the computation of the distribution of the sum of independent negative binomial random variables.
#' @return \code{dSnbinom} gives the density, \code{pSnbinom} gives the distribution function, 
#' \code{qSnbinom} gives the quantile function, and \code{rSnbinom} generates random deviates.
#' @param x vector of (non-negative integer) quantiles.
#' @param n number of observations.
#' @param q vector of quantiles.
#' @param p	vector of probabilities.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param mu alternative parametrization via mean.
#' @param log,log.p logical; if TRUE, probabilities \emph{p} are given as \emph{log(p)}.
#' @param tol Tolerance for recurrence for Furman (2007) method.
#' @param method Can be Furman (default), Vellaisamy&Upadhye, approximate.normal, approximate.negativebinomial or approximate.RandomObservations
#' @param max.iter Number of maximum iterations for Furman method. Can be NULL.
#' @param mean Mean of the distribution for approximate.normal method. If NULL, the theoretical mean will be used.
#' @param sd Standard deviation of the distribution for approximate.normal method. If NULL, the theoretical sd will be used.
#' @param n.random Number of random numbers used to estimate parameters of distribution for approximate.RandomObservations method.
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @param parallel logical; if FALSE (default), parallel computing is not used for Vellaisamy&Upadhye methods.
#' @param verbose Give more information on the method.
#' @description Distribution of the sum of random variable with negative binomial distributions.\cr
#' \code{dSnbinom} returns the density for the sum of random variable with negative binomial distributions.\cr
#' \code{pSnbinom} returns the distribution function for the sum of random variable with negative binomial distributions.\cr
#' \code{qSnbinom} returns the quantile function for the sum of random variable with negative binomial distributions.\cr
#' \code{rSnbinom} returns random numbers for the sum of random variable with negative binomial distributions.\cr
#' If all prob values are the same, exact probabilities are estimated.\cr
#' Estimate using \code{Vellaisamy&Upadhye} method uses parallel computing 
#' depending on value of \code{parallel}. The number of cores in usage can be 
#' defined using \code{options(mc.cores = c)} with \code{c} being the number of cores to be used. By default it will 
#' use all the available cores. Forking will be used in Unix system and no forking on Windows systems.\cr
#' When \code{Furman} method is in use, it will return the progress of Pr(S = x) during recursion 
#' in an attribute if verbose is TRUE (see examples).
#' @family Distributions
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' alpha <- c(1, 2, 5, 1, 2)
#' p <- c(0.1, 0.12, 0.13, 0.14, 0.14)
#' # By default, the Furman method with tol=1E-12 is used
#' dSnbinom(20, size=alpha, prob=p)
#' # Note the attribute is the dynamics of convergence of Pr(X=x)
#' attributes(dSnbinom(20, size=alpha, prob=p, verbose=TRUE))$Pk[, 1]
#' 
#' # Exact probability
#' dSnbinom(20, size=2, mu=c(0.01, 0.02, 0.03), method="vellaisamy&upadhye")
#' # By default, with Furman method and tol=1E-12, when probability 
#' # is very low, it will be biased
#' dSnbinom(20, size=2, mu=c(0.01, 0.02, 0.03), method="Furman")
#' # The solution is to use a tolerance lower than the estimate
#' dSnbinom(20, size=2, mu=c(0.01, 0.02, 0.03), method="Furman", tol=1E-40)
#' # Or a huge number of iterations
#' dSnbinom(20, size=2, mu=c(0.01, 0.02, 0.03), method="Furman", max.iter=10000)
#' 
#' # Test with a single distribution
#' dSnbinom(20, size=1, mu=20)
#' # when only one distribution is available, it is the same as dnbinom()
#' dnbinom(20, size=1, mu=20)
#' 
#' # If a parameter is supplied as only one value, it is supposed to be constant
#' dSnbinom(20, size=1, mu=c(14, 15, 10))
#' 
#' # The functions are vectorized:
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
#' # Note that value for k = 7 is very long
#' k <- 2:7
#' x <- c(3, 5, 8, 10, 15)
#' table1_Vellaisamy <- matrix(NA, ncol=length(x), nrow=length(k))
#' rownames(table1_Vellaisamy) <- paste0("n = ", as.character(k))
#' colnames(table1_Vellaisamy) <- paste0("x = ", as.character(x))
#' table1_approximateObservations <- table1_Vellaisamy
#' table1_Furman3 <- table1_Vellaisamy
#' table1_Furman6 <- table1_Vellaisamy
#' table1_Furman9 <- table1_Vellaisamy
#' table1_Furman12 <- table1_Vellaisamy
#' table1_Vellaisamy_parallel <- table1_Vellaisamy
#' table1_Approximate_Normal <- table1_Vellaisamy
#' 
#' st_Furman3 <- rep(NA, length(k))
#' st_Furman6 <- rep(NA, length(k))
#' st_Furman9 <- rep(NA, length(k))
#' st_Furman12 <- rep(NA, length(k))
#' st_approximateObservations <- rep(NA, length(k))
#' st_Vellaisamy <- rep(NA, length(k))
#' st_Vellaisamy_parallel <- rep(NA, length(k))
#' st_Approximate_Normal <- rep(NA, length(k))
#' 
#' for (n in k) {
#'     print(n)
#'     alpha <- 1:n
#'     p <- (1:n)/10
#'     st_Vellaisamy[which(n == k)] <- 
#'         system.time({
#'         table1_Vellaisamy[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="Vellaisamy&Upadhye", log=FALSE, verbose=FALSE)
#'         })[1]
#'     st_Vellaisamy_parallel[which(n == k)] <- 
#'         system.time({
#'         table1_Vellaisamy_parallel[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      parallel=TRUE, 
#'                                      method="Vellaisamy&Upadhye", log=FALSE, verbose=FALSE)
#'         })[1]
#'     st_approximateObservations[which(n == k)] <- 
#'         system.time({
#'         table1_approximateObservations[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="approximate.RandomObservations", log=FALSE, 
#'                                      verbose=FALSE)
#'         })[1]
#'     st_Furman3[which(n == k)] <- 
#'         system.time({
#'             table1_Furman3[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="Furman", tol=1E-3, log=FALSE, 
#'                                      verbose=FALSE)
#'         })[1]
#'     st_Furman6[which(n == k)] <- 
#'         system.time({
#'             table1_Furman6[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="Furman", tol=1E-6, log=FALSE, 
#'                                      verbose=FALSE)
#'         })[1]
#'     st_Furman9[which(n == k)] <- 
#'         system.time({
#'             table1_Furman9[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="Furman", tol=1E-9, log=FALSE, 
#'                                      verbose=FALSE)
#'         })[1]
#'     st_Furman12[which(n == k)] <- 
#'         system.time({
#'             table1_Furman12[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="Furman", tol=1E-12, log=FALSE, 
#'                                      verbose=FALSE)
#'         })[1]
#'         
#'     st_Approximate_Normal[which(n == k)] <- 
#'         system.time({
#'             table1_Approximate_Normal[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="approximate.normal", tol=1E-12, log=FALSE, 
#'                                      verbose=FALSE)
#'         })[1]
#' }
#' 
#' cbind(table1_Vellaisamy, st_Vellaisamy)
#' cbind(table1_Vellaisamy_parallel, st_Vellaisamy_parallel)
#' cbind(table1_Furman3, st_Furman3)
#' cbind(table1_Furman6, st_Furman6)
#' cbind(table1_Furman9, st_Furman9)
#' cbind(table1_Furman12, st_Furman12)
#' cbind(table1_approximateObservations, st_approximateObservations)
#' cbind(table1_Approximate_Normal, st_Approximate_Normal)
#' 
#' 
#' # Test of different methods
#' n <- 7
#' x <- 8
#' alpha <- 1:n
#' p <- (1:n)/10
#' 
#' # Parallel computing is not always performant
#' # Here it is approximately the same time of execution
#' system.time({print(dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, 
#'          verbose=TRUE, parallel=TRUE))})
#' system.time({print(dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, 
#'          verbose=TRUE, parallel=FALSE))})
#'          
#' # Test of different methods
#' n <- 7
#' x <- 15
#' alpha <- 1:n
#' p <- (1:n)/10
#' 
#' # Parallel computing is sometimes very performant
#' # Here parallel computing is 7 times faster (with 8 cores computer) 
#' #             for vellaisamy&upadhye method
#' system.time(dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, 
#'          verbose=TRUE, parallel=TRUE))
#' system.time(dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, 
#'          verbose=TRUE, parallel=FALSE))
#'          
#' # Test of different methods
#' n <- 2
#' x <- 3
#' alpha <- 1:n
#' p <- (1:n)/10
#' 
#' # Parallel computing is sometimes very performant
#' # Here parallel computing is 7 times faster (with 8 cores computer) 
#' #             for vellaisamy&upadhye method
#' system.time(dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, 
#'          verbose=TRUE, parallel=TRUE))
#' system.time(dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, 
#'          verbose=TRUE, parallel=FALSE))
#'          
#' # Test for different tolerant values
#' n <- 7
#' x <- 8
#' alpha <- 1:n
#' p <- (1:n)/10
#' dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, verbose=TRUE)
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-3, verbose=TRUE))
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-6, verbose=TRUE))
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-9, verbose=TRUE))
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-12, verbose=TRUE))
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="approximate.RandomObservations", 
#'          log=FALSE, verbose=TRUE))
#' 
#' # Test for criteria of convergence
#' Pr_exact <- dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", 
#'                    log=FALSE, verbose=TRUE)
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-12, 
#'                  verbose=TRUE)
#' Pr_exact;as.numeric(Pr_Furman)
#' plot(1:length(attributes(Pr_Furman)$Pk), 
#'      log10(abs(attributes(Pr_Furman)$Pk-Pr_exact)), type="l", xlab="Iterations", 
#'      ylab="Abs log10", bty="n")
#' lines(1:(length(attributes(Pr_Furman)$Pk)-1), 
#'       log10(abs(diff(attributes(Pr_Furman)$Pk))), col="red")
#' legend("bottomleft", legend=c("Log10 Convergence to true value", "Log10 Rate of change"), 
#'        col=c("black", "red"), 
#'        lty=1)
#'        
#' n <- 7
#' x <- 6
#' alpha <- 1:n
#' p <- (1:n)/10
#' Pr_exact <- dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", 
#'                    log=FALSE, verbose=TRUE)
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-12, 
#'                  verbose=TRUE)
#'                  
#'                  
#' # pdf("figure 2.pdf", width=7, height=7, pointsize=14)                 
#' ylab <- as.expression(bquote(.("P(S=6) x 10")^"6"))
#' layout(1:2)
#' par(mar=c(3, 4, 1, 1))
#' plot(1:length(attributes(Pr_Furman)$Pk), 
#'      attributes(Pr_Furman)$Pk*1E6, type="l", xlab="", 
#'      ylab="", bty="n", las=1, xlim=c(0, 80))
#' mtext(text=ylab, side=2, line=2.5)
#' segments(x0=25, x1=80, y0=Pr_exact*1E6, y1=Pr_exact*1E6, lty=3)
#' par(xpd=TRUE)
#' text(x=0, y=6, labels="Exact probability", pos=4)
#' txt <- "        Approximate probability\n    based on Furman (2007)\nrecursive iterations"
#' text(x=20, y=2, labels=txt, pos=4)
#' text(x=75, y=5, labels="A", cex=2)
#' par(mar=c(4, 4, 1, 1))
#' ylab <- as.expression(bquote("log"["10"]*""*"(P"["k+1"]*""*" - P"["k"]*""*")"))
#' plot(1:(length(attributes(Pr_Furman)$Pk)-1), 
#'       log10(diff(attributes(Pr_Furman)$Pk)), col="black", xlim=c(0, 80), type ="l", 
#'       bty="n", las=1, xlab="Iterations", ylab="")
#' mtext(text=ylab, side=2, line=2.5)
#'  peak <- (1:(length(attributes(Pr_Furman)$Pk)-1))[which.max(
#'             log10(abs(diff(attributes(Pr_Furman)$Pk))))]
#' segments(x0=peak, x1=peak, y0=-12, y1=-5, lty=2)
#' text(x=0, y=-6, labels="Positive trend", pos=4)
#' text(x=30, y=-6, labels="Negative trend", pos=4)
#' segments(x0=0, x1=43, y0=-12, y1=-12, lty=4)
#' segments(x0=62, x1=80, y0=-12, y1=-12, lty=4)
#' text(x=45, y=-12, labels="Tolerance", pos=4)
#' text(x=75, y=-7, labels="B", cex=2)
#' # dev.off()
#' 
#' # pdf("figure 2.pdf", width=7, height=7, pointsize=14)                 
#' ylab <- as.expression(bquote(.("P(S=6) x 10")^"6"))
#' layout(1:2)
#' par(mar=c(3, 4, 1, 1))
#' plot(1:length(attributes(Pr_Furman)$Pk), 
#'      attributes(Pr_Furman)$Pk*1E6, type="l", xlab="", 
#'      ylab="", bty="n", las=1, xlim=c(0, 80))
#' mtext(text=ylab, side=2, line=2.5)
#' segments(x0=25, x1=80, y0=Pr_exact*1E6, y1=Pr_exact*1E6, lty=3)
#' par(xpd=TRUE)
#' text(x=0, y=6, labels="Exact probability", pos=4)
#' txt <- "        Approximate probability\n    based on Furman (2007)\nrecursive iterations"
#' text(x=20, y=2, labels=txt, pos=4)
#' text(x=75, y=5, labels="A", cex=2)
#' par(mar=c(4, 4, 1, 1))
#' ylab <- as.expression(bquote("(P"["k+1"]*""*" - P"["k"]*""*") x 10"^"7"))
#' plot(1:(length(attributes(Pr_Furman)$Pk)-1), 
#'       diff(attributes(Pr_Furman)$Pk)*1E7, 
#'       col="black", xlim=c(0, 80), type ="l", 
#'       bty="n", las=1, xlab="Iterations", ylab="")
#' mtext(text=ylab, side=2, line=2.5)
#'  peak <- (1:(length(attributes(Pr_Furman)$Pk)-1))[which.max(diff(attributes(Pr_Furman)$Pk))]
#' segments(x0=peak, x1=peak, y0=0, y1=3.5, lty=2)
#' text(x=-2, y=3.5, labels="Positive trend", pos=4)
#' text(x=30, y=3.5, labels="Negative trend", pos=4)
#' segments(x0=0, x1=22, y0=1E-12, y1=1E-12, lty=4)
#' segments(x0=40, x1=80, y0=1E-12, y1=1E-12, lty=4)
#' text(x=22, y=1E-12+0.2, labels="Tolerance", pos=4)
#' text(x=75, y=3, labels="B", cex=2)
#' # dev.off()
#' 
#' # Test of different methods
#' n <- 2
#' x <- 15
#' alpha <- 1:n
#' p <- (1:n)/10
#' dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, verbose=TRUE)
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-3, verbose=TRUE))
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-6, verbose=TRUE))
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-9, verbose=TRUE))
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-12, verbose=TRUE))
#' dSnbinom(x=x, prob=p, size=alpha, method="approximate.RandomObservations", 
#'          log=FALSE, verbose=TRUE)
#'
#'
#' n <- 50
#' x <- 300
#' alpha <- (1:n)/100
#' p <- (1:n)/1000
#' # Produce an error
#' Pr_exact <- dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", 
#'                    log=FALSE, verbose=TRUE)
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-40, 
#'                  verbose=FALSE)
#' Pr_ApproximateNormal <- dSnbinom(x=x, prob=p, size=alpha, method="approximate.normal", 
#'                         log=FALSE, 
#'                         verbose=TRUE)
#' Pr_ApproximateRandom <- dSnbinom(x=x, prob=p, size=alpha, method="approximate.RandomObservations", 
#'                         log=FALSE, n.random=1E6, 
#'                         verbose=TRUE)
#'                         
#' n <- 500
#' x <- 3000
#' alpha <- (1:n)/100
#' p <- (1:n)/1000
#' # Produce an error
#' Pr_exact <- dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", 
#'                    log=FALSE, verbose=TRUE)
#' # Produce an error
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-40, 
#'                  verbose=FALSE)
#' Pr_ApproximateNormal <- dSnbinom(x=x, prob=p, size=alpha, method="approximate.normal", 
#'                         log=FALSE, 
#'                         verbose=TRUE)
#' Pr_ApproximateNegativeBinomial <- dSnbinom(x=x, prob=p, size=alpha, 
#'                         method="approximate.negativebinomial", 
#'                         log=FALSE, 
#'                         verbose=TRUE)
#' Pr_ApproximateRandom <- dSnbinom(x=x, prob=p, size=alpha, 
#'                         method="approximate.RandomObservations", 
#'                         log=FALSE, n.random=1E6, 
#'                         verbose=TRUE)
#'                
#'
#' # pdf("figure 1.pdf", width=7, height=7, pointsize=14)    
#' 
#' layout(1:2)
#' par(mar=c(3, 4, 1, 1))
#' # Normal approximation
#' alpha <- seq(from=10, to=100, length.out=10)
#' p <- seq(from=0.5, to=0.9, length.out=10)
#' test <- rSnbinom(100000, prob=p, size=alpha)
#' c(mean(test), sum(alpha*(1-p)/p))
#' c(var(test), sum(alpha*(1-p)/p^2))
#' c(sd(test), sqrt(sum(alpha*(1-p)/p^2)))
#' # plot(table(test)/length(test), las=1, bty="n", col="grey", xlab="x", 
#' #       ylab="P(X=x)")
#' plot(1, 1, las=1, bty="n", col="grey", xlab="", 
#'        xlim=c(120, 240), ylim=c(0, 0.025), 
#'        ylab="P(X=x)", type="n", yaxt="n")
#' axis(2, at=c(0, 0.01, 0.02), las=1)
#' p_nb <- dSnbinom(120:240, prob=p, size=alpha, method="Furman", verbose=TRUE)
#' segments(x0=(120:240), x1=(120:240), 
#'          y0=0, y1=as.numeric(p_nb), col="black")
#' p_n <- dSnbinom(120:240, prob=p, size=alpha, method="approximate.normal", verbose=TRUE)
#' points(120:240, p_n, col="black", pch=19, cex=0.5)
#' p_nb <- dSnbinom(120:240, prob=p, size=alpha, method="approximate.negativebinomial", verbose=TRUE)
#' points(120:240, p_nb, col="black", pch=4, cex=1)
#' text(x=ScalePreviousPlot(x = 0.95, y = 0.95)$x, 
#'      y=ScalePreviousPlot(x = 0.95, y = 0.95)$y, labels="A", cex=2)
#'      
#' text(x=ScalePreviousPlot(x = 0.0, y = 0.80)$x, 
#'      y=ScalePreviousPlot(x = 0.0, y = 0.80)$y, labels="|", cex=1)
#' text(x=ScalePreviousPlot(x = 0.05, y = 0.80)$x, 
#'      y=ScalePreviousPlot(x = 0.05, y = 0.80)$y, labels="Exact probability", cex=1, pos=4)
#'      
#' points(x=ScalePreviousPlot(x = 0.0, y = 0.87)$x, 
#'      y=ScalePreviousPlot(x = 0.0, y = 0.87)$y, pch=4, cex=1)
#' text(x=ScalePreviousPlot(x = 0.05, y = 0.87)$x, 
#'      y=ScalePreviousPlot(x = 0.05, y = 0.87)$y, labels="NB approximation", cex=1, pos=4)
#' 
#' points(x=ScalePreviousPlot(x = 0.0, y = 0.95)$x, 
#'      y=ScalePreviousPlot(x = 0.0, y = 0.95)$y, pch=19, cex=0.8)
#' text(x=ScalePreviousPlot(x = 0.05, y = 0.95)$x, 
#'      y=ScalePreviousPlot(x = 0.05, y = 0.95)$y, labels="Normal approximation", cex=1, pos=4)
#'      
#' # Normal approximation
#' n <- 2
#' x <- 15
#' alpha <- 1:n
#' p <- (1:n)/10
#' test <- rSnbinom(100000, prob=p, size=alpha)
#' c(mean(test), sum(alpha*(1-p)/p))
#' c(var(test), sum(alpha*(1-p)/p^2))
#' c(sd(test), sqrt(sum(alpha*(1-p)/p^2)))
#' # plot(table(test)/length(test), las=1, bty="n", col="grey", xlab="x", 
#' #      ylab="P(X=x)", ylim=c(0, 0.05))
#' par(mar=c(4, 4, 1, 1))
#' plot(1, 1, las=1, bty="n", col="black", xlab="x", 
#'        xlim=c(0, 80), ylim=c(0, 0.04), 
#'        ylab="P(X=x)", type="n")
#' p_nb <- as.numeric(dSnbinom(0:80, prob=p, size=alpha))
#' par(xpd=TRUE)
#' segments(x0=0:80, x1=0:80, 
#'          y0=0, y1=as.numeric(p_nb), col="black")
#' p_n <- dSnbinom(x=0:80, prob=p, size=alpha, mu=NULL, method="approximate.normal", 
#'                 verbose=TRUE)
#' points(0:80, p_n, col="black", pch=19, cex=0.5)
#' p_nb <- dSnbinom(x=0:80, prob=p, size=alpha, method="approximate.negativebinomial", 
#'                  verbose=TRUE)
#' points(x=0:80, p_nb, col="black", pch=4, cex=1)

#' text(x=ScalePreviousPlot(x = 0.95, y = 0.95)$x, 
#'         y=ScalePreviousPlot(x = 0.95, y = 0.95)$y, labels="B", cex=2)
#' 
#' # dev.off()
#' 
#' 
#' # Test for criteria of convergence
#' Pr_exact <- dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", 
#'                    log=FALSE, verbose=TRUE)
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-12, 
#'                  verbose=TRUE)
#' Pr_exact;as.numeric(Pr_Furman)
#' plot(1:length(attributes(Pr_Furman)$Pk), 
#'      log10(abs(attributes(Pr_Furman)$Pk-Pr_exact)), type="l", xlab="Iterations", 
#'      ylab="Abs log10", bty="n")
#' lines(1:(length(attributes(Pr_Furman)$Pk)-1), 
#'       log10(abs(diff(attributes(Pr_Furman)$Pk))), col="red")
#' legend("bottomleft", legend=c("Log10 Convergence to true value", "Log10 Rate of change"), 
#'        col=c("black", "red"), 
#'        lty=1)
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
#' # Test for criteria of convergence for x = 50
#' x <- 50
#' # Test for criteria of convergence
#' Pr_exact <- dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", 
#'                    log=FALSE, verbose=TRUE)
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, tol=1E-12, 
#'                  verbose=TRUE)
#' Pr_exact;as.numeric(Pr_Furman)
#' plot(1:length(attributes(Pr_Furman)$Pk), 
#'      log10(abs(attributes(Pr_Furman)$Pk-Pr_exact)), type="l", xlab="Iterations", 
#'      ylab="Abs log10", bty="n")
#' lines(1:(length(attributes(Pr_Furman)$Pk)-1), 
#'       log10(abs(diff(attributes(Pr_Furman)$Pk))), col="red")
#' legend("bottomleft", legend=c("Log10 Convergence to true value", "Log10 Rate of change"), 
#'        col=c("black", "red"), 
#'        lty=1)
#' 
#' # Another example more complicated
#' set.seed(2)
#' mutest <- c(56, 6.75, 1)
#' ktest <- c(50, 50, 50)
#' nr <- 100000
#' test <- rSnbinom(nr, size=ktest, mu=mutest)
#' system.time({pr_vellaisamy <- dSnbinom(x=0:150, size=ktest, mu=mutest, 
#'                           method = "vellaisamy&upadhye", verbose=FALSE, parallel=FALSE)})
#' # Parallel computing is not efficient
#' system.time({pr_vellaisamy <- dSnbinom(x=0:150, size=ktest, mu=mutest, 
#'                           method = "vellaisamy&upadhye", verbose=FALSE, parallel=TRUE)})
#' system.time({pr_furman <- dSnbinom(x=0:150, size=ktest, mu=mutest, prob=NULL, 
#'                       method = "furman", tol=1E-12, verbose=FALSE, log=FALSE)})
#' pr_approximateObservations <- dSnbinom(0:150, size=ktest, mu=mutest, 
#'                                        method = "approximate.randomobservations")
#' 
#' plot(table(test), xlab="N", ylab="Density", las=1, bty="n", ylim=c(0, 4000), xlim=c(0, 150))
#' lines(0:150, pr_vellaisamy*nr, col="red")
#' lines(0:150, pr_furman*nr, col="blue")
#' lines(0:150, pr_approximateObservations*nr, col="green")
#' 
#' dSnbinom(x=42, size=ktest, mu=mutest, prob=NULL, 
#'          method = "vellaisamy&upadhye", verbose=TRUE)
#' as.numeric(dSnbinom(x=42, size=ktest, mu=mutest, prob=NULL, 
#'            method = "Furman", verbose=TRUE))
#' dSnbinom(x=42, size=ktest, mu=mutest, prob=NULL, 
#'          method = "approximate.randomobservations", verbose=TRUE)
#' 
#' x <- 100
#' # Test for criteria of convergence
#' Pr_exact <- dSnbinom(x=x, size=ktest, mu=mutest, method="vellaisamy&upadhye", 
#'                    log=FALSE, verbose=TRUE)
#' Pr_Furman <- dSnbinom(x=x, size=ktest, mu=mutest, method="Furman", log=FALSE, tol=1E-12, 
#'                  verbose=TRUE)
#' Pr_exact;as.numeric(Pr_Furman)
#' plot(1:length(attributes(Pr_Furman)$Pk), 
#'      log10(abs(attributes(Pr_Furman)$Pk-Pr_exact)), type="l", xlab="Iterations", 
#'      ylab="Abs log10", bty="n", ylim=c(-100, 0))
#' lines(1:(length(attributes(Pr_Furman)$Pk)-1), 
#'       log10(abs(diff(attributes(Pr_Furman)$Pk))), col="red")
#' legend("bottomright", legend=c("Log10 Convergence to true value", "Log10 Rate of change"), 
#'        col=c("black", "red"), 
#'        lty=1)
#'        
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
#' 
#' alpha <- c(2.1, 2.05, 2)
#' mu <- c(10, 30, 20)
#' p <- pSnbinom(q=10, size=alpha, mu=mu, lower.tail = TRUE)
#' 
#' alpha <- c(2.1, 2.05, 2)
#' mu <- c(10, 30, 20)
#' q <- qSnbinom(p=0.1, size=alpha, mu=mu, lower.tail = TRUE)
#' 
#' alpha <- c(2.1, 2.05, 2)
#' mu <- c(10, 30, 20)
#' rep <- 100000
#' distEmpirique <- rSnbinom(n=rep, size=alpha, mu=mu)
#' tabledistEmpirique <- rep(0, 301)
#' names(tabledistEmpirique) <- as.character(0:300)
#' tabledistEmpirique[names(table(distEmpirique))] <- table(distEmpirique)/rep
#' 
#' plot(0:300, dSnbinom(0:300, size=alpha, mu=mu), type="h", bty="n", 
#'    xlab="x", ylab="Density", ylim=c(0,0.02))
#' plot_add(0:300, tabledistEmpirique, type="l", col="red")
#' legend(x=200, y=0.02, legend=c("Empirical", "Theoretical"), 
#'    text.col=c("red", "black"), bty="n")
#' }
#' @export

#' @describeIn Snbinom Density for the sum of random variable with negative binomial distributions.
# #' @section Another section after function section:


dSnbinom <- function(x = stop("You must provide a x value")       , 
                     size = NULL                                  , 
                     prob = NULL                                  , 
                     mu = NULL                                    , 
                     log = FALSE                                  ,  
                     tol = 1E-12                                  , 
                     method="Furman"                              ,
                     max.iter=NULL                                ,
                     mean=NULL                                    ,
                     sd=NULL                                      ,
                     n.random = 1E6                               , 
                     parallel = FALSE                             ,
                     verbose = FALSE                              ) {
  
  method <- tolower(method)
  method <- match.arg(arg=method, choices = c("furman", "vellaisamy&upadhye", 
                                              "approximate.randomobservations", 
                                              "approximate.normal", 
                                              "approximate.negativebinomial"))
  
  # prob=NULL; mu=NULL; log = FALSE
  if (is.null(mu) + is.null(size) + is.null(prob) != 1) stop("Two values exactly among mu, size and prob must be provided")
  
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
    if (verbose) message("Exact method because all prob parameters are the same.")
    return(dnbinom(x, size=sum(size), prob=prob[1], log=log))
  }
  if (method == "furman") {
    if (verbose) message("Furman (2007) method")
    alpha <- size
    p <- prob
    
    q <- 1-p
    p1 <- max(p)
    q1 <- 1-p1
    
    R <- sum(log(((q*p1)/(q1*p))^(-alpha)))
    
    delta <- 1
    xi <- NULL
    
    prePr <- rep(-Inf, length(x))
    # preDlt <- rep(-Inf, length(x))
    pr_tot <- NULL
    
    
    found <- rep(FALSE, length(x))
    Salpha <- sum(alpha)
    Prx_S <- dnbinom(x, size=Salpha, prob=p1)
    
    i <- 1
    repeat {
      # if (i == 545) stop()
      xi <- c(xi, sum((alpha*(1-((q1*p)/(q*p1)))^i)/i))
      Ks <- 1:i
      newdelta <- (1/i)*sum(Ks*xi[Ks]*delta[i-Ks+1])
      delta <- c(delta, newdelta)
      
      xf <- x[!found]
      
      Prx <- sapply(xf, function(S) {
        return(newdelta*dnbinom(S, size=Salpha+length(delta)-1, prob=p1))
      })
      Prx_S[!found] <- Prx_S[!found] + Prx
      Pr <- Prx_S
      
      if (log) {
        Pr <- R + log(Pr)
      } else {
        Pr <- exp(R) * Pr
      }
      
      pr_tot <- rbind(pr_tot, matrix(Pr, nrow=1))
      
      if (i>3) {
        if (is.null(max.iter)) {
          found <- found | ((abs(prePr - Pr) <= rate) & (all(abs(Pr - prePr) <= tol)))
          if (any(is.na(found))) break
          if (all(found)) break
        } else {
          if (i>=max.iter) break
        }
        if (verbose) {
          message(paste0("Iteration ", as.character(i), "; Max difference=", 
                         as.character(max(abs(Pr - prePr))), "; Trend of rate change ", 
                         ifelse((abs(prePr - Pr) <= rate)[which.max(abs(Pr - prePr))], "negative", "positive")))
        }
      }
      
      
      rate <- abs(prePr - Pr)
      prePr <- Pr
      # preDlt <- delta[i+1]
      
      i <- i + 1
    }
    
    if (any(is.na(found))) {
      warning("Error during recursions.")
      return() 
    }
    
    if (verbose) message("Loop has been stopped after ", as.character(i), " iterations.")
    
    # attributes(Pr) <- modifyList(as.list(attributes(Pr)), list(delta=delta, xi=xi, Pr=pr_tot))
    if (verbose) attributes(Pr) <- modifyList(as.list(attributes(Pr)), list("Pk"=pr_tot))
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
  
  if (method == "approximate.negativebinomial") {
    if (verbose) message("Approximate method with negative binomial distribution.")
    
    if (is.null(mean)) {
      mean <- sum(size*(1-prob)/prob)
    }
    
    if (is.null(sd)) {
      sd <- sqrt(sum(size*(1-prob)/prob^2))
    }
    
    size <- mean^2 / (sd^2 - mean)
    
    if (verbose) message(paste0("Mean=", specify_decimal(mean)))
    if (verbose) message(paste0("Standard deviation=", specify_decimal(sd)))
    
    p <- dnbinom(x = x, mu=mean, size=size, log=log)
    
    if (verbose) attributes(p) <- modifyList(as.list(attributes(p)), list(mean=mean, sd=sd))
    return(p)
  }
  
  if (method == "approximate.normal") {
    if (verbose) message("Approximate method with normal distribution.")
    
    if (is.null(mean)) {
      mean <- sum(size*(1-prob)/prob)
    }
    if (is.null(sd)) {
      sd <- sqrt(sum(size*(1-prob)/prob^2))
    }
    
    if (verbose) message(paste0("Mean=", specify_decimal(mean)))
    if (verbose) message(paste0("Standard deviation=", specify_decimal(sd)))
    # if (verbose) message(paste0("Index of quality of approximation=", specify_decimal(pnorm(q=0, mean=mean, sd=sd, lower.tail = FALSE, log.p = FALSE))))
    
    p <- pnorm(x+0.5, mean=mean, sd=sd, log.p=FALSE) - pnorm(x-0.5,  mean=mean, sd=sd, log.p=FALSE)
    if (log) p <- log(p)
    
    if (verbose) attributes(p) <- modifyList(as.list(attributes(p)), list(mean=mean, sd=sd))
    return(p)
  }
  
  if (verbose) message("Vellaisamy & Upadhye (2009) method")
  if ((m > 6) & verbose) message("The Vellaisamy method with more than 6 summed distributions can be very slow and produces out of memory error.")
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
      message(paste0(as.character(nrow(df)), " combinations of ", as.character(m), " values produced a sum of ", as.character(xec), ". A total of ", as.character(nrow(df)*nx), " iterations will be necessary."))
    }
    
    if (parallel) {
      if (verbose) {
        message(paste0("Parallel computing with ", as.character(getOption("mc.cores", parallel::detectCores())), " cores."))
      }
      suppressMessages(L <- universalmclapply(1:nrow(df), FUN = function(rw) {
        n <- df[rw, , drop=TRUE]
        sum(sapply(1:m, FUN = function(j) dnbinom(x=n[j], size = size[j], mu=mu[j], log=TRUE)))
      }, mc.cores = getOption("mc.cores", parallel::detectCores()), 
      clusterExport=list(varlist=c("df", "size", "mu"), envir=environment())))
      L <- unlist(L)
    } else {
      L <- sapply(1:nrow(df), FUN = function(rw) {
        n <- df[rw, , drop=TRUE]
        sum(sapply(1:m, FUN = function(j) dnbinom(x=n[j], size = size[j], mu=mu[j], log=TRUE)))
      }
      )
    }
    
    
    
    TS <- c(TS, sum(exp(L)))
  }
  TS <- unname(TS)
  if (log) TS <- log(TS)
  return(TS)
  
}

#' @export
#' @describeIn Snbinom Distribution function for the sum of random variable with negative binomial distributions.

pSnbinom <- function(q=stop("At least one quantile must be provided"), 
                     size=NULL, 
                     prob=NULL, mu=NULL, lower.tail = TRUE, log.p = FALSE, tol=1E-12, 
                     method="Furman") {
  
  method <- tolower(method)
  method <- match.arg(arg=method, choices = c("furman", "vellaisamy&upadhye", 
                                              "approximate.randomobservations", 
                                              "approximate.normal", 
                                              "approximate.negativebinomial"))
  
  
  # prob=NULL; mu=NULL; log = FALSE; infinite=10
  
  if (is.null(mu) + is.null(size) + is.null(prob) != 1) stop("Two values among mu, size and prob must be provided")
  
  m <- max(c(length(size), length(prob), length(mu)))
  if (!is.null(mu)) mu <- rep(mu, m)[1:m]
  if (!is.null(size)) size <- rep(size, m)[1:m]
  if (!is.null(prob)) prob <- rep(prob, m)[1:m]
  
  if (is.null(prob)) prob <- size/(size+mu)
  if (is.null(mu)) mu <- size/prob - size
  if (is.null(size))  size  <- (prob * mu) / (1 - prob)
  
  #  if (length(prob)<length(size)) prob <- rep(prob, length(size))[1:length(size)]
  #  if (length(size)<length(prob)) size <- rep(size, length(prob))[1:length(prob)]
  
  pp <- vapply(q, FUN=function(qq) {
    
    l <- dSnbinom(0:qq, prob=prob, size=size, mu=NULL, log=FALSE, 
                  tol=tol, method = method)
    p <- sum(l)
    if (!lower.tail) p <- 1-p
    if (log.p) p <- log(p)
    
    return(p)
  }, FUN.VALUE = 10.1)
  
  return(pp)
}

#' @export
#' @describeIn Snbinom Quantile function for the sum of random variable with negative binomial distributions.

qSnbinom <- function(p=stop("At least one probability must be provided"), 
                     size=stop("size parameter is mandatory"), 
                     prob=NULL, mu=NULL, lower.tail = TRUE, log.p = FALSE, 
                     tol=1E-12, method="Furman") {
  
  # prob=NULL; mu=NULL; log = FALSE; infinite=10
  
  method <- tolower(method)
  method <- match.arg(arg=method, choices = c("furman", "vellaisamy&upadhye", 
                                              "approximate.randomobservations", 
                                              "approximate.normal", 
                                              "approximate.negativebinomial"))
  
  
  if (is.null(mu) + is.null(size) + is.null(prob) != 1) stop("Two values among mu, size and prob must be provided")
  
  m <- max(c(length(size), length(prob), length(mu)))
  if (!is.null(mu)) mu <- rep(mu, m)[1:m]
  if (!is.null(size)) size <- rep(size, m)[1:m]
  if (!is.null(prob)) prob <- rep(prob, m)[1:m]
  
  if (is.null(prob)) prob <- size/(size+mu)
  if (is.null(mu)) mu <- size/prob - size
  if (is.null(size))  size  <- (prob * mu) / (1 - prob)  
  #  if (length(prob)<length(size)) prob <- rep(prob, length(size))[1:length(size)]
  #  if (length(size)<length(prob)) size <- rep(size, length(prob))[1:length(prob)]
  
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1-p
  
  
  
  
  p1 <- max(p)
  
  rankl <- 10
  limit <- 10
  repeat {
    if (pSnbinom(limit, prob=prob, size=size, mu=NULL, log.p=FALSE, tol=tol, method = method)>=p1) break
    limit <- limit + rankl
    rankl <- rankl + 5
  }
  
  seqp <- pSnbinom(0:limit, prob=prob, size=size, mu=NULL, log.p=FALSE, tol=tol, method = method)
  
  qq <- vapply(p, FUN=function(pp) {
    
    return(which(seqp>pp)[1]-1)
  }, FUN.VALUE = 1.1)
  
  return(qq)
}

#' @export
#' @describeIn Snbinom Random numbers for the sum of random variable with negative binomial distributions.

rSnbinom <- function(n=1, 
                     size=NULL, 
                     prob=NULL, 
                     mu=NULL) {
  
  # prob=NULL; mu=NULL; log = FALSE; infinite=10
  
  if (is.null(mu) + is.null(size) + is.null(prob) != 1) stop("Two values among mu, size and prob must be provided")
  
  m <- max(c(length(size), length(prob), length(mu)))
  if (!is.null(mu)) mu <- rep(mu, m)[1:m]
  if (!is.null(size)) size <- rep(size, m)[1:m]
  if (!is.null(prob)) prob <- rep(prob, m)[1:m]
  
  if (is.null(prob)) prob <- size/(size+mu)
  if (is.null(mu)) mu <- size/prob - size
  if (is.null(size))  size  <- (prob * mu) / (1 - prob)  
  #  if (length(prob)<length(size)) prob <- rep(prob, length(size))[1:length(size)]
  #  if (length(size)<length(prob)) size <- rep(size, length(prob))[1:length(prob)]
  
  m <- matrix(1:m, nrow=1)
  
  rl <- apply(m, MARGIN=2, function(i) rnbinom(n, size=size[i], prob=prob[i]))
  
  if (inherits(rl, "integer")) rl <- as.data.frame(matrix(rl, nrow=1))
  
  return(apply(rl, MARGIN=1, sum))
}

