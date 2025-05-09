#' Distribution of the Sum of Independent Negative Binomial Random Variables.
#' @title Distribution of the sum independent negative binomial random variables. 
#' @author Marc Girondot \email{marc.girondot@@gmail.com} and Jon Barry \email{jon.barry@@cefas.gov.uk}
#' @references Furman, E., 2007. On the convolution of the negative binomial random variables. Statistics & Probability Letters 77, 169-172.
#' @references Vellaisamy, P. & Upadhye, N.S. 2009. On the sums of compound negative binomial and gamma random variables. Journal of Applied Probability, 46, 272-283.
#' @references Girondot M, Barry J. 2023. Computation of the distribution of the sum of independent negative binomial random variables. Mathematical and Computational Applications 2023, 28, 63, doi:10.3390/mca28030063 
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
#' @param tol Tolerance for recurrence for Furman (2007) method. If NULL, will use a saddlepoint estimation.
#' @param method Can be Furman (default), Vellaisamy&Upadhye or exact, approximate.normal, approximate.negativebinomial, approximate.RandomObservations, or saddlepoint.
#' @param normalize If TRUE (default) will normalize the saddlepoint approximation estimate.
#' @param max.iter Number of maximum iterations for Furman method. Can be NULL.
#' @param mean Mean of the distribution for approximate.normal method. If NULL, the theoretical mean will be used.
#' @param sd Standard deviation of the distribution for approximate.normal method. If NULL, the theoretical sd will be used.
#' @param n.random Number of random numbers used to estimate parameters of distribution for approximate.RandomObservations method.
#' @param lower.tail logical; if TRUE (default), probabilities are P\[X <= x\], otherwise, P\[X > x\].
#' @param parallel logical; if FALSE (default), parallel computing is not used for Vellaisamy&Upadhye methods.
#' @param verbose Give more information on the method.
#' @description Distribution of the sum of random variable with negative binomial distributions.\cr
#' Technically the sum of random variable with negative binomial distributions is a convolution of negative 
#' binomial random variables.\cr
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
#' # By default, the Furman method with tol=1E-40 is used
#' dSnbinom(20, size=alpha, prob=p)
#' # Note the attribute is the dynamics of convergence of Pr(X=x)
#' attributes(dSnbinom(20, size=alpha, prob=p, verbose=TRUE))$Pk[, 1]
#' 
#' mutest <- c(0.01, 0.02, 0.03)
#' sizetest <- 2
#' x <- 20
#' # Exact probability
#' dSnbinom(x, size=sizetest, mu=mutest, method="vellaisamy&upadhye")
#' dSnbinom(x, size=sizetest, mu=mutest, method="vellaisamy&upadhye", log=TRUE)
#' # With Furman method and tol=1E-12, when probability 
#' # is very low, it will be biased
#' dSnbinom(x, size=sizetest, mu=mutest, method="Furman", tol=1E-12)
#' # The solution is to use a tolerance lower than the estimate
#' dSnbinom(x, size=sizetest, mu=mutest, method="Furman", tol=1E-45)
#' # Here the estimate used a first estimation by saddlepoint approximation
#' dSnbinom(x, size=sizetest, mu=mutest, method="Furman", tol=NULL)
#' # Or a huge number of iterations; but it is not the best solution
#' dSnbinom(x, size=sizetest, mu=mutest, method="Furman", 
#'          tol=1E-12, max.iter=10000)
#' # With the saddle point approximation method
#' dSnbinom(x, size=sizetest, mu=mutest, method="saddlepoint", log=FALSE)
#' dSnbinom(x, size=sizetest, mu=mutest, method="saddlepoint", log=TRUE)
#' 
#' # Another example
#' sizetest <- c(1, 1, 0.1)
#' mutest <- c(2, 1, 10)
#' x <- 5
#' (exact <- dSnbinom(x=x, size=sizetest, mu=mutest, method="Vellaisamy&Upadhye"))
#' (sp <- dSnbinom(x=x, size=sizetest, mu=mutest, method="saddlepoint"))
#' paste0("Saddlepoint approximation: Error of ", specify_decimal(100*abs(sp-exact)/exact, 2), "%")
#' (furman <- dSnbinom(x=x, size=sizetest, mu=mutest, method="Furman"))
#' paste0("Inversion of mgf: Error of ", specify_decimal(100*abs(furman-exact)/exact, 2), "%")
#' (na <- dSnbinom(x=x, size=sizetest, mu=mutest, method="approximate.normal")) 
#' paste0("Gaussian approximation: Error of ", specify_decimal(100*abs(na-exact)/exact, 2), "%")
#' (nb <- dSnbinom(x=x, size=sizetest, mu=mutest, method="approximate.negativebinomial"))
#' paste0("NB approximation: Error of ", specify_decimal(100*abs(nb-exact)/exact, 2), "%")
#' 
#' plot(0:20, dSnbinom(0:20, size=sizetest, mu=mutest, method="furman"), bty="n", type="h", 
#'      xlab="x", ylab="Density", ylim=c(0, 0.2), las=1)
#' points(x=0:20, y=dSnbinom(0:20, size=sizetest, mu=mutest, 
#'                            method="saddlepoint"), pch=1, col="blue")
#' points(x=0:20, y=dSnbinom(0:20, size=sizetest, mu=mutest, 
#'                            method="approximate.negativebinomial"), 
#'                            col="red")
#' points(x=0:20, y=dSnbinom(0:20, size=sizetest, mu=mutest, 
#'                            method="approximate.normal"), 
#'                            col="green")
#' 
#' # Test with a single distribution
#' dSnbinom(20, size=1, mu=20)
#' # when only one distribution is available, it is the same as dnbinom()
#' dnbinom(20, size=1, mu=20)
#' 
#' # If a parameter is supplied as only one value, it is supposed to be constant
#' dSnbinom(20, size=1, mu=c(14, 15, 10))
#' dSnbinom(20, size=c(1, 1, 1), mu=c(14, 15, 10))
#' 
#' # The functions are vectorized:
#' plot(0:200, dSnbinom(0:200, size=alpha, prob=p, method="furman"), bty="n", type="h", 
#'      xlab="x", ylab="Density")
#' points(0:200, dSnbinom(0:200, size=alpha, prob=p, method="saddlepoint"), 
#'      col="red", pch=3)
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
#' plot_add(0:(length(tabledistEmpirique)-1), tabledistEmpirique, type="l", col="red")
#' legend(x=200, y=0.02, legend=c("Empirical", "Theoretical"), 
#'    text.col=c("red", "black"), bty="n")
#' 
#' 
#' # Example from Vellaisamy, P. & Upadhye, N.S. (2009) - Table 1
#' # Note that computing time for k = 7 using exact method is very long
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
#' table1_Furman40 <- table1_Vellaisamy
#' table1_Furman40 <- table1_Vellaisamy
#' table1_FurmanAuto <- table1_Vellaisamy
#' table1_FurmanAuto_iter <- table1_Vellaisamy
#' table1_Vellaisamy_parallel <- table1_Vellaisamy
#' table1_Approximate_Normal <- table1_Vellaisamy
#' table1_saddlepoint <- table1_Vellaisamy
#' 
#' st_Furman3 <- rep(NA, length(k))
#' st_Furman6 <- rep(NA, length(k))
#' st_Furman9 <- rep(NA, length(k))
#' st_Furman12 <- rep(NA, length(k))
#' st_Furman40 <- rep(NA, length(k))
#' st_FurmanAuto <- rep(NA, length(k))
#' st_approximateObservations <- rep(NA, length(k))
#' st_Vellaisamy <- rep(NA, length(k))
#' st_Vellaisamy_parallel <- rep(NA, length(k))
#' st_Approximate_Normal <- rep(NA, length(k))
#' st_saddlepoint <- rep(NA, length(k))
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
#'     st_Furman40[which(n == k)] <- 
#'         system.time({
#'             table1_Furman40[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="Furman", tol=1E-40, log=FALSE, 
#'                                      verbose=FALSE)
#'         })[1]
#' 
#'     st_FurmanAuto[which(n == k)] <- 
#'         system.time({
#'             table1_FurmanAuto[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="Furman", tol=NULL, log=FALSE, 
#'                                      verbose=FALSE)
#'         })[1]
#'
#'     st_Approximate_Normal[which(n == k)] <- 
#'         system.time({
#'             table1_Approximate_Normal[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="approximate.normal", tol=1E-12, log=FALSE, 
#'                                      verbose=FALSE)
#'         })[1]
#'         st_saddlepoint[which(n == k)] <- 
#'         system.time({
#'             table1_saddlepoint[which(n == k), ] <- dSnbinom(x=x, prob=p, size=alpha, 
#'                                      method="saddlepoint", tol=1E-12, log=FALSE, 
#'                                      verbose=FALSE)
#'         })[1]
#'
#' for (xc in x) {
#'  essai <- dSnbinom(x=xc, prob=p, size=alpha, method="Furman", tol=NULL, log=FALSE, verbose=TRUE)
#'  table1_FurmanAuto_iter[which(n == k), which(xc == x)] <- nrow(attributes(essai)[[1]])
#' }
#' }
#' 
#' cbind(table1_Vellaisamy, st_Vellaisamy)
#' cbind(table1_Vellaisamy_parallel, st_Vellaisamy_parallel)
#' cbind(table1_Furman3, st_Furman3)
#' cbind(table1_Furman6, st_Furman6)
#' cbind(table1_Furman9, st_Furman9)
#' cbind(table1_Furman12, st_Furman12)
#' cbind(table1_Furman40, st_Furman40)
#' cbind(table1_FurmanAuto, st_FurmanAuto)
#' cbind(table1_approximateObservations, st_approximateObservations)
#' cbind(table1_Approximate_Normal, st_Approximate_Normal)
#' cbind(table1_saddlepoint, st_saddlepoint)
#' 
#' 
#' # Test of different methods
#' n <- 9
#' x <- 17
#' alpha <- 1:n
#' p <- (1:n)/10
#' 
#' # Parallel computing is not always performant
#' # Here it is very performant
#' system.time({print(dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, 
#'          verbose=TRUE, parallel=TRUE))})
#' system.time({print(dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", log=FALSE, 
#'          verbose=TRUE, parallel=FALSE))})
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
#' # Here parallel computing is 7 times faster (with a 8 cores computer) 
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
#' # Here parallel computing is 7 times faster (with a 8 cores computer) 
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
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, verbose=TRUE))
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Saddlepoint", log=FALSE, verbose=TRUE))
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="approximate.RandomObservations", 
#'          log=FALSE, verbose=TRUE))
#' 
#' # Test for criteria of convergence
#' Pr_exact <- dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", 
#'                    log=FALSE, verbose=TRUE)
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, 
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
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, 
#'                  verbose=TRUE)
#' Pr_saddlepoint <- dSnbinom(x=x, prob=p, size=alpha, method="saddlepoint", log=FALSE,  
#'                  verbose=TRUE)                  
#'                  
#' # pdf("figure.pdf", width=7, height=7, pointsize=14)                 
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
#' as.numeric(dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, verbose=TRUE))
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
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, 
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
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, 
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
#' Pr_ApproximateSaddlepoint <- dSnbinom(x=x, prob=p, size=alpha, 
#'                         method="saddlepoint", 
#'                         log=FALSE, 
#'                         verbose=TRUE)
#'              
#'              
#'              
#' layout(matrix(1:4, ncol=2, byrow=TRUE))
#' par(mar=c(3, 4.5, 1, 1))
#' alpha <- seq(from=10, to=100, length.out=3)
#' p <- seq(from=0.5, to=0.9, length.out=3)
#' 
#' p_nb <- dSnbinom(0:100, prob=p, size=alpha, method="vellaisamy&upadhye", verbose=TRUE)
#' p_Furman <- dSnbinom(0:100, prob=p, size=alpha, method="Furman", verbose=FALSE)
#' p_normal <- dSnbinom(0:100, prob=p, size=alpha, method="approximate.normal", verbose=TRUE)
#' p_aNB <- dSnbinom(0:100, prob=p, size=alpha, method="approximate.negativebinomial", verbose=TRUE)
#' p_SA <- dSnbinom(0:100, prob=p, size=alpha, method="saddlepoint", verbose=TRUE)
#' 
#' lab_PSnx <- bquote(italic("P(S"* "" [n] * "=x)"))
#' 
#' plot(1, 1, las=1, bty="n", col="grey", xlab="", 
#'        xlim=c(10, 70), ylim=c(0, 0.05), 
#'        ylab=lab_PSnx, type="n")
#' par(xpd=FALSE)
#' segments(x0=(0:100), x1=(0:100), 
#'          y0=0, y1=as.numeric(p_nb), col="black")
#'          
#' plot(x=p_nb, y=p_normal, pch=4, cex=0.5, las=1, bty="n", 
#'      xlab=bquote(italic("P" * "" [exact] * "(S"* "" [n] * "=x)")), 
#'      ylab=bquote(italic("P" * "" [approximate] * "(S"* "" [n] * "=x)")),  
#'      xlim=c(0, 0.05), ylim=c(0, 0.05))
#' points(x=p_nb, y=p_aNB, pch=5, cex=0.5)
#' points(x=p_nb, y=p_SA, pch=6, cex=0.5)
#' points(x=p_Furman, y=p_SA, pch=19, cex=0.5)
#' 
#' n <- 2
#' x <- 15
#' alpha <- 1:n
#' p <- (1:n)/10
#' 
#' p_nb <- dSnbinom(0:80, prob=p, size=alpha, method="vellaisamy&upadhye", verbose=TRUE)
#' p_Furman <- dSnbinom(0:80, prob=p, size=alpha, method="Furman", verbose=FALSE)
#' p_normal <- dSnbinom(0:80, prob=p, size=alpha, method="approximate.normal", verbose=TRUE)
#' p_aNB <- dSnbinom(0:80, prob=p, size=alpha, method="approximate.negativebinomial", verbose=TRUE)
#' p_SA <- dSnbinom(0:80, prob=p, size=alpha, method="saddlepoint", verbose=TRUE)
#' 
#' 
#' par(mar=c(4, 4.5, 1, 1))
#' plot(1, 1, las=1, bty="n", col="grey", xlab="x", 
#'        xlim=c(0, 60), ylim=c(0, 0.05), 
#'        ylab=lab_PSnx, type="n")
#' par(xpd=FALSE)
#' segments(x0=(0:80), x1=(0:80), 
#'          y0=0, y1=as.numeric(p_nb), col="black")
#'          
#' plot(x=p_nb, y=p_normal, pch=4, cex=0.5, las=1, bty="n", 
#'      xlab=bquote(italic("P" * "" [exact] * "(S"* "" [n] * "=x)")), 
#'      ylab=bquote(italic("P" * "" [approximate] * "(S"* "" [n] * "=x)")), 
#'      xlim=c(0, 0.05), ylim=c(0, 0.05))
#' points(x=p_nb, y=p_aNB, pch=5, cex=0.5)
#' points(x=p_nb, y=p_SA, pch=6, cex=0.5)
#' points(x=p_Furman, y=p_SA, pch=19, cex=0.5)
#'
#' # pdf("figure 1.pdf", width=7, height=7, pointsize=14)    
#' 
#' layout(1:2)
#' par(mar=c(3, 4, 1, 1))
#' alpha <- seq(from=10, to=100, length.out=3)
#' p <- seq(from=0.5, to=0.9, length.out=3)
#' 
#' p_nb <- dSnbinom(0:100, prob=p, size=alpha, method="vellaisamy&upadhye", verbose=TRUE)
#' p_Furman <- dSnbinom(0:100, prob=p, size=alpha, method="Furman", verbose=FALSE)
#' p_normal <- dSnbinom(0:100, prob=p, size=alpha, method="approximate.normal", verbose=TRUE)
#' p_aNB <- dSnbinom(0:100, prob=p, size=alpha, method="approximate.negativebinomial", verbose=TRUE)
#' p_SA <- dSnbinom(0:100, prob=p, size=alpha, method="saddlepoint", verbose=TRUE)
#' 
#' lab_PSnx <- bquote(italic("P(S"* "" [n] * "=x)"))
#' 
#' plot(1, 1, las=1, bty="n", col="grey", xlab="", 
#'        xlim=c(10, 70), ylim=c(0, 0.09), 
#'        ylab="", type="n", yaxt="n")
#' axis(2, at=seq(from=0, to=0.05, by=0.01), las=1)
#' mtext(lab_PSnx, side = 2, adj=0.3, line=3)
#' par(xpd=FALSE)
#' segments(x0=(0:100), x1=(0:100), 
#'          y0=0, y1=as.numeric(p_nb), col="black")
#' errr <- (abs((100*(p_Furman-p_nb)/p_nb)))/1000+0.055
#' errr <- ifelse(is.infinite(errr), NA, errr)
#' lines(x=(0:100), y=errr, lty=5, col="red", lwd=2)
#' errr <- (abs((100*(p_normal-p_nb)/p_nb)))/1000+0.055
#' lines(x=(0:100), y=errr, lty=2, col="blue", lwd=2)
#' errr <- (abs((100*(p_aNB-p_nb)/p_nb)))/1000+0.055
#' lines(x=(0:100), y=errr, lty=3, col="purple", lwd=2)
#' errr <- (abs((100*(p_SA-p_nb)/p_nb)))/1000+0.055
#' lines(x=(0:100), y=errr, lty=4, col="green", lwd=2)
#' axis(2, at=seq(from=0, to=40, by=10)/1000+0.055, las=1, 
#'      labels=as.character(seq(from=0, to=40, by=10)))
#' mtext("|% error|", side = 2, adj=0.9, line=3)
#' par(xpd=TRUE)
#' legend(x=30, y=0.1, legend=c("Inversion of mgf", "Saddlepoint", "Normal", "Negative binomial"), 
#'        lty=c(5, 4, 2, 3), bty="n", cex=0.8, col=c("red", "green", "blue", "purple"), lwd=2)
#' legend(x=10, y=0.05, legend=c("Exact"), lty=c(1), bty="n", cex=0.8)
#' par(xpd=TRUE)
#' text(x=ScalePreviousPlot(x = 0.95, y = 0.1)$x, 
#'      y=ScalePreviousPlot(x = 0.95, y = 0.1)$y, labels="A", cex=2)
#'      
#'      
#' # When normal approximation will fail
#' n <- 2
#' x <- 15
#' alpha <- 1:n
#' p <- (1:n)/10
#' 
#' p_nb <- dSnbinom(0:80, prob=p, size=alpha, method="vellaisamy&upadhye", verbose=TRUE)
#' p_Furman <- dSnbinom(0:80, prob=p, size=alpha, method="Furman", verbose=FALSE)
#' p_normal <- dSnbinom(0:80, prob=p, size=alpha, method="approximate.normal", verbose=TRUE)
#' p_aNB <- dSnbinom(0:80, prob=p, size=alpha, method="approximate.negativebinomial", verbose=TRUE)
#' p_SA <- dSnbinom(0:80, prob=p, size=alpha, method="saddlepoint", verbose=TRUE)
#' 
#' par(mar=c(4, 4, 1, 1))
#' plot(1, 1, las=1, bty="n", col="grey", xlab="x", 
#'        xlim=c(0, 60), ylim=c(0, 0.09), 
#'        ylab="", type="n", yaxt="n")
#' axis(2, at=seq(from=0, to=0.05, by=0.01), las=1)
#' mtext(lab_PSnx, side = 2, adj=0.3, line=3)
#' par(xpd=FALSE)
#' segments(x0=(0:80), x1=(0:80), 
#'          y0=0, y1=as.numeric(p_nb), col="black")
#' errr <- (abs((100*(p_Furman-p_nb)/p_nb)))/1000+0.055
#' errr <- ifelse(is.infinite(errr), NA, errr)
#' lines(x=(0:80), y=errr, lty=5, col="red", lwd=2)
#' errr <- (abs((100*(p_normal-p_nb)/p_nb)))/1000+0.055
#' lines(x=(0:80), y=errr, lty=2, col="blue", lwd=2)
#' errr <- (abs((100*(p_aNB-p_nb)/p_nb)))/1000+0.055
#' lines(x=(0:80), y=errr, lty=3, col="purple", lwd=2)
#' errr <- (abs((100*(p_SA-p_nb)/p_nb)))/1000+0.055
#' lines(x=(0:80), y=errr, lty=4, col="green", lwd=2)
#' axis(2, at=seq(from=0, to=40, by=10)/1000+0.055, las=1, 
#'      labels=as.character(seq(from=0, to=40, by=10)))
#' mtext("|% error|", side = 2, adj=0.9, line=3)
#' legend(x=30, y=0.055, 
#'        legend=c("Exact", "Inversion of mgf", "Saddlepoint", "Normal", "Negative binomial"), 
#'        lty=c(1, 5, 4, 2, 3), bty="n", cex=0.8, col=c("black", "red", "green", "blue", "purple"), 
#'        lwd=c(1, 2, 2, 2, 2))
#' par(xpd=TRUE)
#' text(x=ScalePreviousPlot(x = 0.95, y = 0.1)$x, 
#'      y=ScalePreviousPlot(x = 0.95, y = 0.1)$y, labels="B", cex=2)
#'      
#' 
#' # dev.off()
#' 
#' 
#' # Test for criteria of convergence
#' Pr_exact <- dSnbinom(x=x, prob=p, size=alpha, method="vellaisamy&upadhye", 
#'                    log=FALSE, verbose=TRUE)
#' Pr_Furman <- dSnbinom(x=x, prob=p, size=alpha, method="Furman", log=FALSE, 
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
#'                       method = "furman", verbose=FALSE, log=FALSE)})
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
#' Pr_Furman <- dSnbinom(x=x, size=ktest, mu=mutest, method="Furman", log=FALSE, 
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
#'    
#'    
#' # Test if saddlepoint approximation must be normalized
#' # Yes, it must be
#' n <- 7
#' alpha <- 1:n
#' p <- (1:n)/10
#' dSnbinom(x=10, prob=p, size=alpha, method="saddlepoint", log=FALSE,  
#'                  verbose=TRUE)
#' dSnbinom(x=10, prob=p, size=alpha, method="saddlepoint", log=FALSE,  
#'                  verbose=TRUE, normalize=FALSE)
#'                  
#' # Test for saddlepoint when x=0
#' n <- 7
#' alpha <- 1:n
#' p <- (1:n)/10
#' dSnbinom(x=0, prob=p, size=alpha, method="saddlepoint", log=FALSE,  
#'                  verbose=TRUE)
#' dSnbinom(x=1, prob=p, size=alpha, method="saddlepoint", log=FALSE,  
#'                  verbose=TRUE)
#' dSnbinom(x=c(0, 1), prob=p, size=alpha, method="saddlepoint", log=FALSE,  
#'                  verbose=TRUE)
#' dSnbinom(x=c(0, 1), prob=p, size=alpha, method="saddlepoint", log=FALSE,  
#'                  verbose=FALSE)
#'                  
#' # Test when prob are all the same
#' p <- rep(0.2, 7)
#' n <- 7
#' alpha <- 1:n
#' dSnbinom(x=0:10, prob=p, size=alpha, method="saddlepoint", log=FALSE,  
#'                  verbose=TRUE)
#' dSnbinom(x=0:10, prob=p, size=alpha, method="furman", log=FALSE,  
#'                  verbose=TRUE)
#' dSnbinom(x=0:10, prob=p, size=alpha, method="exact", log=FALSE,  
#'                  verbose=TRUE)
#'                  
#' # Test when n=1
#' p <- 0.2
#' n <- 1
#' alpha <- 1:n
#' dSnbinom(x=0:10, prob=p, size=alpha, method="saddlepoint", log=FALSE,  
#'                  verbose=TRUE)
#' dSnbinom(x=0:10, prob=p, size=alpha, method="furman", log=FALSE,  
#'                  verbose=TRUE)
#' dSnbinom(x=0:10, prob=p, size=alpha, method="exact", log=FALSE,  
#'                  verbose=TRUE)
#'                  
#' }
#' @export

#' @describeIn Snbinom Density for the sum of random variable with negative binomial distributions.
# #' @section Another section after function section:


dSnbinom <- function(x = stop("You must provide at least one x value")       , 
                     size = NULL                                             , 
                     prob = NULL                                             , 
                     mu = NULL                                               , 
                     log = FALSE                                             ,  
                     tol = NULL                                              , 
                     method="Furman"                                         ,
                     normalize=TRUE                                          ,
                     max.iter=NULL                                           ,
                     mean=NULL                                               ,
                     sd=NULL                                                 ,
                     n.random = 1E6                                          , 
                     parallel = FALSE                                        ,
                     verbose = FALSE                                         ) {
  
  method <- tolower(method)
  method <- match.arg(arg=method, choices = c("furman", 
                                              "vellaisamy&upadhye", "exact", 
                                              "approximate.randomobservations", 
                                              "approximate.normal", 
                                              "approximate.negativebinomial", 
                                              "saddlepoint"))
  
  # if (method == "convolution") method <- "furman"
  if (method == "exact") method <- "vellaisamy&upadhye"
  
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
    if (method == "vellaisamy&upadhye") {
      return(dnbinom(x, size=sum(size), prob=prob[1], log=log))
    } else {
      if (verbose) message("Exact method could be used because all prob parameters are the same.")
    }
  }
  
  #SaddlePoint approximation####
  
  if (method == "saddlepoint") {
    
    K = function(t, phi, mu) {
      core = phi*(log(phi) - log(phi+mu*(1-exp(t))))
      Kt = sum(core)
      Kt
    }
    
    Kdash1 = function(t, phi, mu) {
      ## 1st derivative
      num = phi*mu*exp(t)
      den = phi + mu*(1-exp(t))
      core = num / den
      Kdash1t = sum(core)
      Kdash1t
    }
    
    
    Kdash2 = function(t, phi, mu) {
      # 2nd derivative
      num = phi*mu*(phi+mu)*exp(t)
      den = (phi + mu*(1-exp(t)))^2
      core = num / den
      Kdash2t = sum(core)
      Kdash2t
    }
    
    saddlefun = function(t, phis, mus, x) {
      ## Used by optimise
      f = abs(Kdash1(t, phis, mus) - x)^2
      f
    }
    
    
    dsaddle = function(x, size, mu, tol, verbose=FALSE) {
      ## Saddle point approximation to density
      if (x == 0) {
        if (verbose) message("Exact method is used for x = 0.")
        p <- sum(dnbinom(x=0, size=size, mu=mu, log=TRUE))
        p <- exp(p)
        return(p)
      }
      
      upper.bound = min(log(size/mu + 1))
      lower.bound <- -100
      repp <- 0
      repeat {
        ## I'm not sure how to set the lower bound in the general case
        arse = optimise(saddlefun, phis=size, mus=mu, x=x, 
                        lower=lower.bound, upper=upper.bound, 
                        tol=tol)
        repp <- repp + 1
        if ((arse$minimum == lower.bound) & repp < 10)
          lowerbound <- lowerbound / 10
        else break
      }
      
      if (verbose & (repp == 10)) warning("The lower bound reached its minimum; use results with caution.")
      
      sx = arse$minimum
      
      ## Generate the density now that we have xs
      numfx = exp(K(sx, size, mu) - sx*x)
      denfx = sqrt(2*pi*Kdash2(sx, size, mu))
      fx = numfx / denfx
      fx
    }
    
    if (verbose) message("Saddlepoint approximation method")
    if (is.null(tol)) {
      tol <- 1E-10
    }
    
    if (normalize & any(x != 0)) {
      if (verbose) message(paste0("Tolerance for normalization=", as.character(tol)))
      # Prepare normalization
      mean <- sum(size*(1-prob)/prob)
      sd <- sqrt(sum(size*(1-prob)/prob^2))
      Max <- max(c(floor(mean+20*sd), x))+1
      dstot <- sapply(X = 0:Max, FUN=function(y) dsaddle(x = y, 
                                                         size = size, 
                                                         mu = mu, tol=tol, verbose=FALSE))
      Max <- Max + 1
      repeat {
        ds <- dsaddle(x = Max, size = size, mu = mu, tol=tol, verbose=FALSE)
        dstot <- c(dstot, ds)
        if (ds < tol) break
        Max <- Max + 1
        if (verbose) message(paste0("Normalization: P(Sn=", as.character(Max), ")=", as.character(ds)))
      }
      
      if (verbose) message(paste0("Sum for normalization = ", as.character(sum(dstot))))
      ds <- dstot[x+1]/sum(dstot)
      
    } else {
      ds <- sapply(X = x, FUN=function(y) dsaddle(x = y, 
                                                  size = size, 
                                                  mu = mu, tol=tol, verbose=verbose))
      if (verbose & any(x != 0)) warning("Saddlepoint approximation was not normalized. Use this output with caution.")
    }
    if (!log) return(ds) else return(log(ds))
  }
  
  #Furman - Convolution####
  
  if (method == "furman") {
    if (verbose) message("Furman (2007) method by inversion of moment generating function")
    if (is.null(tol)) {
      tol <- min(dSnbinom(x=x, size = size, mu=mu, log=FALSE, 
                          normalize=FALSE, 
                          tol=1E-10, method="saddlepoint"))*1E-10
    }
    
    alpha <- size
    p <- prob
    
    q <- 1-p
    p1 <- max(p)
    q1 <- 1-p1
    
    # R <- sum(log(((q*p1)/(q1*p))^(-alpha)))
    R <- sum((-alpha)*log(((q*p1)/(q1*p))))
    
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
      # xi <- c(xi, sum((alpha*(1-((q1*p)/(q*p1)))^i)/i))
      # xi <- c(xi, sum((alpha/i*((1-((q1*p)/(q*p1)))^i))/i))
      # 
      # (alpha*(1-((q1*p)/(q*p1)))^i)/i
      # exp(log((alpha*(1-((q1*p)/(q*p1)))^i)/i))
      # exp(log((alpha*(1-((q1*p)/(q*p1)))^i))-log(i))
      # 
      # exp((log(alpha)+i*log((1-((q1*p)/(q*p1)))))-log(i))
      # exp((log(alpha)+i*log((q*p1-q1*p)/(q*p1)))-log(i))
      # exp(log(alpha)+i*(log(q*p1-q1*p)-log(q*p1))-log(i))
      # New expression to prevent the use of ^
      xi <- c(xi, sum(exp(log(alpha)+i*(log(q*p1-q1*p)-log(q)-log(p1))-log(i))))
      
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
  
  #approximate.randomobservations####
  
  if (method == "approximate.randomobservations") {
    if (verbose) message("Approximate method with probabilities of observations.")
    test <- rSnbinom(n=n.random, size = size, mu=mu)
    return(sapply(X = x, FUN=function(y) sum(test==y)/n.random))
  }
  
  #approximate.negativebinomial####
  
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
  
  #approximate.normal####
  
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
  
  #Exact####
  
  if (verbose) message("Exact Vellaisamy & Upadhye (2009) method")
  if ((m > 7) & verbose) message("The Vellaisamy method with more than 7 summed distributions can be slow and produces out of memory error.")
  
  
  TS <- NULL
  
  for (xec in x) {
    
    # créer df avec xec balls dans m boxes
    
    Nballs <- xec
    # Number of boxes
    Nboxes <- m
    
    # The number of different ways to distribute n indistinguishable balls into
    # k distinguishable boxes is C(n+k-1,k-1).
    # nb<-choose(N+nbjour-1,nbjour-1)=dim(tb)[1]
    # divers<-matrix(rep(0, nbjour*nb), ncol=nbjour)
    
    # generate all possible positions of the boundaries
    xx <- combn(Nballs+Nboxes-1, Nboxes-1)
    # compute the number of balls in each box
    a <- cbind(0, diag(Nboxes)) - cbind(diag(Nboxes), 0)
    df <- t(a %*% rbind(0, xx, Nballs+Nboxes) - 1)
    
    if (verbose) {
      message(paste0(as.character(nrow(df)), " combinations of ", as.character(xec), " objects in ", as.character(m), " categories. The number of required iterations will be ", as.character(nrow(df)*m), "."))
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
                     prob=NULL, mu=NULL, lower.tail = TRUE, log.p = FALSE, tol=NULL, 
                     method="Furman", normalize=TRUE) {
  
  method <- tolower(method)
  method <- match.arg(arg=method, choices = c("furman", 
                                              "vellaisamy&upadhye", "exact", 
                                              "approximate.randomobservations", 
                                              "approximate.normal", 
                                              "approximate.negativebinomial", 
                                              "saddlepoint"))
  
  # if (method == "convolution") method <- "furman"
  if (method == "exact") method <- "vellaisamy&upadhye"
  
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
                  tol=tol, method = method, normalize=normalize)
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
                     tol=NULL, method="Furman") {
  
  # prob=NULL; mu=NULL; log = FALSE; infinite=10
  
  method <- tolower(method)
  method <- match.arg(arg=method, choices = c("furman", 
                                              "vellaisamy&upadhye", "exact", 
                                              "approximate.randomobservations", 
                                              "approximate.normal", 
                                              "approximate.negativebinomial", 
                                              "saddlepoint"))
  
 # if (method == "convolution") method <- "furman"
  if (method == "exact") method <- "vellaisamy&upadhye"
  
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

