#' universalmclapply runs the function FUN on X using parallel computing
#' @title Run the function FUN on X using parallel computing
#' @author Marc Girondot
#' @return The results of the function FUN applied to X
#' @param X A vector (atomic or list) or an expressions vector. Other objects (including classed objects) will be coerced by as.list.
#' @param FUN The function to be applied to each element of X
#' @param ... Optional arguments to FUN
#' @param mc.preschedule if set to TRUE then the computation is first divided to (at most) as many jobs are there are cores and then the jobs are started, each job possibly covering more than one value. If set to FALSE then one job is forked for each value of X. The former is better for short computations or large number of values in X, the latter is better for jobs that have high variance of completion time and not too many values of X compared to mc.cores.
#' @param mc.cores The number of cores to use, i.e. at most how many child processes will be run simultaneously. 
#' @param clusterExport List of clusterExport parameters as list
#' @param clusterEvalQ List of clusterEvalQ parameters as list
#' @param forking If TRUE will use forking
#' @param progressbar If pbapply package is installed, show a progressbar
#' @description Return the results of the function FUN applied to X. It uses forking in unix system and not in windows system.
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' x <- 1:1000
#' funx <- function(y) {
#'   mint <- rep(NA, length(y))
#'   for (i in seq_along(y)) {
#'     k <- rnorm(runif(n = 1, 50, 50), mean=10, sd=2)
#'     mint[i] <- mean(k)
#'   }
#'   mint
#' }
#' # Note that parallel computing is not always the best solution !
#' (tp <- system.time({
#'    m <- lapply(X=x, FUN=funx)
#' }))
#' (tp <- system.time({
#'    m <- universalmclapply(X=x, FUN=funx, forking=FALSE)
#' }))
#' (tp <- system.time({
#'    m <- universalmclapply(X=x, FUN=funx, forking=TRUE)
#' }))
#' 
#' ### An example using clusterExport
#' # Here no error is generated because environment was exported
#' # However forking is not possible in windows and non parallel code is ran
#' pp <- runif(100)
#' x <- 1:100
#' funx1 <- function(y) {pp[y]*10}
#' u <- universalmclapply(x, FUN=funx1, forking=TRUE)
#' 
#' # Here an error is generated because environment was not exported when parLapplyLB is used
#' pp <- runif(100)
#' x <- 1:100
#' u <- universalmclapply(x, FUN=funx1, forking=FALSE)
#' 
#' # here no error is generated because the variable pp is exported
#' pp <- runif(100)
#' x <- 1:100
#' u <- universalmclapply(x, FUN=funx1, forking=FALSE, 
#'                        clusterExport=list(varlist=c("pp"), envir=environment()))
#'                        
#' ### An example using clusterEvalQ
#' asc("a") # asc() is a function from packages HelpersMG
#' funx2 <- function(y) {asc("a")*10}
#' # In unix, the loaded packages are visible from all cores
#' x <- 1:100
#' u <- universalmclapply(x, FUN=funx2, forking=TRUE)
#' # In windows, the loaded packages are not visible from all cores
#' x <- 1:100
#' u <- universalmclapply(x, FUN=funx2, forking=FALSE)
#' # In windows, the loaded packages are not visible from all cores
#' x <- 1:100
#' u <- universalmclapply(x, FUN=funx2, forking=FALSE, 
#' clusterEvalQ=list(expr=expression(library(HelpersMG)))
#' )
#' 
#' ### If package pbapply is available, progress bar can be shown
#' m <- universalmclapply(X=x, FUN=funx, forking=FALSE, progressbar=TRUE)
#' m <- universalmclapply(X=x, FUN=funx, forking=TRUE, progressbar=TRUE)
#' }
#' @export


#

universalmclapply <- function(X, FUN, ..., 
                              mc.cores=getOption("mc.cores", parallel::detectCores()), 
                              mc.preschedule = TRUE, 
                              clusterExport=list(), 
                              clusterEvalQ=list(), 
                              forking=ifelse(.Platform$OS.type=="windows", 
                                             FALSE, TRUE), 
                              progressbar=FALSE) {
  m <- NULL

if (forking) {
  if (progressbar & ("pbapply" %in% rownames(installed.packages()))) {
    m <- do.call(pbapply::pblapply, modifyList(list(...), list(X=X, FUN=FUN, cl=mc.cores)))
  } else {
    m <- do.call(parallel::mclapply, modifyList(list(...), list(X=X, FUN=FUN, 
                                                                mc.cores=mc.cores, 
                                        mc.preschedule = mc.preschedule)))
  }
} else {
  
  cl <- parallel::makeCluster(mc.cores)
  if (!identical(list(), clusterExport)) {
  do.call(parallel::clusterExport, args = modifyList(clusterExport, list(cl=cl)))
  }
  if (!identical(list(), clusterEvalQ)) {
  do.call(parallel::clusterEvalQ, args=modifyList(clusterEvalQ, list(cl=cl)))
  }
  if (progressbar & ("pbapply" %in% rownames(installed.packages()))) {
    tt <- try(expr = {
    m <- do.call(pbapply::pblapply, args=modifyList(list(...), list(X=X, FUN=FUN, cl=cl)))
    }, silent = FALSE)
  } else {
    tt <- try(expr = {
    m <- do.call(parallel::parLapplyLB, args=modifyList(list(...), list(X=X, fun=FUN, cl=cl)))
    }, silent = FALSE)
  }
  parallel::stopCluster(cl)
}
return(m)
}
