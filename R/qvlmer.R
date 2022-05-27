#' qvlmer is Quasi Variances for lmer Model Coefficients
#' @title Quasi Variances for lmer Model Coefficients
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A list of class qv.
#' @param object A object obtained using lmer from package lme4
#' @param factorname Either NULL, or a character vector of length 1
#' @param coef.indices Either NULL, or a numeric vector of length at least 3
#' @param dispersion An optional scalar multiplier for the covariance matrix, to cope with overdispersion for example
#' @param ... Other arguments to pass to qvcalc.default
#' @description Computes a set of quasi variances (and corresponding quasi standard errors) 
#' for estimated model coefficients relating to the levels of a categorical (i.e., factor) 
#' explanatory variable. For details of the method see Firth (2000), Firth (2003) or Firth 
#' and de Menezes (2004). Quasi variances generalize and improve the accuracy of “floating 
#' absolute risk” (Easton et al., 1991). This device for economical model summary was first 
#' suggested by Ridout (1989).\cr
#' Modified from qvcalc.lm() of packages qvcalc by David Firth, d.firth@@warwick.ac.uk
#' @references Easton, D. F, Peto, J. and Babiker, A. G. A. G. (1991) Floating absolute risk: an alternative to relative risk in survival and case-control analysis avoiding an arbitrary reference group. Statistics in Medicine 10, 1025–1035.
#' @references Firth, D. (2000) Quasi-variances in Xlisp-Stat and on the web. Journal of Statistical Software 5.4, 1–13. At http://www.jstatsoft.org
#' @references Firth, D. (2003) Overcoming the reference category problem in the presentation of statistical models. Sociological Methodology 33, 1–18.
#' @references Firth, D. and de Mezezes, R. X. (2004) Quasi-variances. Biometrika 91, 65–80.
#' @references McCullagh, P. and Nelder, J. A. (1989) Generalized Linear Models. London: Chapman and Hall.
#' @references Menezes, R. X. de (1999) More useful standard errors for group and factor effects in generalized linear models. D.Phil. Thesis, Department of Statistics, University of Oxford.
#' @references Ridout, M.S. (1989). Summarizing the results of fitting generalized linear models to data from designed experiments. In: Statistical Modelling: Proceedings of GLIM89 and the 4th International Workshop on Statistical Modelling held in Trento, Italy, July 17–21, 1989 (A. Decarli et al., eds.), pp 262–269. New York: Springer.
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' y <- rnorm(100)
#' G <- as.factor(sample(c("A", "B", "C", "D"), 100, replace = TRUE))
#' R <- as.factor(rep(1:25, 4))
#' library(lme4)
#' m <- lmer(y ~ x + G + (1 | R))
#' qvlmer(m, factorname="G")
#' }
#' @export

qvlmer <- function (object, factorname = NULL, coef.indices = NULL, dispersion = NULL, 
          ...) {
  
  coef.indices.saved <- coef.indices
  if (is.null(factorname) && is.null(coef.indices)) {
    stop("arguments 'factorname' and 'coef.indices' are both NULL")
  }
  
  qvcd <- function (object, factorname = NULL, coef.indices = NULL, labels = NULL, 
            dispersion = NULL, estimates = NULL, modelcall = NULL, ...) 
  {
    covmat <- object
    if (!is.null(labels)) 
      rownames(covmat) <- colnames(covmat) <- labels
    n <- dim(covmat)[1]
    if (n <= 2) 
      stop("qvcalc works only for factors with 3 or more levels")
    simple.contrasts <- function(n, levelnames = 1:n) {
      result <- list()
      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          result[[paste(levelnames[i], levelnames[j], sep = ",")]] <- c(i, 
                                                                        j)
        }
      }
      result
    }
    qvdesign <- function(n) {
      nrows <- choose(n, 2)
      m <- matrix(0, nrows, n)
      indices <- simple.contrasts(n)
      for (i in 1:nrows) {
        m[i, indices[[i]][1]] <- 1
        m[i, indices[[i]][2]] <- 1
      }
      m
    }
    level <- qvdesign(n)
    contrast.variance <- function(contrast, covmat) {
      if (!(is.matrix(covmat) && (dim(covmat)[1] == dim(covmat)[2]))) 
        stop("covmat must be a square matrix")
      n <- dim(covmat)[1]
      if (length(contrast) == n && sum(contrast) == 0) 
        return(as.vector(contrast %*% covmat %*% contrast))
      if (length(contrast) == 2 && all(contrast %in% 1:n)) {
        i <- contrast[1]
        j <- contrast[2]
        return(covmat[i, i] + covmat[j, j] - 2 * covmat[i, 
                                                        j])
      } else stop("invalid contrast")
    }
    simple.contrast.variances <- function(n, covmat) {
      if (!is.null(rownames(covmat))) 
        levelnames <- rownames(covmat) else levelnames <- 1:n
      sapply(simple.contrasts(n, levelnames), function(contrast) {
        contrast.variance(contrast, covmat)
      })
    }
    response <- simple.contrast.variances(n, covmat)
    if (any(response <= 0)) {
      stop("not all contrasts have positive variance")
    } else response <- log(response)
    expLinear <- structure(list(family = "expLinear", link = "exp", 
                                linkfun = function(mu) exp(mu), linkinv = function(eta) log(eta), 
                                variance = function(mu) rep(1, length(mu)), dev.resids = function(y, 
                                                                                                  mu, wt) wt * ((y - mu)^2), aic = function(y, n, mu, 
                                                                                                                                            wt, dev) sum(wt) * (log(dev/sum(wt) * 2 * pi) + 1) + 
                                  2, mu.eta = function(eta) 1/eta, initialize = expression({
                                    n <- rep(1, nobs)
                                    mustart <- y
                                  }), validmu = function(mu) TRUE), class = "family")
    model <- glm(response ~ 0 + level, family = expLinear)
    qv <- coef(model)
    NAs <- rep(NA, length(qv))
    if (!is.null(rownames(covmat))) names(qv) <- rownames(covmat)
    frame <- data.frame(estimate = NAs, SE = sqrt(diag(covmat)), 
                        quasiSE = sqrt(qv), quasiVar = qv, row.names = names(qv))
    if (!is.null(estimates)) frame$estimate <- estimates
    relerrs <- sqrt(exp(-residuals(model))) - 1
    names(relerrs) <- names(response)
    return(structure(list(covmat = covmat, qvframe = frame, dispersion = dispersion, 
                          relerrs = relerrs, factorname = factorname, coef.indices = coef.indices, 
                          modelcall = modelcall), class = "qv"))
  }
  
  
  if (is.null(coef.indices)) {
    term.index <- which(attr(terms(object), "term.labels") == 
                          factorname)
    modelmat <- model.matrix(object)
    
    has.coef <- colnames(modelmat) %in% names(fixef(object))
    
    coef.indices <- which(attr(modelmat, "assign")[has.coef] == 
                            term.index)
    xlevels <- list()
    
    for (r in 1:ncol(object@frame)) 
      if ((is.factor(object@frame[, r])) & (colnames(object@frame)[r] %in% names(attr(object@pp$X, "contrasts")))) {
        l <- list(levels(object@frame[, r]))
        names(l) <- colnames(object@frame)[r]
        xlevels <- c(xlevels, l)
      }
    
    if (length(xlevels[[factorname]]) == length(coef.indices)) {
      contmat <- diag(length(coef.indices))
    } else {
      # contmat <- diag(length(xlevels[[factorname]]))
      # rownames(contmat) <- colnames(contmat) <- xlevels[[factorname]]
      contmat <- eval(call(attr(object@pp$X,"contrasts")[[factorname]], xlevels[[factorname]]))
    }
    
    estimates <- contmat %*% fixef(object)[coef.indices]
    covmat <- vcov(object, dispersion = dispersion)
    covmat <- covmat[coef.indices, coef.indices, drop = FALSE]
    covmat <- contmat %*% covmat %*% t(contmat)
  } else {
    k <- length(coef.indices)
    refPos <- numeric(0)
    if (0 %in% coef.indices) {
      refPos <- which(coef.indices == 0)
      coef.indices <- coef.indices[-refPos]
    }
    covmat <- vcov(object, dispersion = dispersion)
    covmat <- covmat[coef.indices, coef.indices, drop = FALSE]
    estimates <- coef(object)[coef.indices]
    if (length(refPos) == 1) {
      if (length(estimates) != k) 
        estimates <- c(0, estimates)
      covmat <- cbind(0, rbind(0, covmat))
      names(estimates)[1] <- rownames(covmat)[1] <- colnames(covmat)[1] <- "(reference)"
      if (refPos != 1) {
        if (refPos == k) {
          perm <- c(2:k, 1)
        } else {
          perm <- c(2:refPos, 1, (refPos + 1):k)
        }
        estimates <- estimates[perm]
        covmat <- covmat[perm, perm, drop = FALSE]
      }
    }
  }
  
  # print(covmat)
  # print(dim(covmat))
  
  return(qvcd(object=as.matrix(covmat), factorname = factorname, coef.indices = coef.indices.saved, 
                        labels = NULL, dispersion = dispersion, estimates = estimates, 
                        modelcall = object@call, ...))
}
