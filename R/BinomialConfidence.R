.BinomialConfidence <- 
function (x, n, conf.level = 0.95, sides = c("two.sided", "left", "right"),
          method = c("wilson", "wald", "waldcc", "agresti-coull", "jeffreys",
                     "modified wilson", "wilsoncc","modified jeffreys",
                     "clopper-pearson", "arcsine", "logit", "witting", "pratt", 
                     "midp", "lik", "blaker"),
          rand = 123, tol = 1e-05, std_est = TRUE, silent=TRUE)
#   Arguments
# x	number of successes.
# n	number of trials.
# conf.level  confidence level, defaults to 0.95.
# sides	 a character string specifying the side of the confidence interval, must be one of "two.sided" (default), "left" or "right". You can specify just the initial letter. "left" would be analogue to a hypothesis of "greater" in a t.test.
# method	character string specifing which method to use; this can be one out of: "wald", "wilson", "wilsoncc", "agresti-coull", "jeffreys", "modified wilson", "modified jeffreys", "clopper-pearson", "arcsine", "logit", "witting", "pratt", "midp", "lik" and "blaker". Defaults to "wilson". Abbreviation of method is accepted. See details.
# rand	seed for random number generator; see details.
# tol	tolerance for method "blaker".
# std_est	logical, specifying if the standard point estimator for the proportion value x/n should be returned (TRUE, default) or the method-specific internally used alternative point estimate (FALSE).

  {
    if (missing(method)) 
      method <- "wilsoncc"
    if (!silent) message(paste0("Method ", method))
    if (missing(sides)) 
      sides <- "two.sided"
    Recycle <- function (...) 
    {
      lst <- list(...)
      maxdim <- max(lengths(lst))
      res <- lapply(lst, rep, length.out = maxdim)
      attr(res, "maxdim") <- maxdim
      return(res)
    }
    
    iBinomCI <- function(x, n, conf.level = 0.95, sides = c("two.sided", 
                                                            "left", "right"), method = c("wilson", "wilsoncc", "wald", 
                                                                                         "waldcc", "agresti-coull", "jeffreys", "modified wilson", 
                                                                                         "modified jeffreys", "clopper-pearson", "arcsine", "logit", 
                                                                                         "witting", "pratt", "midp", "lik", "blaker"), rand = 123, 
                         tol = 1e-05, std_est = TRUE) {
      if (length(x) != 1) 
        stop("'x' has to be of length 1 (number of successes)")
      if (length(n) != 1) 
        stop("'n' has to be of length 1 (number of trials)")
      if (length(conf.level) != 1) 
        stop("'conf.level' has to be of length 1 (confidence level)")
      if (conf.level < 0.5 | conf.level > 1) 
        stop("'conf.level' has to be in [0.5, 1]")
      method <- match.arg(arg = method, choices = c("wilson", 
                                                    "wald", "waldcc", "wilsoncc", "agresti-coull", "jeffreys", 
                                                    "modified wilson", "modified jeffreys", "clopper-pearson", 
                                                    "arcsine", "logit", "witting", "pratt", "midp", "lik", 
                                                    "blaker"))
      sides <- match.arg(sides, choices = c("two.sided", "left", 
                                            "right"), several.ok = FALSE)
      if (sides != "two.sided") 
        conf.level <- 1 - 2 * (1 - conf.level)
      alpha <- 1 - conf.level
      kappa <- qnorm(1 - alpha/2)
      p.hat <- x/n
      q.hat <- 1 - p.hat
      est <- p.hat
      switch(method, wald = {
        term2 <- kappa * sqrt(p.hat * q.hat)/sqrt(n)
        CI.lower <- max(0, p.hat - term2)
        CI.upper <- min(1, p.hat + term2)
      }, waldcc = {
        term2 <- kappa * sqrt(p.hat * q.hat)/sqrt(n)
        term2 <- term2 + 1/(2 * n)
        CI.lower <- max(0, p.hat - term2)
        CI.upper <- min(1, p.hat + term2)
      }, wilson = {
        if (!std_est) {
          x.tilde <- x + kappa^2/2
          n.tilde <- n + kappa^2
          p.tilde <- x.tilde/n.tilde
          est <- p.tilde
        }
        term1 <- (x + kappa^2/2)/(n + kappa^2)
        term2 <- kappa * sqrt(n)/(n + kappa^2) * sqrt(p.hat * 
                                                        q.hat + kappa^2/(4 * n))
        CI.lower <- max(0, term1 - term2)
        CI.upper <- min(1, term1 + term2)
      }, wilsoncc = {
        if (!std_est) {
          x.tilde <- x + kappa^2/2
          n.tilde <- n + kappa^2
          p.tilde <- x.tilde/n.tilde
          est <- p.tilde
        }
        lci <- (2 * x + kappa^2 - 1 - kappa * sqrt(kappa^2 - 
                                                     2 - 1/n + 4 * p.hat * (n * q.hat + 1)))/(2 * 
                                                                                                (n + kappa^2))
        uci <- (2 * x + kappa^2 + 1 + kappa * sqrt(kappa^2 + 
                                                     2 - 1/n + 4 * p.hat * (n * q.hat - 1)))/(2 * 
                                                                                                (n + kappa^2))
        CI.lower <- max(0, ifelse(p.hat == 0, 0, lci))
        CI.upper <- min(1, ifelse(p.hat == 1, 1, uci))
      }, `agresti-coull` = {
        x.tilde <- x + kappa^2/2
        n.tilde <- n + kappa^2
        p.tilde <- x.tilde/n.tilde
        q.tilde <- 1 - p.tilde
        if (!std_est) est <- p.tilde
        term2 <- kappa * sqrt(p.tilde * q.tilde)/sqrt(n.tilde)
        CI.lower <- max(0, p.tilde - term2)
        CI.upper <- min(1, p.tilde + term2)
      }, jeffreys = {
        if (x == 0) CI.lower <- 0 else CI.lower <- qbeta(alpha/2, 
                                                         x + 0.5, n - x + 0.5)
        if (x == n) CI.upper <- 1 else CI.upper <- qbeta(1 - 
                                                           alpha/2, x + 0.5, n - x + 0.5)
      }, `modified wilson` = {
        if (!std_est) {
          x.tilde <- x + kappa^2/2
          n.tilde <- n + kappa^2
          p.tilde <- x.tilde/n.tilde
          est <- p.tilde
        }
        term1 <- (x + kappa^2/2)/(n + kappa^2)
        term2 <- kappa * sqrt(n)/(n + kappa^2) * sqrt(p.hat * 
                                                        q.hat + kappa^2/(4 * n))
        if ((n <= 50 & x %in% c(1, 2)) | (n >= 51 & x %in% 
                                          c(1:3))) CI.lower <- 0.5 * qchisq(alpha, 2 * 
                                                                              x)/n else CI.lower <- max(0, term1 - term2)
        if ((n <= 50 & x %in% c(n - 1, n - 2)) | (n >= 51 & 
                                                  x %in% c(n - (1:3)))) CI.upper <- 1 - 0.5 * qchisq(alpha, 
                                                                                                     2 * (n - x))/n else CI.upper <- min(1, term1 + 
                                                                                                                                           term2)
      }, `modified jeffreys` = {
        if (x == n) CI.lower <- (alpha/2)^(1/n) else {
          if (x <= 1) CI.lower <- 0 else CI.lower <- qbeta(alpha/2, 
                                                           x + 0.5, n - x + 0.5)
        }
        if (x == 0) CI.upper <- 1 - (alpha/2)^(1/n) else {
          if (x >= n - 1) CI.upper <- 1 else CI.upper <- qbeta(1 - 
                                                                 alpha/2, x + 0.5, n - x + 0.5)
        }
      }, `clopper-pearson` = {
        CI.lower <- qbeta(alpha/2, x, n - x + 1)
        CI.upper <- qbeta(1 - alpha/2, x + 1, n - x)
      }, arcsine = {
        p.tilde <- (x + 0.375)/(n + 0.75)
        if (!std_est) est <- p.tilde
        CI.lower <- sin(asin(sqrt(p.tilde)) - 0.5 * kappa/sqrt(n))^2
        CI.upper <- sin(asin(sqrt(p.tilde)) + 0.5 * kappa/sqrt(n))^2
      }, logit = {
        lambda.hat <- log(x/(n - x))
        V.hat <- n/(x * (n - x))
        lambda.lower <- lambda.hat - kappa * sqrt(V.hat)
        lambda.upper <- lambda.hat + kappa * sqrt(V.hat)
        CI.lower <- exp(lambda.lower)/(1 + exp(lambda.lower))
        CI.upper <- exp(lambda.upper)/(1 + exp(lambda.upper))
      }, witting = {
        set.seed(rand)
        x.tilde <- x + runif(1, min = 0, max = 1)
        pbinom.abscont <- function(q, size, prob) {
          v <- trunc(q)
          term1 <- pbinom(v - 1, size = size, prob = prob)
          term2 <- (q - v) * dbinom(v, size = size, prob = prob)
          return(term1 + term2)
        }
        qbinom.abscont <- function(p, size, x) {
          fun <- function(prob, size, x, p) {
            pbinom.abscont(x, size, prob) - p
          }
          uniroot(fun, interval = c(0, 1), size = size, 
                  x = x, p = p)$root
        }
        CI.lower <- qbinom.abscont(1 - alpha, size = n, x = x.tilde)
        CI.upper <- qbinom.abscont(alpha, size = n, x = x.tilde)
      }, pratt = {
        if (x == 0) {
          CI.lower <- 0
          CI.upper <- 1 - alpha^(1/n)
        } else if (x == 1) {
          CI.lower <- 1 - (1 - alpha/2)^(1/n)
          CI.upper <- 1 - (alpha/2)^(1/n)
        } else if (x == (n - 1)) {
          CI.lower <- (alpha/2)^(1/n)
          CI.upper <- (1 - alpha/2)^(1/n)
        } else if (x == n) {
          CI.lower <- alpha^(1/n)
          CI.upper <- 1
        } else {
          z <- qnorm(1 - alpha/2)
          A <- ((x + 1)/(n - x))^2
          B <- 81 * (x + 1) * (n - x) - 9 * n - 8
          C <- (0 - 3) * z * sqrt(9 * (x + 1) * (n - x) * 
                                    (9 * n + 5 - z^2) + n + 1)
          D <- 81 * (x + 1)^2 - 9 * (x + 1) * (2 + z^2) + 
            1
          E <- 1 + A * ((B + C)/D)^3
          CI.upper <- 1/E
          A <- (x/(n - x - 1))^2
          B <- 81 * x * (n - x - 1) - 9 * n - 8
          C <- 3 * z * sqrt(9 * x * (n - x - 1) * (9 * 
                                                     n + 5 - z^2) + n + 1)
          D <- 81 * x^2 - 9 * x * (2 + z^2) + 1
          E <- 1 + A * ((B + C)/D)^3
          CI.lower <- 1/E
        }
      }, midp = {
        f.low <- function(pi, x, n) {
          1/2 * dbinom(x, size = n, prob = pi) + pbinom(x, 
                                                        size = n, prob = pi, lower.tail = FALSE) - 
            (1 - conf.level)/2
        }
        f.up <- function(pi, x, n) {
          1/2 * dbinom(x, size = n, prob = pi) + pbinom(x - 
                                                          1, size = n, prob = pi) - (1 - conf.level)/2
        }
        CI.lower <- 0
        CI.upper <- 1
        if (x != 0) {
          CI.lower <- uniroot(f.low, interval = c(0, p.hat), 
                              x = x, n = n)$root
        }
        if (x != n) {
          CI.upper <- uniroot(f.up, interval = c(p.hat, 
                                                 1), x = x, n = n)$root
        }
      }, lik = {
        CI.lower <- 0
        CI.upper <- 1
        z <- qnorm(1 - alpha * 0.5)
        tol = .Machine$double.eps^0.5
        BinDev <- function(y, x, mu, wt, bound = 0, tol = .Machine$double.eps^0.5, 
                           ...) {
          ll.y <- ifelse(y %in% c(0, 1), 0, dbinom(x, wt, 
                                                   y, log = TRUE))
          ll.mu <- ifelse(mu %in% c(0, 1), 0, dbinom(x, 
                                                     wt, mu, log = TRUE))
          res <- ifelse(abs(y - mu) < tol, 0, sign(y - 
                                                     mu) * sqrt(-2 * (ll.y - ll.mu)))
          return(res - bound)
        }
        if (x != 0 && tol < p.hat) {
          CI.lower <- if (BinDev(tol, x, p.hat, n, -z, 
                                 tol) <= 0) {
            uniroot(f = BinDev, interval = c(tol, if (p.hat < 
                                                      tol || p.hat == 1) 1 - tol else p.hat), bound = -z, 
                    x = x, mu = p.hat, wt = n)$root
          }
        }
        if (x != n && p.hat < (1 - tol)) {
          CI.upper <- if (BinDev(y = 1 - tol, x = x, mu = ifelse(p.hat > 
                                                                 1 - tol, tol, p.hat), wt = n, bound = z, tol = tol) < 
                          0) {
            CI.lower <- if (BinDev(tol, x, if (p.hat < 
                                               tol || p.hat == 1) 1 - tol else p.hat, n, 
                                   -z, tol) <= 0) {
              uniroot(f = BinDev, interval = c(tol, p.hat), 
                      bound = -z, x = x, mu = p.hat, wt = n)$root
            }
          } else {
            uniroot(f = BinDev, interval = c(if (p.hat > 
                                                 1 - tol) tol else p.hat, 1 - tol), bound = z, 
                    x = x, mu = p.hat, wt = n)$root
          }
        }
      }, blaker = {
        acceptbin <- function(x, n, p) {
          p1 <- 1 - pbinom(x - 1, n, p)
          p2 <- pbinom(x, n, p)
          a1 <- p1 + pbinom(qbinom(p1, n, p) - 1, n, p)
          a2 <- p2 + 1 - pbinom(qbinom(1 - p2, n, p), n, 
                                p)
          return(min(a1, a2))
        }
        CI.lower <- 0
        CI.upper <- 1
        if (x != 0) {
          CI.lower <- qbeta((1 - conf.level)/2, x, n - 
                              x + 1)
          while (acceptbin(x, n, CI.lower + tol) < (1 - 
                                                    conf.level)) CI.lower = CI.lower + tol
        }
        if (x != n) {
          CI.upper <- qbeta(1 - (1 - conf.level)/2, x + 
                              1, n - x)
          while (acceptbin(x, n, CI.upper - tol) < (1 - 
                                                    conf.level)) CI.upper <- CI.upper - tol
        }
      })
      ci <- c(est = est, lwr.ci = max(0, CI.lower), upr.ci = min(1, 
                                                                 CI.upper))
      if (sides == "left") 
        ci[3] <- 1
      else if (sides == "right") 
        ci[2] <- 0
      return(ci)
    }
    lst <- list(x = x, n = n, conf.level = conf.level, sides = sides, 
                method = method, rand = rand, std_est = std_est)
    maxdim <- max(unlist(lapply(lst, length)))
    lgp <- lapply(lst, rep, length.out = maxdim)
    lgn <- Recycle(x = if (is.null(names(x))) 
      paste("x", seq_along(x), sep = ".")
      else names(x), n = if (is.null(names(n))) 
        paste("n", seq_along(n), sep = ".")
      else names(n), conf.level = conf.level, sides = sides, method = method, 
      std_est = std_est)
    xn <- apply(as.data.frame(lgn[sapply(lgn, function(x) length(unique(x)) != 
                                           1)]), 1, paste, collapse = ":")
    res <- t(sapply(1:maxdim, function(i) iBinomCI(x = lgp$x[i], 
                                                   n = lgp$n[i], conf.level = lgp$conf.level[i], sides = lgp$sides[i], 
                                                   method = lgp$method[i], rand = lgp$rand[i], std_est = lgp$std_est[i])))
    # colnames(res)[1] <- c("est")
    rownames(res) <- rep("", nrow(res))
    colnames(res) <- c("PointEst", "Lower", "Upper")
    return(res)
  }


.Arrows <- function (x0, y0, x1, y1, code = 2, arr.length = 0.4, arr.width = arr.length/2, 
    arr.adj = 1, arr.type = "curved", segment = TRUE, col = "black", 
    lcol = col, lty = 1, arr.col = lcol, lwd = 1, arr.lwd = lwd, 
    ...) 
{
    if (arr.type == "simple") {
        arrows(x0, y0, x1, y1, code = code, length = arr.length/2.54, 
            lty = lty, col = col, lwd = lwd, ...)
        return()
    }
    if (arr.type == "T") {
        arrows(x0, y0, x1, y1, code = code, length = arr.length/(2 * 
            2.54), lty = lty, angle = 90, col = col, lwd = lwd, 
            ...)
        return()
    }
    if (segment) 
        segments(x0, y0, x1, y1, col = lcol, lty = lty, lwd = lwd, 
            ...)
    user <- par("usr")
    pin <- par("pin")
    pin <- pin/max(pin)
    sy <- (user[4] - user[3])/pin[2]
    sx <- (user[2] - user[1])/pin[1]
    angle <- atan((y1 - y0)/(x1 - x0) * sx/sy)/pi * 180
    angle[is.nan(angle)] <- 0
    angle[x1 < x0] <- 180 + angle[x1 < x0]
    xx <- x1
    yy <- y1
    if (code == 3) 
      getFromNamespace(".Arrowhead", ns="HelpersMG")(x0 = xx, y0 = yy, angle = angle, lcol = lcol, 
            arr.col = arr.col, arr.adj = arr.adj, lty = lty, 
            arr.length = arr.length, arr.width = arr.width, arr.type = arr.type, 
            arr.lwd = arr.lwd)
    if (code != 2) {
        angle <- 180 + angle
        xx <- x0
        yy <- y0
    }
    getFromNamespace(".Arrowhead", ns="HelpersMG")(x0 = xx, y0 = yy, angle = angle, lcol = lcol, arr.col = arr.col, 
        arr.adj = arr.adj, lty = lty, arr.length = arr.length, 
        arr.width = arr.width, arr.type = arr.type, arr.lwd = arr.lwd)
}

.Arrowhead <- function (x0, y0, angle = 0, arr.length = 0.4, arr.width = arr.length/2, 
    arr.adj = 0.5, arr.type = "curved", lcol = "black", lty = 1, 
    arr.col = lcol, arr.lwd = 2, npoint = 5) 
{
  
    if (arr.type == "curved") {
        rad <- 0.7
        len <- 0.25 * pi
        mid <- c(0, rad)
        x <- seq(1.5 * pi + len, 1.5 * pi, length.out = npoint)
        rr <- cbind(mid[1] - rad * cos(x), mid[2] + rad * sin(x))
        mid <- c(0, -rad)
        x <- rev(x)
        rr <- rbind(rr, cbind(mid[1] - rad * cos(x), mid[2] - 
            rad * sin(x)))
        mid <- c(rr[nrow(rr), 1], 0)
        rd <- rr[1, 2]
        x <- seq(pi/2, 3 * pi/2, length.out = 3 * npoint)
        rr <- rbind(rr, cbind(mid[1] - rd * 0.25 * cos(x), mid[2] - 
            rd * sin(x)))
        rr[, 1] <- rr[, 1] * 2.6
        rr[, 2] <- rr[, 2] * 3.45
    }
    else if (arr.type == "triangle") {
        x <- c(-0.2, 0, -0.2)
        y <- c(-0.1, 0, 0.1)
        rr <- 6.22 * cbind(x, y)
    }
    else if (arr.type %in% c("circle", "ellipse")) {
        if (arr.type == "circle") 
            arr.width = arr.length
        rad <- 0.1
        mid <- c(-rad, 0)
        x <- seq(0, 2 * pi, length.out = 15 * npoint)
        rr <- 6.22 * cbind(mid[1] + rad * sin(x), mid[2] + rad * 
            cos(x))
    }
    if (arr.adj == 0.5) 
        rr[, 1] <- rr[, 1] - min(rr[, 1])/2
    if (arr.adj == 0) 
        rr[, 1] <- rr[, 1] - min(rr[, 1])
    user <- par("usr")
    pcm <- par("pin") * 2.54
    sy <- (user[4] - user[3])/pcm[2]
    sx <- (user[2] - user[1])/pcm[1]
    nr <- max(length(x0), length(y0), length(angle), length(arr.length), 
        length(arr.width), length(lcol), length(lty), length(arr.col))
    if (nr > 1) {
        x0 <- rep(x0, length.out = nr)
        y0 <- rep(y0, length.out = nr)
        angle <- rep(angle, length.out = nr)
        arr.length <- rep(arr.length, length.out = nr)
        arr.width <- rep(arr.width, length.out = nr)
        lcol <- rep(lcol, length.out = nr)
        lty <- rep(lty, length.out = nr)
        arr.col <- rep(arr.col, length.out = nr)
    }
    RR <- rr
    for (i in 1:nr) {
        dx <- rr[, 1] * arr.length[i]
        dy <- rr[, 2] * arr.width[i]
        angpi <- angle[i]/180 * pi
        cosa <- cos(angpi)
        sina <- sin(angpi)
        RR[, 1] <- cosa * dx - sina * dy
        RR[, 2] <- sina * dx + cosa * dy
        RR[, 1] <- x0[i] + RR[, 1] * sx
        RR[, 2] <- y0[i] + RR[, 2] * sy
        polygon(RR, col = arr.col[i], border = lcol[i], lty = lty[i], 
            lwd = arr.lwd)
    }
}
