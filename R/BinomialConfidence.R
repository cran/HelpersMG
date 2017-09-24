#' @title Confidence Intervals for Binomial Probabilities
#' @author Rollin Brant, Modified by Frank Harrell and Brad Biggerstaff 
#' @return a matrix or data.frame containing the computed intervals and, optionally, x and n.
#' @param x Vector containing the number of "successes" for binomial variates
#' @param n Vector containing the numbers of corresponding observations
#' @param alpha Probability of a type I error, so confidence coefficient = 1-alpha
#' @param method Character string specifing which method to use. The "all" method only works when x and n are length 1. The "exact" method uses the F distribution to compute exact (based on the binomial cdf) intervals; the "wilson" interval is score-test-based; and the "asymptotic" is the text-book, asymptotic normal interval. Following Agresti and Coull, the Wilson interval is to be preferred and so is the default.
#' @param include.x Logical flag to indicate whether x should be included in the returned matrix or data frame
#' @param include.n Logical flag to indicate whether n should be included in the returned matrix or data frame
#' @param return.df Logical flag to indicate that a data frame rather than a matrix be returned
#' @description Produces 1-alpha confidence intervals for binomial probabilities.\cr
#' @examples
#' \dontrun{
#' HelpersMG:::.BinomialConfidence(0:10,10,include.x=TRUE,include.n=TRUE)
#' HelpersMG:::.BinomialConfidence(46,50,method="all")
#' }
#' @export

.BinomialConfidence <- 
function (x, n, alpha = 0.05, method = c("wilson", "exact", "asymptotic", 
    "all"), include.x = FALSE, include.n = FALSE, return.df = FALSE) 
{
    method <- match.arg(method)
    bc <- function(x, n, alpha, method) {
        nu1 <- 2 * (n - x + 1)
        nu2 <- 2 * x
        ll <- if (x > 0) 
            x/(x + qf(1 - alpha/2, nu1, nu2) * (n - x + 1))
        else 0
        nu1p <- nu2 + 2
        nu2p <- nu1 - 2
        pp <- if (x < n) 
            qf(1 - alpha/2, nu1p, nu2p)
        else 1
        ul <- ((x + 1) * pp)/(n - x + (x + 1) * pp)
        zcrit <- -qnorm(alpha/2)
        z2 <- zcrit * zcrit
        p <- x/n
        cl <- (p + z2/2/n + c(-1, 1) * zcrit * sqrt((p * (1 - 
            p) + z2/4/n)/n))/(1 + z2/n)
        if (x == 1) 
            cl[1] <- -log(1 - alpha)/n
        if (x == (n - 1)) 
            cl[2] <- 1 + log(1 - alpha)/n
        asymp.lcl <- x/n - qnorm(1 - alpha/2) * sqrt(((x/n) * 
            (1 - x/n))/n)
        asymp.ucl <- x/n + qnorm(1 - alpha/2) * sqrt(((x/n) * 
            (1 - x/n))/n)
        res <- rbind(c(ll, ul), cl, c(asymp.lcl, asymp.ucl))
        res <- cbind(rep(x/n, 3), res)
        switch(method, wilson = res[2, ], exact = res[1, ], asymptotic = res[3, 
            ], all = res, res)
    }
    if ((length(x) != length(n)) & length(x) == 1) 
        x <- rep(x, length(n))
    if ((length(x) != length(n)) & length(n) == 1) 
        n <- rep(n, length(x))
    if ((length(x) > 1 | length(n) > 1) & method == "all") {
        method <- "wilson"
        warning("method=all will not work with vectors...setting method to wilson")
    }
    if (method == "all" & length(x) == 1 & length(n) == 1) {
        mat <- bc(x, n, alpha, method)
        dimnames(mat) <- list(c("Exact", "Wilson", "Asymptotic"), 
            c("PointEst", "Lower", "Upper"))
        if (include.n) 
            mat <- cbind(N = n, mat)
        if (include.x) 
            mat <- cbind(X = x, mat)
        if (return.df) 
            mat <- as.data.frame(mat)
        return(mat)
    }
    mat <- matrix(ncol = 3, nrow = length(x))
    for (i in 1:length(x)) mat[i, ] <- bc(x[i], n[i], alpha = alpha, 
        method = method)
    dimnames(mat) <- list(rep("", dim(mat)[1]), c("PointEst", 
        "Lower", "Upper"))
    if (include.n) 
        mat <- cbind(N = n, mat)
    if (include.x) 
        mat <- cbind(X = x, mat)
    if (return.df) 
        mat <- as.data.frame(mat, row.names = NULL)
    class(mat) <- "binconf"
    mat
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
