#' IC_clean_data checks and corrects the dataframe to be used with IC_threshold_matrix
#' @title Clean the dataframe before to be used with IC_threshold_matrix
#' @author Marc Girondot \email{marc.girondot@@gmail.com}
#' @return A dataframe
#' @param data The data.frame to be cleaned
#' @param use an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param method a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param variable.retain a vector with the name of columns to keep
#' @param test.partial.correlation should the partial correlations be tested ?
#' @param progress Show a progress bar
#' @param debug if TRUE, information about progression of cleaning are shown
#' @description This function must be used if missing values are present in the dataset.\cr
#' It ensures that all correlations and partial correlations can be calculated. 
#' The columns of the dataframe are removed one per one until all can be calculated without error. 
#' It is possible to say that one or more columns must be retained because they are of particular importance in the analysis. 
#' The use and method parameters are used by cor() function. The function uses by default a parallel computing in Unix or MacOSX systems. 
#' If progress is TRUE and the package pbmcapply is present, a progress bar is displayed. If debug is TRUE, some informations are shown during the process.
#' \code{https://fr.wikipedia.org/wiki/Iconographie_des_corrélations}
#' @family Iconography of correlations
#' @references Lesty, M., 1999. Une nouvelle approche dans le choix des régresseurs de la régression multiple en présence d’interactions et de colinéarités. Revue de Modulad 22, 41-77.
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' # based on https://fr.wikipedia.org/wiki/Iconographie_des_corrélations
#' es <- structure(list(Student = c("e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8"), 
#'                      Mass = c(52, 59, 55, 58, 66, 62, 63, 69), 
#'                      Age = c(12, 12.5, 13, 14.5, 15.5, 16, 17, 18), 
#'                      Assiduity = c(12, 9, 15, 5, 11, 15, 12, 9), 
#'                      Note = c(5, 5, 9, 5, 13.5, 18, 18, 18)), 
#'                      row.names = c(NA, -8L), class = "data.frame")
#' es
#' 
#' df_clean <- IC_clean_data(es, debug = TRUE)
#' cor_matrix <- IC_threshold_matrix(data=df_clean, threshold = NULL, progress=FALSE)
#' cor_threshold <- IC_threshold_matrix(data=df_clean, threshold = 0.3)
#' plot(cor_threshold, show.legend.strength=FALSE, show.legend.direction = FALSE)
#' cor_threshold_Note <- IC_correlation_simplify(matrix=cor_threshold, variable="Note")
#' plot(cor_threshold_Note, show.legend.strength=FALSE, show.legend.direction = FALSE)
#' 
#' cor_threshold <- IC_threshold_matrix(data=df_clean, threshold = 0.6)
#' plot(cor_threshold, 
#' layout=matrix(data=c(53, 53, 55, 55, 
#'                      55, 53, 55, 53), ncol=2, byrow=FALSE), 
#' show.legend.direction = FALSE,
#' show.legend.strength = FALSE, xlim=c(-2, 2), ylim=c(-2, 2))
#' }
#' @export


IC_clean_data <- function(data=stop("A dataframe object is required"), 
                       use = c("pairwise.complete.obs", "everything", 
                                     "all.obs", "complete.obs", "na.or.complete"), 
                       method=c("pearson", "kendall", "spearman"), 
                       variable.retain=NULL,
                       test.partial.correlation=TRUE, 
                       progress=TRUE, debug=FALSE) {
  
  # data=NULL
  # use="pairwise.complete.obs"
  # method="pearson"
  # variable.retain=NULL
  # test.partial.correlation=TRUE
  # progress=TRUE
  # debug=FALSE
  
  method <- method[1]
  use <- use[1]
  
  if (!(is.element("ppcor", installed.packages()[,1]))) {
    ginv <- function (X, tol = sqrt(.Machine$double.eps)) 
    {
      if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
      if (!is.matrix(X)) 
        X <- as.matrix(X)
      Xsvd <- svd(X)
      if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
      Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
      if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
      else if (!any(Positive)) 
        array(0, dim(X)[2L:1L])
      else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                                                   t(Xsvd$u[, Positive, drop = FALSE]))
    }
    
    pcor <- function (x, method = c("pearson", "kendall", "spearman")) 
    {
      method <- match.arg(method)
      if (is.data.frame(x)) 
        x <- as.matrix(x)
      if (!is.matrix(x)) 
        stop("supply a matrix-like 'x'")
      if (!(is.numeric(x) || is.logical(x))) 
        stop("'x' must be numeric")
      stopifnot(is.atomic(x))
      n <- dim(x)[1]
      gp <- dim(x)[2] - 2
      cvx <- cov(x, method = method)
      if (det(cvx) < .Machine$double.eps) {
        warning("The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero.")
        icvx <- ginv(cvx)
      }
      else icvx <- solve(cvx)
      pcor <- -cov2cor(icvx)
      diag(pcor) <- 1
      if (method == "kendall") {
        statistic <- pcor/sqrt(2 * (2 * (n - gp) + 5)/(9 * (n - 
                                                              gp) * (n - 1 - gp)))
        p.value <- 2 * pnorm(-abs(statistic))
      }
      else {
        statistic <- pcor * sqrt((n - 2 - gp)/(1 - pcor^2))
        p.value <- 2 * pt(-abs(statistic), (n - 2 - gp))
      }
      diag(statistic) <- 0
      diag(p.value) <- 0
      list(estimate = pcor, p.value = p.value, statistic = statistic, 
           n = n, gp = gp, method = method)
    }
  } else {
    pcor <- getFromNamespace("pcor", ns="ppcor")
  }


  zc_y <- data
  # J'enlève les colonnes qui n'ont que des NA
   if (debug) cat(paste0("The data have ", ncol(zc_y), " variables\n"))
   if (debug) cat("I remove the columns with only NA\n")
  if (debug) {cat(colnames(zc_y)[apply(X = zc_y, MARGIN = 2, FUN = function(x) all(is.na(x)))]); cat("\n")}
  zc_y <- zc_y[, !apply(X = zc_y, MARGIN = 2, FUN = function(x) all(is.na(x)))]
   if (debug) cat(paste0("The data have ", ncol(zc_y), " variables\n"))
  
  # J'enlève les colonnes qui ne présente pas de variabilité
  if (debug) cat("I remove the columns with no variability\n")
  if (debug) {cat(colnames(zc_y)[apply(X = zc_y, MARGIN = 2, FUN = function(x) all(na.omit(x)==na.omit(x)[1]))]); cat("\n")}
  zc_y <- zc_y[, !apply(X = zc_y, MARGIN = 2, FUN = function(x) all(na.omit(x)==na.omit(x)[1]))]
   if (debug) cat(paste0("The data have ", ncol(zc_y), " variables\n"))
  # J'enlève les colonnes qui ne sont pas numériques
  if (debug) cat("I remove the non-numeric columns\n")
  # Là c'est plutôt un inherits
  zc_y <- zc_y[, (sapply(zc_y, class) == "numeric") | (sapply(zc_y, class) == "integer")]
   if (debug) cat(paste0("The data have ", ncol(zc_y), " variables\n"))
  # J'enlève les colonnes qui empêchent de calculer les corrélations
  if (debug) cat("I remove the columns producing an error for correlations estimations\n")

  repeat {
    c <- suppressWarnings(expr=cor(zc_y, use=use, method=method))
    k <- apply(c, MARGIN=1, FUN=function(x) sum(ifelse(is.na(x), 1, 0)))
    # names(k) <- colnames(zc_y)
    # km1 <- max(k)
    if (!is.null(variable.retain)) k[variable.retain] <- 0
    km <- max(k)
    if (km == 0) break
    if (debug) cat(paste0("Remove: ", colnames(zc_y)[which.max(k)], "\n"))
    zc_y <- zc_y[, -which.max(k)]
  }
  
  if (debug) cat(paste0("The data have ", ncol(zc_y), " variables\n"))
  
  if (test.partial.correlation) {
    if (debug) cat("I remove the columns producing an error for partial correlations estimations\n")
  # J'enlève les colonnes qui empêchent les calculs des corrélations partielles
  
    lapp <- lapply
  if (is.element("parallel", installed.packages()[,1]))  lapp <- getFromNamespace("mclapply", ns="parallel")
  if (progress & (is.element("pbmcapply", installed.packages()[,1])))  lapp <- getFromNamespace("pbmclapply", ns="pbmcapply")
  # if (debug) lapp <- lapply
  
  options(mc.cores = getFromNamespace("detectCores", ns="parallel")())
  
  cname <- colnames(zc_y)
  lc <- length(cname)
  llc <- 1:lc
  # df <- expand.grid(e=1:(lc-2), f=2:(lc-1), g=3:lc)
  # df <- df[(df[, 1]<df[, 2]) & (df[, 2]<df[, 3]), ]
  
  # lc <- 5
  
  df <- expand.grid(e=1:(lc-2), f=2:(lc-1), g=3:lc)
  df <- df[(df[, 1]<df[, 2]) & (df[, 2]<df[, 3]), ]
  
  # pc <- array(data=Inf, dim=c(lc, lc, lc), dimnames = list(cname, cname, cname))
  

  z <- lapp(1:nrow(df), FUN = function(x) {
    e <- df[x, 1]
    f <- df[x, 2]
    g <- df[x, 3]
    
    dfg <- na.omit(zc_y[, c(e, f, g)])
    pc <- try(suppressWarnings(pcor(dfg, method=method)), 
              silent=TRUE)
    outg <- (inherits(pc, "try-error"))
    if (!outg) outg <- any(!is.finite(pc$estimate))
    # if (outg==0) {
    #   outg <- c(outg, pc$estimate[1, 2], pc$estimate[1, 3], pc$estimate[2, 3])
    # } else {
    #   outg <- c(outg, NA, NA, NA)
    # }
      return(!outg)
    }
  )
  
  
  # Si TRUE, c'est bon
  # FALSE c'est une erreur
  zul <- unlist(z)
  
  asupprimer <- NULL
    repeat {

    if (all(zul)) break
    
    zdf <- as.data.frame(table(as.vector(as.matrix(df[!zul, ]))), stringsAsFactors = FALSE)
    
    if (!is.null(variable.retain)) {
      zdf <- cbind(zdf, name=colnames(zc_y)[as.numeric(zdf$Var1)])
      zdf[zdf$name == variable.retain, "Freq"] <- 0
    }
    as <- as.numeric(zdf$Var1[which.max(zdf$Freq)])
    
    if (debug) cat(paste0("Remove: ", colnames(zc_y)[as], "\n"))
    
    zul[which(df[,1]==as)] <- TRUE
    zul[which(df[,2]==as)] <- TRUE
    zul[which(df[,3]==as)] <- TRUE

    asupprimer <- c(asupprimer, as)
    }
  
  if (!is.null(asupprimer)) zc_y <- zc_y[, -asupprimer]
  if (debug) cat(paste0("The data have ", ncol(zc_y), " variables\n"))
  }
  
return(zc_y)
}
