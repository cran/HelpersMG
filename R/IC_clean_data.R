#' IC_clean_data checks and corrects the dataframe to be used with IC_threshold_matrix
#' @title Clean the dataframe before to be used with IC_threshold_matrix
#' @author Marc Girondot
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
#' es <- matrix(c("e1", "52", "12", "12", "5",
#' "e2", "59", "12.5", "9", "5",
#' "e3", "55", "13", "15", "9",
#' "e4", "58", "14.5", "5", "5",
#' "e5", "66", "15.5", "11", "13.5",
#' "e6", "62", "16", "15", "18",
#' "e7", "63", "17", "12", "18",
#' "e8", "69", "18", "9", "18"), ncol=5, byrow = TRUE)
#' colnames(es) <- c("Élève", "Poids", "Âge", "Assiduité", "Note")
#' es <- as.data.frame(es, stringsasFactor=FALSE)
#' es[, 2] <- as.numeric(as.character(es[, 2]))
#' es[, 3] <- as.numeric(as.character(es[, 3]))
#' es[, 4] <- as.numeric(as.character(es[, 4]))
#' es[, 5] <- as.numeric(as.character(es[, 5]))
#' 
#' es
#' 
#' df <- IC_clean_data(es, debug = TRUE)
#' cor_matrix <- IC_threshold_matrix(data=df, threshold = NULL, progress=FALSE)
#' cor_threshold <- IC_threshold_matrix(data=df, threshold = 0.3)
#' par(mar=c(1,1,1,1))
#' set.seed(4)
#' plot(cor_threshold)
#' cor_threshold_Note <- IC_correlation_simplify(matrix=cor_threshold, variable="Note")
#' plot(cor_threshold_Note)
#' }
#' @export


IC_clean_data <- function(data=stop("A dataframe object is required"), 
                       use = c("pairwise.complete.obs", "everything", 
                                     "all.obs", "complete.obs", "na.or.complete"), 
                       method=c("pearson", "kendall", "spearman"), 
                       variable.retain=NULL,
                       test.partial.correlation=TRUE, 
                       progress=TRUE, debug=FALSE) {
  
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
  zc_y <- zc_y[, !apply(X = zc_y, MARGIN = 2, FUN = function(x) all(is.na(x)))]
   if (debug) cat(paste0("The data have ", ncol(zc_y), " variables\n"))
  # J'enlève les colonnes qui ne présente pas de variabilité
  if (debug) cat("I remove the columns with no variability\n")
  zc_y <- zc_y[, !apply(X = zc_y, MARGIN = 2, FUN = function(x) all(na.omit(x)==na.omit(x)[1]))]
   if (debug) cat(paste0("The data have ", ncol(zc_y), " variables\n"))
  # J'enlève les colonnes qui ne sont pas numériques
  if (debug) cat("I remove the non-numeric columns\n")
  zc_y <- zc_y[, (sapply(zc_y, class) == "numeric") | (sapply(zc_y, class) == "integer")]
   if (debug) cat(paste0("The data have ", ncol(zc_y), " variables\n"))
  # J'enlève les colonnes qui empêchent de calculer les corrélations
  if (debug) cat("I remove the columns producing an error for correlations estimations\n")

  repeat {
    c <- suppressWarnings(expr=cor(zc_y, use=use, method=method))
    k <- apply(c, MARGIN=1, FUN=function(x) sum(ifelse(is.na(x), 1, 0)))
    # names(k) <- colnames(zc_y)
    km1 <- max(k)
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
  df <- expand.grid(e=1:(lc-1), f=2:(lc))
  df <- df[(df[, 1]<df[, 2]), ]

  z <- lapp(1:nrow(df), FUN = function(x) {
    e <- df[x, 1]
    f <- df[x, 2]
      outg <- NULL
      for (g in llc[c(-e, -f)]) {
        dfg <- na.omit(zc_y[, c(e, f, g)])
        pc <- try(suppressWarnings(pcor(dfg, method=method)), 
                  silent=TRUE)
        if (class(pc)=="try-error") {
          outg <- c(outg, e, f, g)
        } else {
          if (is.na(pc$estimate[1, 2])) {
            outg <- c(outg, e, f, g)
          }
        }
      }
      return(outg)
    }
  )
  asupprimer <- NULL
    repeat {
    zul <- unlist(z)
    
    if (identical(zul, integer(0)) | is.null(zul)) break
    
    zdf <- as.data.frame(table(zul), stringsAsFactors = FALSE)
    zdf <- cbind(zdf, name=colnames(zc_y)[as.numeric(zdf$z)])
    if (!is.null(variable.retain)) zdf[zdf$name == variable.retain, "Freq"] <- 0
    as <- as.numeric(zdf$z[which.max(zdf$Freq)])
    
    if (debug) cat(paste0("Remove: ", colnames(zc_y)[as], "\n"))
    
    z[which(df[,1]==as)] <- NULL
    z[which(df[,2]==as)] <- NULL
    
    zbis <- lapply(z, FUN=function(x) x[x != as])
    
    z <- zbis
    
    asupprimer <- c(asupprimer, as)
    }
  
  if (!is.null(asupprimer)) zc_y <- zc_y[, -asupprimer]
  if (debug) cat(paste0("The data have ", ncol(zc_y), " variables\n"))
  }
  
return(zc_y)
}
