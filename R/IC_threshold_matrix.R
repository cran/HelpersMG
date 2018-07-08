#' IC_threshold_matrix calculates correlation matrix thresholed by partial correlation
#' @title Calculate correlation matrix
#' @author Marc Girondot
#' @return A list
#' @param data A dataframe or an IconoCorel object from a previous run of IC_threshold_matrix
#' @param threshold threshold for partial and full correlations
#' @param use an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param method a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param model a character string indicating if linear model uses all variables at a time (AAT) or one at a time (OAT).
#' @param progress show a progress bar
#' @param debug display information about progression of computing
#' @description This function calculates the matrix of correlations thresholded using partial correlation.\cr
#' If the threshold is not given, the object that is produced can be used later for thresholding.\cr
#' For model OAT: a correlation is retained if it is higher that the threshold and if all partial correlationw of the two variables and any third one are all lower than the threshold.\cr
#' For model AAT: a correlation is retained if it is higher than the threshold and the partial correlation is lower than the threshold. In this case, no missing value is accepted.\cr
#' The use and method parameters are used by cor() function. The function uses by default a parallel computing in Unix or MacOSX systems. 
#' If progress is TRUE and the package pbmcapply is present, a progress bar is displayed. If debug is TRUE, some informations are shown during the process but parallel computing is not used.\cr
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
#' cor_threshold <- IC_threshold_matrix(data=cor_matrix, threshold = 0.3)
#' par(mar=c(1,1,1,1))
#' set.seed(4)
#' plot(cor_threshold)
#' cor_threshold_Note <- IC_correlation_simplify(matrix=cor_threshold, variable="Note")
#' plot(cor_threshold_Note)
#' 
#' # Using the model All at a time
#' 
#' cor_threshold_AAT <- IC_threshold_matrix(data=df, threshold = 0.3, model="AAT")
#' par(mar=c(1,1,1,1))
#' set.seed(4)
#' plot(cor_threshold_AAT, show.legend.strength="bottomleft")
#' 
#' 
#' 
#' ############
#' dta <- structure(list(Élève = structure(1:8, .Label = c("e1", "e2", 
#' "e3", "e4", "e5", "e6", "e7", "e8"), class = "factor"), Poids = c(52L, 
#' 59L, 55L, 58L, 66L, 62L, 63L, 69L), Âge = c(12, 12.5, 13, 14.5, 
#' 15.5, 16, 17, 18), Assiduité = c(12L, 9L, 15L, 5L, 11L, 15L, 
#' 12L, 9L), Note = c(5, 5, 9, 5, 13.5, 18, 18, 18), e1 = c(1L, 
#' 0L, 0L, 0L, 0L, 0L, 0L, 0L), e2 = c(0L, 1L, 0L, 0L, 0L, 0L, 0L, 
#' 0L), e3 = c(0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L), e4 = c(0L, 0L, 0L, 
#' 1L, 0L, 0L, 0L, 0L), e5 = c(0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L), 
#'     e6 = c(0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L), e7 = c(0L, 0L, 0L, 
#'     0L, 0L, 0L, 1L, 0L), e8 = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L
#'     )), .Names = c("Élève", "Poids", "Âge", "Assiduité", 
#' "Note", "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8"), class = "data.frame", row.names = c(NA, 
#' -8L))
#'
#' dta0 <- dta[, 2:ncol(dta)]
#' ic0 <- IC_threshold_matrix(data = dta0)
#' cor_threshold <- IC_threshold_matrix(data=ic0, threshold = 0.3)
#' par(mar=c(1,1,1,1))
#' set.seed(4)
#' library("igraph")
#' 
#' plot(cor_threshold, vertex.color="red", show.legend.strength = FALSE)
#' plot(IC_correlation_simplify(matrix=cor_threshold), 
#'      show.legend.strength = FALSE, show.legend.direction = FALSE)
#' 
#' }
#' @export


IC_threshold_matrix <- function(data=stop("A dataframe or an IconoCorel object is required"), 
                             threshold=NULL, 
                             use = c("pairwise.complete.obs", "everything", 
                                     "all.obs", "complete.obs", "na.or.complete"), 
                             method=c("pearson", "kendall", "spearman"),
                             model=c("OAT", "ATT"),
                             progress=TRUE, debug=FALSE) {
  
 # print(class(data))
  if (class(data) == "IconoCorel") {
    cor_mat <- data$correlation
    cor_mat_threshold <- data$thresholded_correlation 
    cor_seuil_binary <- data$thresholded_correlation_binary 
    cor_seuil <- data$threshold_matrix
    
  } else {
    
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
  model <- model[1]
  method <- method[1]
  use <- use[1]
  cor_mat <- cor(data, method=method, use=use)
  
  if (model == "OAT") {
  
  cor_seuil <- cor_mat
  cor_seuil[] <- 0
  

  
  cname <- colnames(data)
  lc <- length(cname)
  llc <- 1:lc
  df <- expand.grid(e=1:(lc-1), f=2:lc)
  df <- df[(df[, 1]<df[, 2]), ]
  
  lapp <- lapply
  if (is.element("parallel", installed.packages()[,1])) lapp <- getFromNamespace("mclapply", ns="parallel")
  if (progress & (is.element("pbmcapply", installed.packages()[,1]))) lapp <- getFromNamespace("pbmclapply", ns="pbmcapply")
  if (debug) lapp <- lapply
  
    options(mc.cores = getFromNamespace("detectCores", ns="parallel")())
  
  z <- lapp(1:nrow(df), FUN = function(x) {
    e <- df[x, 1]
    f <- df[x, 2]
      if (debug) print(paste("index:", x, "     e:", e, "  f: ", f))
      
      k <- 1
      signe <- TRUE
      
      for (g in llc[c(-e, -f)]) {
        
        pc <- pcor(na.omit(data[, c(e, f, g)]), method=method)$estimate[1, 2]
        if (debug & is.na(pc)) print(paste("index:", x, "     e:", e, "  f: ", f, "  g:", g, "  partial correlation=", pc))
         if(abs(pc)<abs(k)) k <- abs(pc)
          if (sign(pc) != sign(cor_mat[e, f])) signe <- FALSE
      }
      if (!signe) k <- 0
      return(k)
  } 
  )
    z <- unlist(z)
    e <- df[, 1]
    f <- df[, 2]
    for (i in seq_along(z)) {
      cor_seuil[e[i], f[i]] <- z[i]
      cor_seuil[f[i], e[i]] <- z[i]
    }
    
  } else {
    
    cor_seuil <- pcor(data, method=method)$estimate
    diag(cor_seuil) <- 0
  }

  }

  if (!is.null(threshold)) {
    # Pour éviter les redondances, le lien AB est tracé si et seulement si la 
    # corrélation totale r(A,B) est supérieure au seuil en valeur absolue, 
    # et si les corrélations partielles r(A,B), par rapport à une variable Z, 
    # sont supérieures au seuil, en valeur absolue, et de même signe que la 
    # corrélation totale, pour tout Z parmi les variables disponibles, y 
    # compris les « instants ».
    # cor_seuil_binary <- (abs(cor_seuil)>threshold) & (abs(cor_mat)>threshold) & (sign(cor_mat) == sign(cor_seuil))
    cor_seuil_binary <- (abs(cor_seuil)>threshold) & (abs(cor_mat)>threshold)
    cor_mat_threshold <- cor_mat * ifelse(cor_seuil_binary, 1, 0)
  } else {
    cor_seuil_binary <- NULL
    cor_mat_threshold <- NULL
  }
  

  
  out <- list(correlation=cor_mat, thresholded_correlation=cor_mat_threshold, 
              thresholded_correlation_binary=cor_seuil_binary, 
              threshold_matrix=cor_seuil, 
              model=model, 
              use=use, 
              method=method, 
              threshold=threshold)
  class(out) <- "IconoCorel"
  
  return(out)
}