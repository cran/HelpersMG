#' IC_threshold_matrix calculates correlation matrix thresholed by partial correlation
#' @title Calculate correlation matrix
#' @author Marc Girondot
#' @return A list
#' @param data A dataframe or an IconoCorel object from a previous run of IC_threshold_matrix
#' @param threshold threshold for partial and full correlations
#' @param use an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param method a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param model a character string indicating if linear model uses all variables at a time (AAT) or one at a time (OAT).
#' @param significance.level if FALSE, does not use significance level; or use this significance level.
#' @param correction.multiple.comparisons "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", or "none".
#' @param progress show a progress bar
#' @param debug display information about progression of computing
#' @description This function calculates the matrix of correlations thresholded using partial correlation.\cr
#' If the threshold is not given, the object that is produced can be used later for thresholding.\cr
#' For model OAT: The link between A and B is “remarkable” if and only if the total correlation between them is higher than a given threshold and if the partial correlation between A and B in respect to any other variable C is also higher in absolute values than this threshold and with the same sign as the total correlation.
#' For model AAT: A correlation is retained if it is higher than the threshold and the partial correlation is lower than the threshold. In this case, no missing value is accepted.\cr
#' The use and method parameters are used by cor() function. The function uses by default a parallel computing in Unix or MacOSX systems. 
#' If progress is TRUE and the package pbmcapply is present, a progress bar is displayed. If debug is TRUE, some informations are shown during the process but parallel computing is not used.\cr
#' \code{https://fr.wikipedia.org/wiki/Iconographie_des_corrélations}
#' @family Iconography of correlations
#' @references Lesty, M., 1999. Une nouvelle approche dans le choix des régresseurs de la régression multiple en présence d’interactions et de colinéarités. Revue de Modulad 22, 41-77.
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' es <- structure(list(Student = c("e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8"), 
#'                  Mass = c(52, 59, 55, 58, 66, 62, 63, 69), 
#'                  Age = c(12, 12.5, 13, 14.5, 15.5, 16, 17, 18), 
#'                  Assiduity = c(12, 9, 15, 5, 11, 15, 12, 9), 
#'                  Note = c(5, 5, 9, 5, 13.5, 18, 18, 18)), 
#'                  row.names = c(NA, -8L), class = "data.frame")
#' 
#' es
#' 
#' df_clean <- IC_clean_data(es, debug = TRUE)
#' cor_matrix <- IC_threshold_matrix(data=df_clean, threshold = NULL, progress=FALSE)
#' cor_threshold <- IC_threshold_matrix(data=df_clean, threshold = 0.3)
#' plot(cor_threshold, show.legend.strength=FALSE, show.legend.direction = FALSE)
#' cor_threshold_Note <- IC_correlation_simplify(matrix=cor_threshold, variable="Note")
#' plot(cor_threshold_Note)
#' 
#' cor_threshold <- IC_threshold_matrix(data=df_clean, threshold = 0.8, progress=FALSE)
#' gr <- plot(cor_threshold, plot=FALSE)
#' ly <- getFromNamespace("layout_nicely", ns="igraph")(gr)
#' plot(cor_threshold, 
#' layout=matrix(data=c(53, 53, 55, 55, 
#'                      55, 53, 55, 53), ncol=2, byrow=FALSE), 
#' show.legend.direction = FALSE,
#' show.legend.strength = FALSE, xlim=c(-2, 2), ylim=c(-2, 2))
#' 
#' # Using significance level
#' 
#' cor_threshold <- IC_threshold_matrix(data=df_clean, threshold = 0.3, 
#'                                      significance.level=0.05)
#' plot(cor_threshold, show.legend.strength=FALSE, show.legend.direction = FALSE)
#' cor_threshold_Note <- IC_correlation_simplify(matrix=cor_threshold, variable="Note")
#' plot(cor_threshold_Note)
#' 
#' # Using the model All at a time
#' 
#' cor_threshold_AAT <- IC_threshold_matrix(data=df_clean, threshold = 0.3, model="AAT")
#' par(mar=c(1,1,1,1))
#' set.seed(4)
#' plot(cor_threshold_AAT, show.legend.strength="bottomleft")
#' 
#' 
#' 
#' ############
#' dta <- structure(list(Student = c("e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8"), 
#'                      Mass = c(52, 59, 55, 58, 66, 62, 63, 69), 
#'                      Age = c(12, 12.5, 13, 14.5, 15.5, 16, 17, 18), 
#'                      Assiduity = c(12, 9, 15, 5, 11, 15, 12, 9), 
#'                      Note = c(5, 5, 9, 5, 13.5, 18, 18, 18)), 
#'                      row.names = c(NA, -8L), class = "data.frame")
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


IC_threshold_matrix <- function(data=stop("A dataframe or an IconoCorel object is required")     , 
                                threshold=NULL                                                   , 
                                use = c("pairwise.complete.obs", "everything", 
                                        "all.obs", "complete.obs", "na.or.complete")             , 
                                method=c("pearson", "kendall", "spearman")                       ,
                                model=c("OAT", "ATT")                                            ,
                                significance.level=FALSE                                         ,
                                correction.multiple.comparisons="fdr"                            ,
                                progress=TRUE                                                    , 
                                debug=FALSE                                                      ) {
  
  
  
  
  
  # print(class(data))
  if (class(data) == "IconoCorel") {
    cor_mat <- data$correlation
    p_cor_mat <- data$correlation.pvalue
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
    
    model <- match.arg(model, choices = c("OAT", "ATT"))
    method <- match.arg(method, choices = c("pearson", "kendall", "spearman"))
    use <- match.arg(use, choices = c("pairwise.complete.obs", "everything", 
                                      "all.obs", "complete.obs", "na.or.complete"))
    # Là j'ai la matrice des corrélations
    cor_mat <- cor(data, method=method, use=use)
    
    if (use == "pairwise.complete.obs") {
    
    n <- nrow(data)
    daf <- expand.grid(e=1:ncol(data), f=1:ncol(data))
    n <- sapply(1:nrow(daf), FUN=function(x) {
      sum(as.numeric((!is.na(data[, daf[x, 1]])) & (!is.na(data[, daf[x, 2]]))))
    })
    
    n_cor_mat <- cor_mat
    n_cor_mat[] <- NA
    for (x in 1:nrow(daf)) n_cor_mat[daf[x, 1], daf[x, 2]] <- n[x]
    
    statistic <- cor_mat * sqrt((n_cor_mat - 2)/(1 - cor_mat^2))
    p_cor_mat <- 2 * pt(-abs(statistic), (n_cor_mat - 2))
    
    daf <- expand.grid(e=1:(ncol(data)-1), f=2:ncol(data))
    daf <- daf[(daf[, 1]<daf[, 2]), ]
    
    ppp <- NULL
    for (x in 1:nrow(daf)) ppp <- c(ppp, p_cor_mat[daf[x, 1], daf[x, 2]])
    
    ppp <- p.adjust(ppp, method = correction.multiple.comparisons, 
             n=length(ppp))
    
    
    for (x in 1:nrow(daf)) {
      p_cor_mat[daf[x, 1], daf[x, 2]] <- ppp[x]
      p_cor_mat[daf[x, 2], daf[x, 1]] <- ppp[x]
    }
    } else {
      p_cor_mat <- NULL
    }
   
    
    if (model == "OAT") {
      
      cname <- colnames(data)
      lc <- length(cname)
      llc <- 1:lc
      
      lapp <- lapply
      if (is.element("parallel", installed.packages()[,1])) lapp <- getFromNamespace("mclapply", ns="parallel")
      if (progress & (is.element("pbmcapply", installed.packages()[,1]))) lapp <- getFromNamespace("pbmclapply", ns="pbmcapply")
      if (debug) lapp <- lapply
      
      options(mc.cores = getFromNamespace("detectCores", ns="parallel")())
      
      daf <- expand.grid(e=1:(lc-2), f=2:(lc-1), g=3:lc)
      daf <- daf[(daf[, 1]<daf[, 2]) & (daf[, 2]<daf[, 3]), ]
      
      
      
      zu <- lapp(1:nrow(daf), FUN = function(x) {
        
        e <- daf[x, 1]
        f <- daf[x, 2]
        g <- daf[x, 3]
        
        dfg <- na.omit(data[, c(e, f, g)])
        pc <- try(suppressWarnings(pcor(dfg, method=method)), 
                  silent=TRUE)
        outg <- (class(pc)=="try-error")
        if (!outg) outg <- any(!is.finite(pc$estimate))
        if (!outg) {
          return(c(1, pc$estimate[1, 2], pc$estimate[1, 3], pc$estimate[2, 3]))
        } else {
          return(0)
        }
      }
      )
      
      cor_partial <- array(data=NA, dim=c(lc, lc, lc), dimnames = list(cname, cname, cname))
      for(x in which(sapply(zu, FUN = function(z) z[1]==1))) {
        e <- daf[x, 1]
        f <- daf[x, 2]
        g <- daf[x, 3]
        cor_partial[e, f, g] <- cor_partial[f, e, g] <- zu[[x]][2]
        cor_partial[e, g, f] <- cor_partial[g, e, f] <- zu[[x]][3]
        cor_partial[f, g, e] <- cor_partial[g, f, e] <- zu[[x]][4]
      }
      
      
      # The link between A and B is “remarkable” if and only if the total correlation 
      # between them is higher than a given threshold and if the partial correlation 
      # between A and B in respect to any other variable C is also higher in absolute 
      # values than this threshold and with the same sign as the total correlation.
      cor_partial_2D <- cor_mat
      daf <- expand.grid(e=1:(lc-1), f=2:lc)
      daf <- daf[(daf[, 1]<daf[, 2]), ]
      for (x in 1:nrow(daf)) {
        cor_point <- cor_mat[daf[x, 1], daf[x, 2]]
        pcor_x <- na.omit(cor_partial[daf[x, 1], daf[x, 2], ])
        if (sign(cor_point) > 0) {
          pos <- which.min(pcor_x)
          
          cor_partial_2D[daf[x, 1], daf[x, 2]] <- pcor_x[pos]
          cor_partial_2D[daf[x, 2], daf[x, 1]] <- pcor_x[pos]
          
        } else {
          pos <- which.max(pcor_x)
          
          cor_partial_2D[daf[x, 1], daf[x, 2]] <- pcor_x[pos]
          cor_partial_2D[daf[x, 2], daf[x, 1]] <- pcor_x[pos]
          
        }
        
      }
      cor_seuil <- cor_partial_2D
      diag(cor_seuil) <- 0
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
    cor_mat_threshold <- cor_mat * as.numeric(cor_seuil_binary)
    if (!identical(significance.level, FALSE)) {
      if (is.null(p_cor_mat)) {
        message("Threshold by significance is not available for this use parameter.")
      } else {
      cor_seuil_binary_plevel <- (p_cor_mat <= significance.level)
      cor_seuil_binary <- cor_seuil_binary & cor_seuil_binary_plevel
      cor_mat_threshold <- cor_mat * as.numeric(cor_seuil_binary)
      }
    }
  } else {
    cor_seuil_binary <- NULL
    cor_mat_threshold <- NULL
  }
  
  
  
  out <- list(correlation=cor_mat, 
              correlation.pvalue=p_cor_mat, 
              thresholded_correlation=cor_mat_threshold, 
              thresholded_correlation_binary=cor_seuil_binary, 
              threshold_matrix=cor_seuil, 
              model=model, 
              use=use, 
              method=method, 
              threshold=threshold)
  class(out) <- "IconoCorel"
  
  return(out)
}
