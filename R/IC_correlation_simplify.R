#' IC_correlation_simplify simplifies the correlation matrix
#' @title Simplify the correlation matrix
#' @author Marc Girondot
#' @return A list
#' @param matrix The correlation matrix to simplify
#' @param variable a vector with the name of columns to keep
#' @description This function can be used to simplify the network of correlations.\cr
#' If no vector of variables is given, the variables not linked to any other variable are removed. 
#' If a vector of variables is given, only link to these variables are retained.
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


IC_correlation_simplify <- function(matrix, variable=NULL) {
  if (class(matrix) != "IconoCorel") {
    stop("Only an object obtained using IC_threshold_matrix() can be used.")
  }
  m <- matrix$thresholded_correlation
  if (is.null(variable)) {
    garde <- !apply(X = m, MARGIN = 1, FUN=function(x) all(abs(x)<1E-5))
  } else {
    garde <- m[, variable, drop=FALSE] != 0
    garde[colnames(x = m) %in% variable] <- TRUE
    garde <- apply(garde, MARGIN = 1, FUN=any)
  }
  matrix$thresholded_correlation <- m[garde, garde]
  return(matrix)
}
