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


IC_correlation_simplify <- function(matrix, variable=NULL) {
  m <- matrix$thresholded_correlation
  if (is.null(variable)) {
    garde <- !apply(X = m, MARGIN = 1, FUN=function(x) all(abs(x)<1E-5))
  } else {
    garde <- m[, variable] != 0
    garde[colnames(x = m) %in% variable] <- TRUE
  }
  matrix$thresholded_correlation <- m[garde, garde]
  return(matrix)
}
