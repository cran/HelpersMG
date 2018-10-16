#' list.packages lists the installed packages with their locations
#' @title List the installed packages with their locations
#' @author Marc Girondot
#' @return A list with the installed packages and their version.
#' @description List the installed packages with their locations and version.
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' list.packages()
#' }
#' @export


list.packages <- function() {
  ip <- as.data.frame(installed.packages(), stringsAsFactors = FALSE)
  lp <- levels(as.factor(ip[, "LibPath"]))
  grandL <- list()
  for (p in lp) {
    L <- subset(ip, subset=(ip$LibPath==p), select=c("Package", "Version"))
    rownames(L) <- L[, "Package"]
    grandL <- c(grandL, list(L))
  }
  names(grandL) <- lp
  return(grandL)
}
